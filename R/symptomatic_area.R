#'Calculates the percentage of symptomatic area
#'
#'Calculates the percentage of symptomatic leaf area in sample or entire
#'leaf based on provided color palettes samples. A general linear model
#'(binomial family) fitted to the RGB values is used to segment the lesions from
#'the healthy leaf. If a pallet of background is provided, the function takes
#'care of the details to isolate it before computing the number and area of
#'lesions. By using `img_pattern` it is possible to process several images with
#'common pattern names that are stored in the current working directory or in
#'the subdirectory informed in `dir_original`.
#' @param img The image to be analyzed.
#' @param img_pattern A pattern of file name used to identify images to be
#'   processed. For example, if `img_pattern = "im"` all images that the name
#'   matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
#'   Providing any number as pattern (e.g., `img_pattern = "1"`) will select
#'   images that are named as 1.-, 2.-, and so on.
#' @param img_healthy A color palette of healthy areas.
#' @param img_symptoms A color palette of symptomatic areas.
#' @param img_background A color palette of areas with symptoms.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time, especially when `img_pattern` is used is informed. The
#'   number of sections is set up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing?
#' @param show_original Show the symptoms in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param col_leaf Leaf color after image processing. Defaults to `"green"`
#' @param col_symptoms Symptoms color after image processing. Defaults to `"red"`.
#' @param col_background Background color after image processing. Defaults to
#'   `"NULL"`.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named with the prefix provided in `proc_*` where `*` is the image name
#'   given in `img`.
#' @param prefix The prefix to be included in the processed images. Defaults to
#'   `"proc_"`.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return A data frame with the results (healthy and symptomatic area) for each
#'   image.
#' @export
#' @md
#' @importFrom stats binomial glm predict kmeans sd aggregate
#' @importFrom grid grid.raster
#' @importFrom grDevices col2rgb dev.off png hcl.colors
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(image_pliman("sev_leaf.jpg"))
#' healthy <- image_import(image_pliman("sev_healthy.jpg"))
#' symptoms <- image_import(image_pliman("sev_sympt.jpg"))
#' background <- image_import(image_pliman("sev_back.jpg"))
#' image_combine(img, healthy, symptoms,background)
#' symptomatic_area(img = img,
#'                  img_healthy = healthy,
#'                  img_symptoms = symptoms,
#'                  img_background = background,
#'                  show_image = TRUE)
#' }
#'
#'
symptomatic_area <- function(img,
                             img_healthy,
                             img_symptoms,
                             img_background = NULL,
                             img_pattern = NULL,
                             parallel = FALSE,
                             workers = NULL,
                             randomize = TRUE,
                             nrows = 10000,
                             show_image = FALSE,
                             show_original = TRUE,
                             show_background = TRUE,
                             col_leaf = "green",
                             col_symptoms = "red",
                             col_background = NULL,
                             save_image = FALSE,
                             prefix = "proc_",
                             dir_original = NULL,
                             dir_processed = NULL,
                             verbose = TRUE){
  if(!missing(img) & !missing(img_pattern)){
    stop("Only one of `img` or `img_pattern` arguments can be used.", call. = FALSE)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <- paste("./", dir_original, sep = "")
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste("./", sep = "")
  } else{
    diretorio_processada <- paste("./", dir_processed, sep = "")
  }
  help_sympt <-
    function(img, img_healthy, img_symptoms, img_background, randomize,
             nrows, show_image, show_original, show_background, col_background,
             save_image, dir_original, dir_processed){
      # Some parts adapted from
      # https://github.com/AlcineiAzevedo/Segmentacao-conchonilha2
      # Thanks to Alcinei Azevedo for his tips
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img)
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "jpg"
      }
      if(is.character(img_healthy)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        imag <- list.files(diretorio_original, pattern = img_healthy)
        check_names_dir(img_healthy, all_files, diretorio_original)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_healthy <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      if(is.character(img_symptoms)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        imag <- list.files(diretorio_original, pattern = img_symptoms)
        check_names_dir(img_symptoms, all_files, diretorio_original)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_symptoms <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      original <- image_to_mat(img)
      sadio <- image_to_mat(img_healthy)
      sintoma <- image_to_mat(img_symptoms)
      ################## no background #############
      if(is.null(img_background)){
        sadio_sintoma <-
          transform(rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                          sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]),
                    Y = ifelse(CODE == "img_healthy", 1, 0))
        usef_area <- nrow(original$df_in)
        model <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
        # isolate plant
        pred1 <- round(predict(model, newdata = original$df_in, type="response"), 0)
        plant_symp <- matrix(pred1, ncol = ncol(original$R))
        plant_symp <- image_correct(plant_symp, perc = 0.01)
        ID <- c(plant_symp == 0)
        pix_sympt <- length(which(ID == TRUE))
        if(show_original == TRUE){
          if(is.null(col_background)){
            col_background <- col2rgb("green")
          } else{
            col_background <- col2rgb(col_background)
          }
          col_symptoms <- col2rgb(col_symptoms)
          im2 <- img
          im2@.Data[,,1][ID] <- col_symptoms[1]
          im2@.Data[,,2][ID] <- col_symptoms[2]
          im2@.Data[,,3][ID] <- col_symptoms[3]
          if(show_background == FALSE){
            im2@.Data[,,1][!ID] <- col_background[1]
            im2@.Data[,,2][!ID] <- col_background[2]
            im2@.Data[,,3][!ID] <- col_background[3]
          }
        } else{
          if(is.null(col_background)){
            col_background <- col2rgb("green")
          } else{
            col_background <- col2rgb(col_background)
          }
          col_symptoms <- col2rgb(col_symptoms)
          im2 <- img
          im2@.Data[,,1][ID] <- col_symptoms[1]
          im2@.Data[,,2][ID] <- col_symptoms[2]
          im2@.Data[,,3][ID] <- col_symptoms[3]
          im2@.Data[,,1][!ID] <- col_background[1]
          im2@.Data[,,2][!ID] <- col_background[2]
          im2@.Data[,,3][!ID] <- col_background[3]
        }
      } else{
        if(is.character(img_background)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = img_background)
          check_names_dir(img_background, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          img_background <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        fundo <- image_to_mat(img_background)
        # separate image from background
        fundo_resto <-
          transform(rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                          sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,],
                          fundo$df_in[sample(1:nrow(fundo$df_in)),][1:nrows,]),
                    Y = ifelse(CODE == "img_background", 0, 1))
        modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = fundo_resto))
        pred1 <- round(predict(modelo1, newdata = original$df_in, type="response"), 0)
        plant_background <- matrix(pred1, ncol = ncol(original$R))
        plant_background <- image_correct(plant_background, perc = 0.009)
        plant_background[plant_background == 1] <- 2
        sadio_sintoma <-
          transform(rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                          sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]),
                    Y = ifelse(CODE == "img_healthy", 1, 0))
        modelo2 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
        # isolate plant
        ID <- c(plant_background == 2)
        usef_area <- nrow(original$df_in[ID,])
        pred3 <- round(predict(modelo2, newdata = original$df_in[ID,], type="response"), 0)
        pix_sympt <- length(which(pred3 == 0))
        if(show_original == TRUE){
          im2 <- img
          col_symptoms <- col2rgb(col_symptoms)
          im2@.Data[,,1][ID][which(pred3 == 0)] <- col_symptoms[1]
          im2@.Data[,,2][ID][which(pred3 == 0)] <- col_symptoms[2]
          im2@.Data[,,3][ID][which(pred3 == 0)] <- col_symptoms[3]
          if(!is.null(col_background)){
            col_background <- col2rgb(col_background)
            im2@.Data[,,1][!ID] <- col_background[1]
            im2@.Data[,,2][!ID] <- col_background[2]
            im2@.Data[,,3][!ID] <- col_background[3]
          }
        } else{
          if(is.null(col_background)){
            col_background <- col2rgb("white")
          } else{
            col_background <- col2rgb(col_background)
          }
          col_leaf <- col2rgb(col_leaf)
          col_symptoms <- col2rgb(col_symptoms)
          im2 <- img
          im2@.Data[,,1][ID][which(pred3 == 0)] <- col_symptoms[1]
          im2@.Data[,,2][ID][which(pred3 == 0)] <- col_symptoms[2]
          im2@.Data[,,3][ID][which(pred3 == 0)] <- col_symptoms[3]
          im2@.Data[,,1][ID][which(pred3 != 0)] <- col_leaf[1]
          im2@.Data[,,2][ID][which(pred3 != 0)] <- col_leaf[2]
          im2@.Data[,,3][ID][which(pred3 != 0)] <- col_leaf[3]
          im2@.Data[,,1][!ID] <- col_background[1]
          im2@.Data[,,2][!ID] <- col_background[2]
          im2@.Data[,,3][!ID] <- col_background[3]
        }
      }
      if(show_image == TRUE){
        plot(im2)
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada)
        }
        image_export(im2,
                   name = paste0(diretorio_processada, "/",
                                 prefix,
                                 name_ori, ".",
                                 extens_ori))
      }
      symptomatic <- pix_sympt /  usef_area * 100
      healthy <- 100 - symptomatic
      results <- data.frame(healthy = healthy,
                            symptomatic = symptomatic)
      return(results)
    }
  if(missing(img_pattern)){
    help_sympt(img, img_healthy, img_symptoms, img_background, randomize,
               nrows, show_image, show_original, show_background, col_background,
               save_image, dir_original, dir_processed)
  } else{
    if(img_pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      img_pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = img_pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(img_pattern, names_plant)) == 0){
      stop(paste("'", img_pattern, "' pattern not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_sympt", "file_name",
                                "check_names_dir", "file_extension", "image_import",
                                "image_binary", "watershed", "distmap", "computeFeatures.moment",
                                "computeFeatures.shape", "colorLabels", "image_show",
                                "image_to_mat", "image_correct", "image_export"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }

      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_sympt(x,
                               img_healthy, img_symptoms, img_background, randomize,
                               nrows, show_image, show_original, show_background, col_background,
                               save_image, dir_original, dir_processed)
                  })

    } else{
    results <- list()
    pb <- progress(max = length(plants), style = 4)
    for (i in 1:length(plants)) {
      run_progress(pb, actual = i,
                   text = paste("Processing image", names_plant[i]))
      results[[i]] <-
        help_sympt(img  = names_plant[i],
                   img_healthy, img_symptoms, img_background, randomize,
                   nrows, show_image, show_original, show_background, col_background,
                   save_image, dir_original, dir_processed)
    }
    }
    results <- transform(do.call(rbind, results),
                         sample = names_plant)[c(3, 1, 2)]
    return(results)
  }
}
