#'Calculates the percentage of symptomatic area
#'
#'Calculates the percentage of symptomatic leaf area in sample or entire
#'leaf based on provided color palettes samples. A general linear model
#'(binomial family) fitted to the RGB values is used to segment the lesions from
#'the healthy leaf. If a pallet of background is provided, the function takes
#'care of the details to isolate it before computing the number and area of
#'lesions. By using `pattern` it is possible to process several images with
#'common pattern names that are stored in the current working directory or in
#'the subdirectory informed in `dir_original`.
#' @param img The image to be analyzed.
#' @param img_healthy A color palette of healthy areas.
#' @param img_symptoms A color palette of symptomatic areas.
#' @param img_background A color palette of areas with symptoms.
#' @param pattern A pattern of file name used to identify images to be imported.
#'   For example, if `pattern = "im"` all images in the current working
#'   directory that the name matches the pattern (e.g., img1.-, image1.-, im2.-)
#'   will be imported as a list. Providing any number as pattern (e.g., `pattern
#'   = "1"`) will select images that are named as 1.-, 2.-, and so on. An error
#'   will be returned if the pattern matches any file that is not supported
#'   (e.g., img1.pdf).
#' @param img_pattern Deprecated. Use `pattern` instead.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param fill_hull Fill holes in the image? Defaults to `TRUE`. This is useful
#'   to fill holes in leaves, e.g., those caused by insect attack, ensuring the
#'   hole area will be accounted for the leaf, not background.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time, especially when `img_pattern` is used is informed. The
#'   number of sections is set up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time. Defaults
#'   to 50% of available cores.
#' @param nrows The number of lines to be used in training step. Defaults to
#'   3000.
#' @param show_image Show image after processing?
#' @param show_original Show the symptoms in the original image? Defaults to
#'   `TRUE`.
#' @param show_contour Show a contour line around the symptomatic area? Defaults
#'   to `TRUE`. If false, the symptomatic area will be filled with the color
#'   informed in `col_symptoms` argument.
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param col_leaf Leaf color after image processing. Defaults to `"green"`
#' @param col_symptoms Symptoms color after image processing. Defaults to `"red"`.
#' @param col_background Background color after image processing. Defaults to
#'   `"NULL"`.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named with the prefix provided in `proc_*` where
#'   `*` is the image name given in `img`.
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
#' @importFrom graphics lines
#' @importFrom grDevices col2rgb dev.off png hcl.colors jpeg
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' healthy <- image_pliman("sev_healthy.jpg")
#' symptoms <- image_pliman("sev_sympt.jpg")
#' background <- image_pliman("sev_back.jpg")
#' image_combine(img, healthy, symptoms,background)
#' symptomatic_area(img = img,
#'                  img_healthy = healthy,
#'                  img_symptoms = symptoms,
#'                  img_background = background,
#'                  show_image = TRUE)
#' }
symptomatic_area <- function(img,
                             img_healthy,
                             img_symptoms,
                             img_background = NULL,
                             pattern = NULL,
                             img_pattern = NULL,
                             resize = FALSE,
                             fill_hull = TRUE,
                             parallel = FALSE,
                             workers = NULL,
                             nrows = 3000,
                             show_image = FALSE,
                             show_original = TRUE,
                             show_contour = TRUE,
                             show_background = TRUE,
                             col_leaf = "green",
                             col_symptoms = "red",
                             col_background = NULL,
                             save_image = FALSE,
                             prefix = "proc_",
                             dir_original = NULL,
                             dir_processed = NULL,
                             verbose = TRUE){
  # check_ebi()
  if(!missing(img_pattern)){
    warning("Argument 'img_pattern' is deprecated. Use 'pattern' instead.", call. = FALSE)
    pattern <- img_pattern
  }
  if(!missing(img) & !missing(pattern)){
    stop("Only one of `img` or `pattern` arguments can be used.", call. = FALSE)
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
    function(img, img_healthy, img_symptoms, img_background,
             nrows, show_image, show_original, show_background, col_background,
             save_image, dir_original, dir_processed){
      # Some parts adapted from
      # https://github.com/AlcineiAzevedo/Segmentacao-conchonilha2
      # Thanks to Alcinei Azevedo for his tips
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
        if(resize != FALSE){
          img <- image_resize(img, resize)
        }
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
        if(resize != FALSE){
          img_healthy <- image_resize(img_healthy, resize)
        }
      }
      if(is.character(img_symptoms)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        imag <- list.files(diretorio_original, pattern = img_symptoms)
        check_names_dir(img_symptoms, all_files, diretorio_original)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_symptoms <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        if(resize != FALSE){
          img_symptoms <- image_resize(img_symptoms, resize)
        }
      }
      original <-
        data.frame(CODE = "img",
                   R = c(img@.Data[,,1]),
                   G = c(img@.Data[,,2]),
                   B = c(img@.Data[,,3]))
      sadio <-
        data.frame(CODE = "img_healthy",
                   R = c(img_healthy@.Data[,,1]),
                   G = c(img_healthy@.Data[,,2]),
                   B = c(img_healthy@.Data[,,3]))
      sintoma <-
        data.frame(CODE = "img_symptoms",
                   R = c(img_symptoms@.Data[,,1]),
                   G = c(img_symptoms@.Data[,,2]),
                   B = c(img_symptoms@.Data[,,3]))
      ncol_img <- dim(img)[[2]]

      ################## no background #############
      if(is.null(img_background)){
        sadio_sintoma <-
          transform(rbind(sadio[sample(1:nrow(sadio)),][1:nrows,],
                          sintoma[sample(1:nrow(sintoma)),][1:nrows,]),
                    Y = ifelse(CODE == "img_healthy", 1, 0))
        usef_area <- nrow(original)
        model <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
        # isolate plant
        pred1 <- round(predict(model, newdata = original, type = "response"), 0)
        plant_symp <- matrix(pred1, ncol = ncol_img)
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
          if(resize != FALSE){
            img_background <- image_resize(img_background, resize)
          }
        }
        fundo <-
          data.frame(CODE = "img_background",
                     R = c(img_background@.Data[,,1]),
                     G = c(img_background@.Data[,,2]),
                     B = c(img_background@.Data[,,3]))
        # separate image from background
        fundo_resto <-
          transform(rbind(sadio[sample(1:nrow(sadio)),][1:nrows,],
                          sintoma[sample(1:nrow(sintoma)),][1:nrows,],
                          fundo[sample(1:nrow(fundo)),][1:nrows,]),
                    Y = ifelse(CODE == "img_background", 0, 1))
        modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = fundo_resto))
        pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
        ifelse(fill_hull == TRUE,
               plant_background <- EBImage::fillHull(matrix(pred1, ncol = ncol_img)),
               plant_background <- matrix(pred1, ncol = ncol_img))
        plant_background[plant_background == 1] <- 2
        sadio_sintoma <-
          transform(rbind(sadio[sample(1:nrow(sadio)),][1:nrows,],
                          sintoma[sample(1:nrow(sintoma)),][1:nrows,]),
                    Y = ifelse(CODE == "img_healthy", 1, 0))
        modelo2 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
        # isolate plant
        ID <- c(plant_background == 2)
        usef_area <- nrow(original[ID,])
        pred3 <- round(predict(modelo2, newdata = original[ID,], type="response"), 0)
        pix_sympt <- length(which(pred3 == 0))
        if(show_image == TRUE | save_image == TRUE){
          if(show_original == TRUE){
            im2 <- img
            if(isFALSE(show_contour)){
              col_symptoms <- col2rgb(col_symptoms)
              im2@.Data[,,1][ID][which(pred3 == 0)] <- col_symptoms[1]
              im2@.Data[,,2][ID][which(pred3 == 0)] <- col_symptoms[2]
              im2@.Data[,,3][ID][which(pred3 == 0)] <- col_symptoms[3]
            }
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
      }
      if(show_image == TRUE){
        if(show_contour){
          img_contour <- img@.Data[,,1]
          img_contour[ID][which(pred3 == 0)] <- 1
          img_contour[ID][which(pred3 == 1)] <- 0
          img_contour[!ID] <- 0
          countor_points <-
            EBImage::ocontour(
              EBImage::Image(
                EBImage::bwlabel(
                  EBImage::distmap(img_contour)
                )
              )
            )
        }
        plot(im2)
        plot_contour(countor_points)
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada, recursive = TRUE)
        }
        if(show_contour){
          img_contour <- img@.Data[,,1]
          img_contour[ID][which(pred3 == 0)] <- 1
          img_contour[ID][which(pred3 == 1)] <- 0
          img_contour[!ID] <- 0
          countor_points <-
            EBImage::ocontour(
              EBImage::Image(
                EBImage::bwlabel(
                  EBImage::distmap(img_contour)
                )
              )
            )
        }
        jpeg(paste0(diretorio_processada, "/",
                    prefix,
                    name_ori, ".",
                    extens_ori),
             width = dim(im2@.Data)[1],
             height = dim(im2@.Data)[2])
        plot(im2)
        plot_contour(countor_points)
        dev.off()
      }
      symptomatic <- pix_sympt /  usef_area * 100
      healthy <- 100 - symptomatic
      results <- data.frame(healthy = healthy,
                            symptomatic = symptomatic)
      return(results)
    }
  if(missing(pattern)){
    help_sympt(img, img_healthy, img_symptoms, img_background,
               nrows, show_image, show_original, show_background, col_background,
               save_image, dir_original, dir_processed)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("'", pattern, "' pattern not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Extensions that match the name pattern: ", paste(unique(extensions), collapse  = ", "),
           "\nAllowed extensions are .png, .jpeg, .jpg, .tiff", call. = FALSE)
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_sympt"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_sympt(x,
                               img_healthy, img_symptoms, img_background,
                               nrows, show_image, show_original, show_background,
                               col_background, save_image, dir_original,
                               dir_processed)
                  })

    } else{
      pb <- progress(max = length(plants), style = 4)
      foo <- function(plants, ...){
        run_progress(pb, ...)
        help_sympt(img  = plants,
                   img_healthy, img_symptoms, img_background,
                   nrows, show_image, show_original, show_background, col_background,
                   save_image, dir_original, dir_processed)
      }
      results <-
        lapply(seq_along(names_plant), function(i){
          foo(names_plant[i],
              actual = i,
              text = paste("Processing image", names_plant[i]))
        })
    }
    results <- transform(do.call(rbind, results),
                         sample = names_plant)[c(3, 1, 2)]
    return(results)
  }
}
