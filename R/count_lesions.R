#' Counts the number of lesions
#'
#'Counts the number of lesions in a sample or entire leaf based on provided
#'color palettes samples. A general linear model (binomial family) fitted to the
#'RGB values is used to segment the lesions from the healthy leaf. If a pallet
#'of background is provided, the function takes care of the details to isolate
#'it before computing the number and area of lesions. By using `img_pattern` it
#'is possible to process several images with common pattern names that are
#'stored in the current working directory or in the subdirectory informed in
#'`dir_original`.
#' @param img The image to be analyzed.
#' @param img_healthy A color palette of healthy areas.
#' @param img_lesion A color palette of lesioned areas.
#' @param img_background An optional color palette of the image background.
#' @param img_pattern A pattern of file name used to identify images to be
#'   processed. For example, if `img_pattern = "im"` all images that the name
#'   matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
#'   Providing any number as pattern (e.g., `img_pattern = "1"`) will select
#'   images that are nammed as 1.-, 2.-, and so on.
#' @param lower_size Lower limit for size for the image analysis. Leaf images
#'   often contain dirt and dust. To prevent dust from affecting the image
#'   analysis, the lower limit of analyzed size is set to 0.1, i.e., objects
#'   with lesser than 10% of the mean of all objects are removed. One can set a
#'   known area or use `lower_limit = 0` to select all objects (not advised).
#' @param upper_size Upper limit for size for the image analysis. Defaults to
#'   `NULL`, i.e., no upper limit used.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing?
#' @param show_original Show the symptoms in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param col_leaf Leaf color after image processing. Defaults to `"green"`
#' @param col_lesions Symptoms color after image processing. Defaults to `"red"`.
#' @param col_background Background color after image processing. Defaults to
#'   `"NULL"`.
#' @param text_col,text_size,text_digits The color, size and significant digits
#'   used in the text. The shows the pattern `o|a`, where `o` and `a` are the
#'   object id and its area, respectively.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named as `proc_*` where `*` is the image name
#'   given in `img`.
#' @param prefix The prefix to be included in the processed images. Defaults to
#'   `"proc_"`.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory.
#' @return A data frame with the results for each image.
#' @export
#' @md
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(system.file("tmp_images/sev3.png", package = "pliman"))
#' healthy <- image_import(system.file("tmp_images/sev_healthy.png", package = "pliman"))
#' lesions <- image_import(system.file("tmp_images/sev_sympt.png", package = "pliman"))
#' image_show(img)
#' image_show(healthy)
#' image_show(lesions)
#' count_lesions(img = img,
#'               img_healthy = healthy,
#'               img_lesion = lesions,
#'               show_image = TRUE)
#' }
#'
count_lesions <- function(img,
                          img_healthy,
                          img_lesion,
                          img_background = NULL,
                          img_pattern = NULL,
                          lower_size = NULL,
                          upper_size = NULL,
                          randomize = TRUE,
                          nrows = 10000,
                          show_image = FALSE,
                          show_original = TRUE,
                          show_background = TRUE,
                          col_leaf = "green",
                          col_lesions = "red",
                          col_background = NULL,
                          text_col = "black",
                          text_size = 1,
                          text_digits = 2,
                          save_image = FALSE,
                          prefix = "proc_",
                          dir_original = NULL,
                          dir_processed = NULL){
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
  help_count <-
    function(img, img_healthy, img_lesion, img_background, randomize,
             nrows, show_image, show_original, show_background, col_background,
             save_image, dir_original, dir_processed){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img)
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      if(is.character(img_healthy)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        imag <- list.files(diretorio_original, pattern = img_healthy)
        check_names_dir(img_healthy, all_files, diretorio_original)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_healthy <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      if(is.character(img_lesion)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        imag <- list.files(diretorio_original, pattern = img_lesion)
        check_names_dir(img_lesion, all_files, diretorio_original)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_lesion <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      original <- image_to_mat(img)
      sadio <- image_to_mat(img_healthy)
      sintoma <- image_to_mat(img_lesion)
      ################## no background #############
      if(is.null(img_background)){
        sadio_sintoma <-
          rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]) %>%
          transform(Y = ifelse(CODE == "img_healthy", 1, 0))
        sadio_sintoma$CODE <- NULL
        usef_area <- nrow(original$df_in)
        model <-
          glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma) %>%
          suppressWarnings()
        # isolate plant
        pred1 <- predict(model, newdata = original$df_in, type="response") %>% round(0)
        plant_symp <- matrix(pred1, ncol = ncol(original$R))
        plant_symp <- image_correct(plant_symp, perc = 0.01)
        ID <- c(plant_symp == 0)
        mpred1 <- bwlabel(plant_symp == 0)
        shape_leaf <-
          cbind(data.frame(computeFeatures.shape(mpred1)),
                data.frame(computeFeatures.moment(mpred1))[,1:2]
          )
        if(!is.null(lower_size)){
          shape_leaf <- shape_leaf[shape_leaf$s.area > lower_size, ]
        } else{
          shape_leaf <- shape_leaf[shape_leaf$s.area > 2, ]
        }
        if(!is.null(upper_size)){
          shape_leaf <- shape_leaf[shape_leaf$s.area < upper_size, ]
          shape_template <- shape_template[shape_template$s.area < upper_size, ]
        }
        shape_leaf$id <- 1:nrow(shape_leaf)
        shape_leaf <- shape_leaf[, c(9, 7, 8, 1, 2:6)]
        if(show_original == TRUE){
          if(is.null(col_background)){
            col_background <- col2rgb("green")
          } else{
            col_background <- col2rgb(col_background)
          }
          col_lesions <- col2rgb(col_lesions)
          im2 <- img
          im2@.Data[,,1][ID] <- col_lesions[1]
          im2@.Data[,,2][ID] <- col_lesions[2]
          im2@.Data[,,3][ID] <- col_lesions[3]
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
          col_lesions <- col2rgb(col_lesions)
          im2 <- img
          im2@.Data[,,1][ID] <- col_lesions[1]
          im2@.Data[,,2][ID] <- col_lesions[2]
          im2@.Data[,,3][ID] <- col_lesions[3]
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
          rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,],
                fundo$df_in[sample(1:nrow(fundo$df_in)),][1:nrows,]) %>%
          transform(Y = ifelse(CODE == "img_background", 0, 1))
        # fundo_resto$CODE <- NULL
        modelo1 <-
          glm(Y ~ R + G + B, family = binomial("logit"), data = fundo_resto) %>%
          suppressWarnings()
        pred1 <- predict(modelo1, newdata = original$df_in, type="response") %>% round(0)
        plant_background <- matrix(pred1, ncol = ncol(original$R))
        plant_background <- image_correct(plant_background, perc = 0.009)
        # image_show(plant_background)
        plant_background[plant_background == 1] <- 2
        sadio_sintoma <-
          rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]) %>%
          transform(Y = ifelse(CODE == "img_healthy", 1, 0))
        sadio_sintoma$CODE <- NULL
        modelo2 <-
          glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma) %>%
          suppressWarnings()
        # isolate plant
        ID <- c(plant_background == 2)
        usef_area <- nrow(original$df_in[ID,])
        pred2 <- predict(modelo2, newdata = original$df_in[ID,], type="response") %>% round(0)
        pred3 <- predict(modelo2, newdata = original$df_in, type="response") %>% round(0)
        pred3[!ID] <- 1
        leaf_sympts <- matrix(pred3, ncol = ncol(original$R))
        leaf_sympts <- image_correct(leaf_sympts, perc = 0.009)
        plant_background[leaf_sympts == 1] <- 3
        mpred1 <- bwlabel(leaf_sympts == 0)
        shape_leaf <-
          cbind(data.frame(computeFeatures.shape(mpred1)),
                data.frame(computeFeatures.moment(mpred1))[,1:2]
          )
        if(!is.null(lower_size)){
          shape_leaf <- shape_leaf[shape_leaf$s.area > lower_size, ]
        } else{
          shape_leaf <- shape_leaf[shape_leaf$s.area > 2, ]
        }
        if(!is.null(upper_size)){
          shape_leaf <- shape_leaf[shape_leaf$s.area < upper_size, ]
          shape_template <- shape_template[shape_template$s.area < upper_size, ]
        }
        shape_leaf$id <- 1:nrow(shape_leaf)
        shape_leaf <- shape_leaf[, c(9, 7, 8, 1, 2:6)]
        if(show_original == TRUE){
          im2 <- img
          col_lesions <- col2rgb(col_lesions)
          im2@.Data[,,1][ID][which(pred2 == 0)] <- col_lesions[1]
          im2@.Data[,,2][ID][which(pred2 == 0)] <- col_lesions[2]
          im2@.Data[,,3][ID][which(pred2 == 0)] <- col_lesions[3]
          image_show(im2)
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
          col_lesions <- col2rgb(col_lesions)
          col_leaf <- col2rgb(col_leaf)
          im2 <- img
          im2@.Data[,,1][ID][which(pred2 == 0)] <- col_lesions[1]
          im2@.Data[,,2][ID][which(pred2 == 0)] <- col_lesions[2]
          im2@.Data[,,3][ID][which(pred2 == 0)] <- col_lesions[3]
          im2@.Data[,,1][ID][which(pred2 != 0)] <- col_leaf[1]
          im2@.Data[,,2][ID][which(pred2 != 0)] <- col_leaf[2]
          im2@.Data[,,3][ID][which(pred2 != 0)] <- col_leaf[3]
          im2@.Data[,,1][!ID] <- col_background[1]
          im2@.Data[,,2][!ID] <- col_background[2]
          im2@.Data[,,3][!ID] <- col_background[3]
        }


      }
      if(show_image == TRUE){
        image_show(im2)
        text(shape_leaf[,2],
             shape_leaf[,3],
             shape_leaf$id,
             col = text_col,
             cex = text_size)
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada)
        }
        png(paste0(diretorio_processada, "/",
                   prefix,
                   name_ori, ".",
                   extens_ori),
            width = dim(im2@.Data)[1],
            height = dim(im2@.Data)[2])
        image_show(im2)
        text(shape_leaf[,2],
             shape_leaf[,3],
             shape_leaf$id,
             col = text_col,
             cex = text_size)
        dev.off()
      }
      stats <-
        data.frame(area = c(n = length(shape_leaf$s.area),
                            min(shape_leaf$s.area),
                            mean(shape_leaf$s.area),
                            max(shape_leaf$s.area),
                            sd(shape_leaf$s.area),
                            sum(shape_leaf$s.area),
                            sum(shape_leaf$s.area) /usef_area * 100),
                   perimeter = c(NA,
                                 min(shape_leaf$s.perimeter),
                                 mean(shape_leaf$s.perimeter),
                                 max(shape_leaf$s.perimeter),
                                 sd(shape_leaf$s.perimeter),
                                 sum(shape_leaf$s.perimeter),
                                 NA)) %>%
        transform(statistics = c("n", "min", "mean", "max", "sd", "sum", "prop"))
      stats <- stats[c(3, 1, 2)]
      results <- list(results = shape_leaf,
                      statistics = stats)
      return(results)
    }

  if(missing(img_pattern)){
    help_count(img, img_healthy, img_lesion, img_background, randomize,
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
    results <- list()
    pb <- progress(max = length(plants), style = 4)
    for (i in 1:length(plants)) {
      run_progress(pb, actual = i,
                   text = paste("Processing image", names_plant[i]))
      results[[i]] <-
        help_count(img  = names_plant[i],
                   img_healthy, img_lesion, img_background, randomize,
                   nrows, show_image, show_original, show_background, col_background,
                   save_image, dir_original, dir_processed)
    }
    names(results) <- names_plant
    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                results[[i]][["statistics"]] %>%
                  transform(id =  names(results[i])) %>%
                  .[,c(4, 1, 2, 3)]
              })
      )
    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                results[[i]][["results"]] %>%
                  transform(img =  names(results[i])) %>%
                  .[, c(10, 1:9)]
              })
      )
    return(list(statistics = stats,
                results = results))
  }
}

