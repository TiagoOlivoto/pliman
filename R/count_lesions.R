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
#'   images that are named as 1.-, 2.-, and so on.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time, especially when `img_pattern` is used is informed. The
#'   number of sections is set up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param invert Inverts the binary image, if desired. This is useful to process
#'   images with black background. Defaults to `FALSE`.
#' @param index,my_index A character value specifying the target mode for
#'   conversion to binary image when `img_healthy` and `img_lesion` are not
#'   declared. Defaults to `"NB"` (normalized blue). See [image_index()] for
#'   more details.
#' @param lower_size Lower limit for size for the image analysis. Leaf images
#'   often contain dirt and dust. To prevent dust from affecting the image
#'   analysis, the lower limit of analyzed size is set to 0.1, i.e., objects
#'   with lesser than 10% of the mean of all objects are removed. One can set a
#'   known area or use `lower_limit = 0` to select all objects (not advised).
#' @param upper_size Upper limit for size for the image analysis. Defaults to
#'   `NULL`, i.e., no upper limit used.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param lesion_size The size of the lesion. Used to automatically set up
#'   `tolerance` and `extension` parameters. One of the following. `"small"` (2-5
#'   mm in diameter, e.g, rust pustules), `"medium"` (0.5-1.0 cm in diameter,
#'   e.g, wheat leaf spot), `"large"` (1-2 cm in diameter, and  `"elarge"` (2-3
#'   cm in diameter, e.g, target spot of soybean).
#' @param segment If `TRUE` (Default) implements the Watershed Algorithm to
#'   segment lesions connected by a fairly few pixels that could be considered
#'   as two distinct lesions. If `FALSE`, lesions that are connected by any
#'   pixel are considered unique lesions. For more details see
#'   [EBImage::watershed()].
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest. Defaults to `NULL`, i.e., starting values are set up according to the argument `lesion_size`.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Defaults to 20. Higher value smooths out small
#'   objects.
#' @param show_segmentation Shows the object segmentation colored with random
#'   permutations. Defaults to `TRUE`.
#' @param show_image Show image after processing?
#' @param show_original Show the symptoms in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param col_leaf Leaf color after image processing. Defaults to `"green"`
#' @param col_lesions Symptoms color after image processing. Defaults to `"red"`.
#' @param col_background Background color after image processing. Defaults to
#'   `"NULL"`.
#' @param marker,marker_col,marker_size The type, color and size of the object
#'   marker. Defaults to `NULL`, which shows a red point when `show_segmentation
#'   = FALSE`. To force a marker to be used with segmented objects, set up to
#'   `marker = "point"` (to show a point) or `marker = "text"` to enumerate the
#'   objects.
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
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return A list with the following objects:
#'  * `results` A data frame with the results (area, perimeter, radius) for
#'  object.
#'  * `statistics` A data frame with the summary statistics for the image.
#'  * `count` (If `img_pattern` is used), summarizing the count number for each
#'  image.
#' @export
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(image_pliman("sev_leaf_nb.jpg"))
#' healthy <- image_import(image_pliman("sev_healthy.jpg"))
#' lesions <- image_import(image_pliman("sev_sympt.jpg"))
#' image_combine(img, healthy, lesions, ncol = 3)
#' a <-
#' count_lesions(img = img,
#'               img_healthy = healthy,
#'               img_lesion = lesions,
#'               lesion_size = "elarge", # extra large lesions
#'               show_image = TRUE,
#'               show_segmentation = FALSE,
#'               marker = "text")
#' }
#'
count_lesions <- function(img,
                          img_healthy = NULL,
                          img_lesion = NULL,
                          img_background = NULL,
                          img_pattern = NULL,
                          parallel = FALSE,
                          workers = NULL,
                          resize = FALSE,
                          invert = FALSE,
                          index = "NB",
                          my_index = NULL,
                          lower_size = NULL,
                          upper_size = NULL,
                          randomize = TRUE,
                          nrows = 10000,
                          lesion_size = "medium",
                          segment = TRUE,
                          tolerance = NULL,
                          extension = NULL,
                          show_segmentation = TRUE,
                          show_image = FALSE,
                          show_original = TRUE,
                          show_background = TRUE,
                          col_leaf = NULL,
                          col_lesions = NULL,
                          col_background = NULL,
                          marker = NULL,
                          marker_col = NULL,
                          marker_size = NULL,
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
  help_count <-
    function(img, img_healthy, img_lesion, img_background, resize, invert,
             index, my_index, lesion_size, tolerance, extension,
             randomize, nrows, show_image, show_original, show_background,
             col_leaf, col_lesions, col_background,
             save_image, dir_original, dir_processed){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "jpg"
      }
      backg <- !is.null(col_background)
      col_background <- col2rgb(ifelse(is.null(col_background), "white", col_background))
      col_lesions <- col2rgb(ifelse(is.null(col_lesions), "red", col_lesions))
      col_leaf <- col2rgb(ifelse(is.null(col_leaf), "green", col_leaf))
      if(!is.null(img_healthy) && !is.null(img_lesion)){
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
            transform(rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                            sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]),
                      Y = ifelse(CODE == "img_healthy", 1, 0))
          sadio_sintoma$CODE <- NULL
          usef_area <- nrow(original$df_in)
          model <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
          # isolate plant
          pred1 <- round(predict(model, newdata = original$df_in, type="response"), 0)
          plant_symp <- matrix(pred1, ncol = ncol(original$R))
          plant_symp <- 1 - image_correct(plant_symp, perc = 0.01)
          ID <- c(plant_symp == 0)
          ID2 <- c(plant_symp == 1)
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE),
                            header = T, sep = ";")
          parms2 <- parms[parms$object_size == lesion_size,]
          res <- length(plant_symp)
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          ifelse(segment == FALSE,
                 nmask <- bwlabel(plant_symp),
                 nmask <- watershed(distmap(plant_symp),
                                    tolerance = tol,
                                    ext = ext))
          if(show_original == TRUE & show_segmentation == FALSE){
            im2 <- img
            im2@.Data[,,1][!ID] <- col_lesions[1]
            im2@.Data[,,2][!ID] <- col_lesions[2]
            im2@.Data[,,3][!ID] <- col_lesions[3]
            if(backg){
              im3 <- colorLabels(nmask)
              im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- img@.Data[,,1][which(im3@.Data[,,1]==0)]
              im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- img@.Data[,,2][which(im3@.Data[,,2]==0)]
              im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- img@.Data[,,3][which(im3@.Data[,,3]==0)]
            }
          }
          if(show_original == TRUE & show_segmentation == TRUE){
            im2 <- colorLabels(nmask)
            if(backg){
              im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
              im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
              im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
            } else{
              im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- img@.Data[,,1][which(im2@.Data[,,1]==0)]
              im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- img@.Data[,,2][which(im2@.Data[,,2]==0)]
              im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- img@.Data[,,3][which(im2@.Data[,,3]==0)]
            }
          }
          if(show_original == FALSE){
            if(show_segmentation == TRUE){
              im2 <- colorLabels(nmask)
              im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_leaf[1]
              im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_leaf[2]
              im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_leaf[3]
            } else{
              im2 <- img
              im2@.Data[,,1][!ID] <- col_lesions[1]
              im2@.Data[,,2][!ID] <- col_lesions[2]
              im2@.Data[,,3][!ID] <- col_lesions[3]
              im2@.Data[,,1][ID] <- col_leaf[1]
              im2@.Data[,,2][ID] <- col_leaf[2]
              im2@.Data[,,3][ID] <- col_leaf[3]
            }
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
          modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"),
                                          data = fundo_resto))
          pred1 <- round(predict(modelo1, newdata = original$df_in, type="response"), 0)
          plant_background <- matrix(pred1, ncol = ncol(original$R))
          plant_background <- image_correct(plant_background, perc = 0.009)
          plant_background[plant_background == 1] <- 2
          sadio_sintoma <-
            transform(rbind(sadio$df_in[sample(1:nrow(sadio$df_in)),][1:nrows,],
                            sintoma$df_in[sample(1:nrow(sintoma$df_in)),][1:nrows,]),
                      Y = ifelse(CODE == "img_healthy", 1, 0))
          sadio_sintoma$CODE <- NULL
          modelo2 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"),
                                          data = sadio_sintoma))
          # isolate plant
          ID <- c(plant_background == 2)
          usef_area <- nrow(original$df_in[ID,])
          pred2 <- round(predict(modelo2, newdata = original$df_in[ID,], type="response"), 0)
          pred3 <- round(predict(modelo2, newdata = original$df_in, type="response"), 0)
          pred3[!ID] <- 1
          leaf_sympts <- matrix(pred3, ncol = ncol(original$R))
          leaf_sympts <- 1 - image_correct(leaf_sympts, perc = 0.009)
          plant_background[leaf_sympts == 1] <- 3
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE),
                            header = T, sep = ";")
          parms2 <- parms[parms$object_size == lesion_size,]
          res <- length(leaf_sympts)
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          ifelse(segment == FALSE,
                 nmask <- bwlabel(leaf_sympts),
                 nmask <- watershed(distmap(leaf_sympts),
                                    tolerance = tol,
                                    ext = ext))
          if(show_original == TRUE & show_segmentation == TRUE){
            im2 <- colorLabels(nmask)
            if(backg){
              im2@.Data[,,1][!ID] <- col_background[1]
              im2@.Data[,,2][!ID] <- col_background[2]
              im2@.Data[,,3][!ID] <- col_background[3]
              im2@.Data[,,1][ID][which(pred2 != 0)] <- img@.Data[,,1][ID][which(pred2 != 0)]
              im2@.Data[,,2][ID][which(pred2 != 0)] <- img@.Data[,,2][ID][which(pred2 != 0)]
              im2@.Data[,,3][ID][which(pred2 != 0)] <- img@.Data[,,3][ID][which(pred2 != 0)]
            } else{
              im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- img@.Data[,,1][which(im2@.Data[,,1]==0)]
              im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- img@.Data[,,2][which(im2@.Data[,,2]==0)]
              im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- img@.Data[,,3][which(im2@.Data[,,3]==0)]
            }
          }
          if(show_original == TRUE & show_segmentation == FALSE){
            im2 <- img
            im2@.Data[,,1][ID][which(pred2 == 0)] <- col_lesions[1]
            im2@.Data[,,2][ID][which(pred2 == 0)] <- col_lesions[2]
            im2@.Data[,,3][ID][which(pred2 == 0)] <- col_lesions[3]
            if(backg){
              im2@.Data[,,1][!ID] <- col_background[1]
              im2@.Data[,,2][!ID] <- col_background[2]
              im2@.Data[,,3][!ID] <- col_background[3]
            }
          }
          if(show_original == FALSE){
            if(show_segmentation == TRUE){
              im2 <- colorLabels(nmask)
              im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
              im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
              im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
              im2@.Data[,,1][ID][which(pred2 != 0)] <- col_leaf[1]
              im2@.Data[,,2][ID][which(pred2 != 0)] <- col_leaf[2]
              im2@.Data[,,3][ID][which(pred2 != 0)] <- col_leaf[3]
            } else{
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
        }
      } else{
        imgs <- img[[1]][["image"]][,,1:3]
        img2 <- image_binary(imgs,
                             index = index,
                             my_index = my_index,
                             invert = invert,
                             resize = FALSE,
                             show_image = FALSE)[[1]]
        img2@.Data[which(imgs[,,1]==1)] <- FALSE

        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        res <- length(img2)
        usef_area <- res
        parms2 <- parms[parms$object_size == lesion_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        nmask <- watershed(distmap(img2),
                           tolerance = tol,
                           ext = ext)
        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
        if(show_original == TRUE & show_segmentation == FALSE){
          im2 <- imgs
            im2@.Data[,,1][which(imgs[,,1]==1)] <- col_background[1]
            im2@.Data[,,2][which(imgs[,,2]==1)] <- col_background[2]
            im2@.Data[,,3][which(imgs[,,3]==1)] <- col_background[3]
            im2@.Data[,,1][ID] <- col_lesions[1]
            im2@.Data[,,2][ID] <- col_lesions[2]
            im2@.Data[,,3][ID] <- col_lesions[3]
            image_show(im2)
        }
        if(show_original == TRUE & show_segmentation == TRUE){
          im2 <- colorLabels(nmask)
            im2@.Data[,,1][which(imgs[,,1]==1)] <- col_background[1]
            im2@.Data[,,2][which(imgs[,,2]==1)] <- col_background[2]
            im2@.Data[,,3][which(imgs[,,3]==1)] <- col_background[3]
            im2@.Data[,,1][ID2] <- imgs@.Data[,,1][ID2]
            im2@.Data[,,2][ID2] <- imgs@.Data[,,2][ID2]
            im2@.Data[,,3][ID2] <- imgs@.Data[,,3][ID2]
        }
        if(show_original == FALSE){
          if(show_segmentation == TRUE){
            im2 <- colorLabels(nmask)
            im2@.Data[,,1][ID2] <- col_leaf[1]
            im2@.Data[,,2][ID2] <- col_leaf[2]
            im2@.Data[,,3][ID2] <- col_leaf[3]
            im2@.Data[,,1][which(imgs[,,1]==1)] <- col_background[1]
            im2@.Data[,,2][which(imgs[,,2]==1)] <- col_background[2]
            im2@.Data[,,3][which(imgs[,,3]==1)] <- col_background[3]
          } else{
            im2 <- imgs
            im2@.Data[,,1][ID2] <- col_leaf[1]
            im2@.Data[,,2][ID2] <- col_leaf[2]
            im2@.Data[,,3][ID2] <- col_leaf[3]
            im2@.Data[,,1][ID] <- col_lesions[1]
            im2@.Data[,,2][ID] <- col_lesions[2]
            im2@.Data[,,3][ID] <- col_lesions[3]
            im2@.Data[,,1][which(imgs[,,1]==1)] <- col_background[1]
            im2@.Data[,,2][which(imgs[,,2]==1)] <- col_background[2]
            im2@.Data[,,3][which(imgs[,,3]==1)] <- col_background[3]
          }
        }
      }

      shape <-
        cbind(data.frame(computeFeatures.shape(nmask)),
              data.frame(computeFeatures.moment(nmask))[,1:2]
        )
      ifelse(!is.null(lower_size),
             shape <- shape[shape$s.area > lower_size, ],
             shape <- shape[shape$s.area > mean(shape$s.area) * 0.1, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$s.area < upper_size, ]
      }
      shape$id <- 1:nrow(shape)
      shape <- shape[, c(9, 7, 8, 1, 2:6)]
      show_mark <- !is.null(marker) && show_segmentation == TRUE | show_segmentation == FALSE
      marker <- ifelse(is.null(marker), "point", marker)
      marker_col <- ifelse(is.null(marker_col), "white", marker_col)
      marker_size <- ifelse(is.null(marker_size), 0.9, marker_size)
      if(show_image == TRUE){
        if(marker == "text"){
          image_show(im2)
          if(show_mark){
            text(shape[,2],
                 shape[,3],
                 shape[,1],
                 col = marker_col,
                 cex = marker_size)
          }
        } else{
          image_show(im2)
          if(show_mark){
            points(shape[,2],
                   shape[,3],
                   col = marker_col,
                   pch = 16,
                   cex = marker_size)
          }
        }
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
        if(marker == "text"){
          image_show(im2)
          if(show_mark){
            text(shape[,2],
                 shape[,3],
                 shape[,1],
                 col = marker_col,
                 cex = marker_size)
          }
        } else{
          image_show(im2)
          if(show_mark){
            points(shape[,2],
                   shape[,3],
                   col = marker_col,
                   pch = 16,
                   cex = marker_size)
          }
        }
        dev.off()
      }
      stats <-
        transform(data.frame(area = c(n = length(shape$s.area),
                                      min(shape$s.area),
                                      mean(shape$s.area),
                                      max(shape$s.area),
                                      sd(shape$s.area),
                                      sum(shape$s.area),
                                      sum(shape$s.area) /usef_area * 100),
                             perimeter = c(NA,
                                           min(shape$s.perimeter),
                                           mean(shape$s.perimeter),
                                           max(shape$s.perimeter),
                                           sd(shape$s.perimeter),
                                           sum(shape$s.perimeter),
                                           NA)),
                  statistics = c("n", "min", "mean", "max", "sd", "sum", "prop"))
      stats <- stats[c(3, 1, 2)]
      shape <- shape[,c(1:6, 8:9, 7)]
      shape <- transform(shape, radius_ratio = s.radius.max / s.radius.min)
      colnames(shape) <- c("id", "x", "y", "area", "perimeter", "radius_mean",
                           "radius_min", "radius_max", "radius_sd", "radius_ratio")
      results <- list(results = shape,
                      statistics = stats)
      class(results) <- "plm_count"
      if(verbose == TRUE){
        print(results$statistics, row.names = FALSE)
      }
      return(results)
    }
  if(missing(img_pattern)){
    help_count(img, img_healthy, img_lesion, img_background, resize, invert, index, my_index,
               lesion_size, tolerance, extension, randomize, nrows, show_image,
               show_original, show_background, col_leaf, col_lesions, col_background,
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
                    varlist = c("names_plant", "help_count", "file_name",
                                "check_names_dir", "file_extension", "image_import",
                                "image_binary", "watershed", "distmap", "computeFeatures.moment",
                                "computeFeatures.shape", "colorLabels", "image_show",
                                "image_to_mat", "image_correct", "bwlabel"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_count(x,
                               img_healthy, img_lesion, img_background, resize, invert, index, my_index,
                               lesion_size, tolerance, extension, randomize, nrows, show_image,
                               show_original, show_background, col_leaf, col_lesions, col_background,
                               save_image, dir_original, dir_processed)
                  })

    } else{
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        if(verbose == TRUE){
        run_progress(pb, actual = i,
                     text = paste("Processing image", names_plant[i]))
        }
        results[[i]] <-
          help_count(img  = names_plant[i],
                     img_healthy, img_lesion, img_background, randomize,
                     nrows, show_image, show_original, show_background, col_background,
                     save_image, dir_original, dir_processed)
      }
    }
    names(results) <- names_plant
    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                  transform(results[[i]][["statistics"]],
                            id =  names(results[i]))[,c(4, 1, 2, 3)]
              })
      )
    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                  transform(results[[i]][["results"]],
                            img =  names(results[i]))[, c(10, 1:9)]
              })
      )
    return(list(statistics = stats,
                results = results))
  }
}

