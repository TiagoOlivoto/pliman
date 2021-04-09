#' Computes number of objects in an image
#'
#'Counts the number of objects in an image. See more at details.
#'
#'Counts the number of objects in an image. A binary image is first generated
#'to segment the foreground and background. The argument channel is useful to
#'choose a proper RGB band to segment the image Then, the number of objects in
#'the foreground is counted. By setting up arguments such as lower_size,
#'upper_size is possible to set a threshold for lower and upper sizes of the
#'objects, respectively.  Change tolerance and extension values to better set up
#'watershed-based object detection.
#'
#'If color palettes samples are provided, a general linear model (binomial
#'family) fitted to the RGB values is used to segment fore- and background. By
#'using `img_pattern` it is possible to process several images with common
#'pattern names that are stored in the current working directory or in the
#'subdirectory informed in '`dir_original`.
#' @param img The image to be analyzed.
#' @param foreground A color palette of the foreground (optional.
#' @param background A color palette of the background (optional).
#' @param img_pattern A pattern of file name used to identify images to be
#'   processed. For example, if `img_pattern = "im"` all images that the name
#'   matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
#'   Providing any number as pattern (e.g., `img_pattern = "1"`) will select
#'   images that are nammed as 1.-, 2.-, and so on.
#' @param channel A character value specifying the target mode for conversion to
#'   binary image. One of `gray`, `grey`, `red`, `green`, or `blue`.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest. Tolerance should be chosen according to the range of
#'   x. Defaults to `3`.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Defaults to 20. Higher value smoothes out small
#'   objects.
#' @param lower_size,upper_size Lower and upper limits for size for the image
#'   analysis. Plant images often contain dirt and dust. To prevent dust from
#'   affecting the image analysis, objects with lesser than 1% of the mean of
#'   all objects are removed. Upper limit is set to `NULL`, i.e., no upper
#'   limitused. One can set a known area or use `lower_limit = 0` to select all
#'   objects (not advised). Objects that matches the size of a given range of
#'   sizes can be selected by setting up the two arguments. For example, if
#'   `lower_size = 120` and `upper_size = 140`, objects with size greater than
#'   or equal 120 and less than or equal 140 will be considered.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing?
#' @param show_original Show the count objects in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param show_segmentation Shows the object segmentation colored with random
#'   permutations.
#' @param col_foreground,col_background Foreground and background color after
#'   image processing. Defaults to `"black"`, and `"white"`, respectively.
#' @param marker,marker_col,marker_size The type, color and size of the object
#'   marker. Defaults to a red point. Use `marker = "text"` to enumerate the
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
#' @return A data frame with the results for each image.
#' @export
#' @importFrom EBImage channel combine watershed distmap otsu getFrames colorLabels resize
#' @md
#' @examples
#' \donttest{
#' library(pliman)
#' }
#'
count_objects <- function(img,
                          foreground = NULL,
                          background = NULL,
                          img_pattern = NULL,
                          channel = "blue",
                          tolerance = NULL,
                          extension = NULL,
                          lower_size = NULL,
                          upper_size = NULL,
                          randomize = TRUE,
                          nrows = 10000,
                          show_image = TRUE,
                          show_original = TRUE,
                          show_background = TRUE,
                          show_segmentation = FALSE,
                          col_foreground = "black",
                          col_background = "white",
                          marker = "point",
                          marker_col = "red",
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
    function(img, foreground, background, tolerance, extension, randomize,
             nrows, show_image, show_original, show_background, marker, marker_col, marker_size, save_image, prefix,
             dir_original, dir_processed){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img)
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- import_image(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      if(!is.null(foreground) && !is.null(background)){
        if(is.character(foreground)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = foreground)
          check_names_dir(foreground, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          foreground <- import_image(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        if(is.character(background)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = background)
          check_names_dir(background, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          background <- import_image(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        original <- image_to_mat(img, randomize = randomize, nrows = nrows)
        foreground <- image_to_mat(foreground, randomize = randomize, nrows = nrows)
        background <- image_to_mat(background, randomize = randomize, nrows = nrows)
        back_fore <-
          rbind(foreground$df_man,
                background$df_man) %>%
          transform(Y = ifelse(CODE == "background", 0, 1))
        modelo1 <-
          glm(Y ~ R + G + B, family = binomial("logit"), data = back_fore) %>%
          suppressWarnings()
        pred1 <- predict(modelo1, newdata = original$df_in, type="response") %>% round(0)
        foreground_background <- matrix(pred1, ncol = ncol(original$R))
        foreground_background <- correct_image(foreground_background, perc = 0.01)
        ID <- c(foreground_background == 1)
        ID2 <- c(foreground_background == 0)
        tol <- ifelse(is.null(tolerance), round(nrow(foreground_background) / 300, 0), tolerance)
        ext <- ifelse(is.null(extension), round(nrow(foreground_background) / 300, 0), extension)
        nmask <- watershed(distmap(foreground_background),
                           tolerance = tol,
                           ext = ext)
      } else{
        img2 <- channel(img, channel)
        threshold <- otsu(img2)
        img2 <- combine(mapply(function(frame, th) frame < th, getFrames(img2), threshold, SIMPLIFY=FALSE))
        tol <- ifelse(is.null(tolerance), round(nrow(img2) / 300, 0), tolerance)
        ext <- ifelse(is.null(extension), round(nrow(img2) / 300, 0), extension)
        nmask <- watershed(distmap(img2),
                           tolerance = tol,
                           ext = ext)
        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
        feat <- computeFeatures.moment(nmask)
      }
      if(show_original == TRUE & show_segmentation == FALSE){
        im2 <- img
      } else{
        col_foreground <- col2rgb(col_foreground)
        col_background <- col2rgb(col_background)
        if(show_segmentation == TRUE){
          im2 <- colorLabels(nmask)
          im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
          im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
          im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
        } else{
          im2 <- img
          im2@.Data[,,1][ID] <- col_foreground[1]
          im2@.Data[,,2][ID] <- col_foreground[2]
          im2@.Data[,,3][ID] <- col_foreground[3]
          im2@.Data[,,1][ID2] <- col_background[1]
          im2@.Data[,,2][ID2] <- col_background[2]
          im2@.Data[,,3][ID2] <- col_background[3]
        }
      }
      shape <-
        cbind(data.frame(computeFeatures.shape(nmask)),
              data.frame(computeFeatures.moment(nmask))[,1:2]
        )
      ifelse(!is.null(lower_size),
             shape <- shape[shape$s.area > lower_size, ],
             shape <- shape[shape$s.area > mean(shape$s.area) * 0.01, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$s.area < upper_size, ]
      }
      # return(shape)
      shape$id <- 1:nrow(shape)
      shape <- shape[, c(9, 7, 8, 1, 2:6)]
      if(show_image == TRUE){
        if(marker == "text"){
          marker_size <- ifelse(is.null(marker_size), 0.9, marker_size)
          show_image(im2)
          text(shape[,2],
               shape[,3],
               shape[,1],
               col = marker_col,
               cex = marker_size)
        } else{
          marker_size <- ifelse(is.null(marker_size), 0.9, marker_size)
          show_image(im2)
          points(shape[,2],
                 shape[,3],
                 col = marker_col,
                 pch = 16,
                 cex = marker_size)
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
          marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
          show_image(im2)
          text(feat[,1],
               feat[,2],
               seq(1:nrow(feat)),
               col = marker_col,
               cex = marker_size)
        } else{
          marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
          show_image(im2)
          points(feat[,1],
                 feat[,2],
                 col = marker_col,
                 pch = 16,
                 cex = marker_size)

        }
        dev.off()
      }
      stats <-
        data.frame(area = c(n = length(shape$s.area),
                            min(shape$s.area),
                            mean(shape$s.area),
                            max(shape$s.area),
                            sd(shape$s.area),
                            sum(shape$s.area)),
                   perimeter = c(NA,
                                 min(shape$s.perimeter),
                                 mean(shape$s.perimeter),
                                 max(shape$s.perimeter),
                                 sd(shape$s.perimeter),
                                 sum(shape$s.perimeter))) %>%
        transform(statistics = c("n", "min", "mean", "max", "sd", "sum"))
      stats <- stats[c(3, 1, 2)]
      results <- list(results = shape,
                      statistics = stats)
      if(verbose == TRUE){
        cat("--------------------------------------------\n")
        cat("Number of objects:", stats[1,2],"\n")
        cat("--------------------------------------------\n")
        print(stats[-1,], row.names = FALSE)
        cat("--------------------------------------------\n")
      }
      invisible(results)
    }

  if(missing(img_pattern)){
    help_count(img, foreground, background, tolerance, extension, randomize,
               nrows, show_image, show_original, show_background, marker,
               marker_col, marker_size, save_image, prefix,
               dir_original, dir_processed)
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
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    results <- list()
    pb <- progress(max = length(plants), style = 4)
    for (i in 1:length(plants)) {
      run_progress(pb, actual = i,
                   text = paste("Processing image", names_plant[i]))
      results[[i]] <-
        help_count(img  = names_plant[i],
                   foreground, background, tolerance, extension, randomize,
                   nrows, show_image, show_original, show_background, marker,
                   marker_col, marker_size, save_image, prefix,
                   dir_original, dir_processed)
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
    invisible(list(statistics = stats,
                   results = results))
  }
}

