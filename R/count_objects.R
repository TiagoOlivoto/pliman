#' Computes number of objects in an image
#'
#'Counts the number of objects in an image. See more at details.
#'
#'Counts the number of objects in an image. A binary image is first generated to
#'segment the foreground and background. The argument index is useful to choose
#'a proper index to segment the image (see [image_binary()] for more details).
#'Then, the number of objects in the foreground is counted. By setting up
#'arguments such as `lower_size`, `upper_size` is possible to set a threshold
#'for lower and upper sizes of the objects, respectively.  Change `tolerance`
#'and `extension` values to better set up watershed-based object detection. If
#'color palettes samples are provided, a general linear model (binomial family)
#'fitted to the RGB values is used to segment fore- and background.
#'
#'By using `img_pattern` it is possible to process several images with common
#'pattern names that are stored in the current working directory or in the
#'subdirectory informed in `dir_original`'. To speed up the computation time,
#'one can set `parallel = TRUE`.
#' @param img The image to be analyzed.
#' @param foreground A color palette of the foreground (optional).
#' @param background A color palette of the background (optional).
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
#' @param fill_hull Fill holes in the binary image? Defaults to `FALSE`. This is
#'   useful to fill holes in objects that have portions with a color similar to
#'   the background. IMPORTANT: Objects touching each other can be combined into
#'   one single object, which may underestimate the number of objects in an
#'   image.
#' @param invert Inverts the binary image, if desired. This is useful to process
#'   images with black background. Defaults to `FALSE`.
#' @param index,my_index A character value specifying the target mode for
#'   conversion to binary image when `foreground` and `background` are not
#'   declared. Defaults to `"NB"` (normalized blue). See [image_index()] for
#'   more details.
#' @param object_size The size of the object. Used to automatically set up
#'   `tolerance` and `extension` parameters. One of the following. `"small"`
#'   (e.g, wheat grains), `"medium"` (e.g, soybean grains), `"large"`(e.g,
#'   peanut grains), and `"elarge"` (e.g, soybean pods)`.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Defaults to 20. Higher value smooths out small
#'   objects.
#' @param lower_size,upper_size Lower and upper limits for size for the image
#'   analysis. Plant images often contain dirt and dust. To prevent dust from
#'   affecting the image analysis, objects with lesser than 10% of the mean of
#'   all objects are removed. Upper limit is set to `NULL`, i.e., no upper
#'   limit used. One can set a known area or use `lower_limit = 0` to select all
#'   objects (not advised). Objects that matches the size of a given range of
#'   sizes can be selected by setting up the two arguments. For example, if
#'   `lower_size = 120` and `upper_size = 140`, objects with size greater than
#'   or equal 120 and less than or equal 140 will be considered.
#' @param topn_lower,topn_upper Select the top `n` objects based on its area.
#'   `topn_lower` selects the `n` elements with the smallest area whereas
#'   `topn_upper` selects the `n` objects with the largest area.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing?
#' @param show_original Show the count objects in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param show_segmentation Shows the object segmentation colored with random
#'   permutations. Defaults to `TRUE`.
#' @param col_foreground,col_background Foreground and background color after
#'   image processing. Defaults to `NULL`, in which `"black"`, and `"white"` are
#'   used, respectively.
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
#' @import EBImage
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(image_pliman("soybean_touch.jpg"))
#' count_objects(img)
#'
#' # Enumerate the objects in the original image
#' count_objects(img,
#'               show_segmentation = FALSE,
#'               marker = "text",
#'               marker_col = "white")
#' }
count_objects <- function(img,
                          foreground = NULL,
                          background = NULL,
                          img_pattern = NULL,
                          parallel = FALSE,
                          workers = NULL,
                          resize = FALSE,
                          fill_hull = FALSE,
                          invert = FALSE,
                          index = "NB",
                          my_index = NULL,
                          object_size = "medium",
                          tolerance = NULL,
                          extension = NULL,
                          lower_size = NULL,
                          upper_size = NULL,
                          topn_lower = NULL,
                          topn_upper = NULL,
                          randomize = TRUE,
                          nrows = 10000,
                          show_image = TRUE,
                          show_original = TRUE,
                          show_background = TRUE,
                          show_segmentation = TRUE,
                          col_foreground = NULL,
                          col_background = NULL,
                          marker = NULL,
                          marker_col = NULL,
                          marker_size = NULL,
                          save_image = FALSE,
                          prefix = "proc_",
                          dir_original = NULL,
                          dir_processed = NULL,
                          verbose = TRUE){
  if(!object_size %in% c("small", "medium", "large", "elarge")){
    stop("'object_size' must be one of 'small', 'medium', 'large', or 'elarge'")
  }
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
    function(img, foreground, background, resize, fill_hull, tolerance, extension,
             randomize, nrows, show_image, show_original, show_background, marker,
             marker_col, marker_size, save_image, prefix,
             dir_original, dir_processed, verbose){
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
      if(resize != FALSE){
        img <- image_resize(img, resize)
      }
      if(!is.null(foreground) && !is.null(background)){
        if(is.character(foreground)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = foreground)
          check_names_dir(foreground, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          foreground <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        if(is.character(background)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = background)
          check_names_dir(background, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          background <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        original <- image_to_mat(img)
        foreground <- image_to_mat(foreground)
        background <- image_to_mat(background)
        back_fore <-
          transform(rbind(foreground$df_in[sample(1:nrow(foreground$df_in)),][1:nrows,],
                          background$df_in[sample(1:nrow(background$df_in)),][1:nrows,]),
                    Y = ifelse(CODE == "background", 0, 1))
        modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = back_fore))
        pred1 <- round(predict(modelo1, newdata = original$df_in, type="response"), 0)
        foreground_background <- matrix(pred1, ncol = ncol(original$R))
        foreground_background <- image_correct(foreground_background, perc = 0.02)
        ID <- c(foreground_background == 1)
        ID2 <- c(foreground_background == 0)
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        res <- length(foreground_background)
        parms2 <- parms[parms$object_size == object_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        nmask <- watershed(distmap(foreground_background),
                           tolerance = tol,
                           ext = ext)
      } else{
        img2 <- image_binary(img,
                             index = index,
                             my_index = my_index,
                             invert = invert,
                             fill_hull = fill_hull,
                             resize = FALSE,
                             show_image = FALSE)[[1]]
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        res <- length(img2)
        parms2 <- parms[parms$object_size == object_size,]
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
      }
      backg <- !is.null(col_background)
      col_background <- col2rgb(ifelse(is.null(col_background), "white", col_background))
      col_foreground <- col2rgb(ifelse(is.null(col_foreground), "black", col_foreground))
      if(show_original == TRUE & show_segmentation == FALSE){
        im2 <- img
        if(backg){
          im3 <- colorLabels(nmask)
          im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- col_background[1]
          im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- col_background[2]
          im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- col_background[3]
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
      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        stop("Only one of 'lower_*' or 'topn_*' can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$s.area > lower_size, ],
             shape <- shape[shape$s.area > mean(shape$s.area) * 0.1, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$s.area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$s.area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$s.area, decreasing = TRUE),][1:topn_upper,]
      }
      shape$id <- 1:nrow(shape)
      shape <- shape[, c(9, 7, 8, 1, 2:6)]
      show_mark <- !is.null(marker) && show_segmentation == TRUE | show_segmentation == FALSE
      marker <- ifelse(is.null(marker), "text", marker)
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
          marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
          image_show(im2)
          text(shape[,2],
               shape[,3],
               shape[,1],
               col = marker_col,
               cex = marker_size)
        } else{
          marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
          image_show(im2)
          text(shape[,2],
               shape[,3],
               col = marker_col,
               pch = 16,
               cex = marker_size)

        }
        dev.off()
      }
      stats <-
        transform(data.frame(area = c(n = length(shape$s.area),
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
                                           sum(shape$s.perimeter))),
                  statistics = c("n", "min", "mean", "max", "sd", "sum"))
      stats <- stats[c(3, 1, 2)]
      shape <- shape[,c(1:6, 8:9, 7)]
      colnames(shape) <- c("id", "x", "y", "area", "perimeter", "radius_mean",
                           "radius_min", "radius_max", "radius_sd")
      results <- list(results = shape,
                      statistics = stats)
      class(results) <- "plm_count"
      if(verbose == TRUE){
        cat("\n--------------------------------------------\n")
        cat("Number of objects:", stats[1,2],"\n")
        cat("--------------------------------------------\n")
        print(stats[-1,], row.names = FALSE)
        cat("\n")
      }
      invisible(results)
    }
  if(missing(img_pattern)){
    help_count(img, foreground, background, resize, fill_hull, tolerance , extension, randomize,
               nrows, show_image, show_original, show_background, marker,
               marker_col, marker_size, save_image, prefix,
               dir_original, dir_processed, verbose)
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
                                "computeFeatures.shape", "colorLabels", "image_show"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }

      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_count(x,
                               foreground, background, resize, fill_hull,
                               tolerance , extension, randomize,
                               nrows, show_image, show_original, show_background, marker,
                               marker_col, marker_size, save_image, prefix,
                               dir_original, dir_processed, verbose =  FALSE)
                  })

    } else{
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        run_progress(pb, actual = i,
                     text = paste("Processing image", names_plant[i]))
        results[[i]] <-
          help_count(img  = names_plant[i],
                     foreground, background, resize, fill_hull, tolerance,
                     extension, randomize, nrows, show_image, show_original,
                     show_background, marker, marker_col, marker_size, save_image,
                     prefix, dir_original, dir_processed, verbose)
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
    summ <- stats[stats$statistics == "n",c(1,3)]
    if(verbose == TRUE){
      names(summ) <- c("Image", "Objects")
      cat("--------------------------------------------\n")
      print(summ)
      cat("--------------------------------------------\n")
      message("Done!")

    }
    invisible(list(statistics = stats,
                   count = summ,
                   results = results))
  }
}

