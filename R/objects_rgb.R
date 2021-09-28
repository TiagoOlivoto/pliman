#' Get Red Green and Blue for image objects
#'
#' Get the Red Green and Blue (RGB) for objects in an image. This function will
#' be deprecated in a near future. Use [analyze_objects()] with the argument
#' `object_index` instead.
#'
#'A binary image is first generated to segment the foreground and background.
#'The argument index is useful to choose a proper index to segment the image
#'(see [image_binary()] for more details). Then, the number of objects in the
#'foreground is counted. Change `tolerance` and `extension` values to better set
#'up watershed-based object detection. If color palettes samples are provided, a
#'general linear model (binomial family) fitted to the RGB values is used to
#'segment fore- and background. For each segmented object, the RGB values are
#'obtained with. Users can also compute an index for each object using the
#'argument `object_index`, useful to classify objects based on its RGB values.
#'
#'By using `img_pattern` it is possible to process several images with common
#'pattern names that are stored in the current working directory or in the
#'subdirectory informed in `dir_original`'. To speed up the computation time,
#'one can set `parallel = TRUE`.
#'
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
#'   number of sections is set up to 70% of available cores.
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
#' @param object_index The same as `index`, used to compute the index for each
#'   object in the image.
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
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing? Defaults to `TRUE`.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named as `proc_*` where `*` is the image name
#'   given in `img`.
#' @param prefix The prefix to be included in the processed images. Defaults to
#'   `"proc_"`.
#' @param marker,marker_col,marker_size,marker_digits The marker, color, size
#'   and significant digits of the object marker. Defaults to `marker =
#'   "index"`, which shows the object index. Set to `marker = "id"` to show the
#'   object id.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory.
#' @param verbose If `FALSE`, runs the code silently.
#'
#' @return A list with the following objects.
#'  * `objects` A data frame with the measures for each object.
#'  * `rgb` A data frame with the Red, Green and Blue values for each object
#'  * `indexes` A data frame with the index computed according to the argument
#'  `object_index`.
#' @export
#'
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' # objects_rgb() will be deprecated in the future.
#' # Use analyze_objects() with the argument `object_index` instead.
#' }
objects_rgb <- function(img,
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
                        object_index = "B",
                        object_size = "large",
                        tolerance = NULL,
                        extension = NULL,
                        lower_size = NULL,
                        upper_size = NULL,
                        topn_lower = NULL,
                        topn_upper = NULL,
                        nrows = 10000,
                        show_image = TRUE,
                        save_image = FALSE,
                        prefix = "proc_",
                        marker = NULL,
                        marker_col = NULL,
                        marker_size = NULL,
                        marker_digits = NULL,
                        dir_original = NULL,
                        dir_processed = NULL,
                        verbose = TRUE){
  check_ebi()
  message("`objects_rgb()` will be deprecated in the future. Use `analyze_objects()` with the argument `object_index` instead.")
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
    function(img, foreground, background, resize, tolerance, extension,show_image,
             save_image, dir_original, dir_processed, verbose){
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
      if(resize != FALSE){
        img <- EBImage::resize(img, resize)
      }
      if(!is.null(foreground) && !is.null(background)){
        if(is.character(foreground)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = foreground)
          check_names_dir(foreground, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          foreground <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
          if(resize != FALSE){
            foreground <- EBImage::resize(foreground, resize)
          }
        }
        if(is.character(background)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = background)
          check_names_dir(background, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          background <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
          if(resize != FALSE){
            background <- EBImage::resize(background, resize)
          }
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
        nmask <- EBImage::watershed(EBImage::distmap(foreground_background),
                           tolerance = tol,
                           ext = ext)
      } else{
        img2 <- image_binary(img,
                             index = index,
                             invert = invert,
                             fill_hull = fill_hull,
                             show_image = FALSE,
                             resize = FALSE)[[1]]
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        res <- length(img2)
        parms2 <- parms[parms$object_size == object_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        nmask <- EBImage::watershed(EBImage::distmap(img2),
                           tolerance = tol,
                           ext = ext)
      }
      data_mask <- nmask@.Data
      get_rgb <- function(img, data_mask, index){
        data.frame(object = index,
                   R = img@.Data[,,1][which(data_mask == index)],
                   G = img@.Data[,,2][which(data_mask == index)],
                   B = img@.Data[,,3][which(data_mask == index)])
      }
      rgb_objects <-
        do.call(rbind,
                lapply(1:max(data_mask), function(i){
                  get_rgb(img, data_mask, i)
                }))
      shape <-
        cbind(data.frame(EBImage::computeFeatures.shape(nmask)),
              data.frame(EBImage::computeFeatures.moment(nmask))[,1:2]
        )
      shape$id <- 1:nrow(shape)
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
      shape <- shape[, c(9, 7, 8, 1, 2:6)][,c(1:6, 8:9, 7)]
      colnames(shape) <- c("id", "x", "y", "area", "perimeter", "radius_mean",
                           "radius_min", "radius_max", "radius_sd")
      rgb_objects <- subset(rgb_objects, object %in% shape$id)
      indexes <- transform(rgb_objects, index = eval(parse(text = object_index)))
      indexes <- aggregate(index ~ object, data = indexes, FUN = mean)
      marker_col <- ifelse(is.null(marker_col), "white", marker_col)
      marker_size <- ifelse(is.null(marker_size), 0.9, marker_size)
      if(show_image == TRUE){
        plot(img)
        if(!is.null(marker)){
          if(marker != "index"){
            text(shape[,2],
                 shape[,3],
                 shape[,1],
                 col = marker_col,
                 cex = marker_size)
          } else{
            text(shape[,2],
                 shape[,3],
                 round(indexes[,2], 2),
                 col = marker_col,
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
            width = dim(img@.Data)[1],
            height = dim(img@.Data)[2])
        plot(img)
        if(!is.null(marker)){
          if(marker != "index"){
            text(shape[,2],
                 shape[,3],
                 shape[,1],
                 col = marker_col,
                 cex = marker_size)
          } else{
            text(shape[,2],
                 shape[,3],
                 round(indexes[,2], 2),
                 col = marker_col,
                 cex = marker_size)

          }
        }
        dev.off()
      }
      return(structure(list(objects = shape,
                            rgb = rgb_objects,
                            indexes = indexes),
                       class = "objects_rgb"))
    }
  if(missing(img_pattern)){
    help_count(img, foreground, background, resize, tolerance, extension,show_image,
               save_image, dir_original, dir_processed, verbose)
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
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff.\nExtensions found:", paste(extensions, sep = ", "))
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_count", "file_name",
                                "check_names_dir", "file_extension", "image_import",
                                "image_binary", "image_resize", "detectCores", "makeCluster", "clusterExport",
                                "stopCluster", "parLapply"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_count(x,
                               foreground, background, resize, tolerance, extension,
                               show_image, save_image, dir_original, dir_processed, verbose =  FALSE)
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
                     foreground, background, resize, tolerance, extension,
                     show_image, save_image, dir_original, dir_processed,
                     verbose)
      }
    }
    names(results) <- names_plant
    objects <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["objects"]],
                          img =  names(results[i]))[,c(10, 1:9)]
              })
      )
    indexes <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["indexes"]],
                          img =  names(results[i]))[, c(3, 1:2)]
              })
      )
    invisible(list(objects = objects,
                   indexes = indexes))
  }
}




