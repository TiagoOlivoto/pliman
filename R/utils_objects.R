#' Utilities for working with image objects
#'
#' * `object_id()` get the object identification in an image.
#' * `object_coord()` get the object coordinates and (optionally) draw a
#' bounding rectangle around multiple objects in an image.
#' * `object_contour()` returns the coordinates (`x` and `y`) for the contours
#' of each object in the image.
#' * `object_isolate()` isolates an object from an image.
#' @name utils_objects
#'
#' @inheritParams analyze_objects
#' @param img An image of class `Image` or a list of `Image` objects.
#' @param center If `TRUE` returns the object contours centered on the origin.
#' @param id
#' * For `object_coord()`, a vector (or scalar) of object `id` to compute the
#' bounding rectangle. Object ids can be obtained with [object_id()]. Set `id =
#' "all"` to compute the coordinates for all objects in the image. If `id =
#' NULL` (default) a bounding rectangle is drawn including all the objects.
#' * For `object_isolate()`, a scalar that identifies the object to be extracted.
#'
#' @param dir_original The directory containing the original images. Defaults
#'    to `NULL`, which means that the current working directory will be
#'    considered.
#' @param index The index to produce a binary image used to compute bounding
#'   rectangle coordinates. See [image_binary()] for more details.
#' @param invert Inverts the binary image, if desired. Defaults to `FALSE`.
#' @param filter Performs median filtering in the binary image? See more at
#'   [image_filter()]. Defaults to `FALSE`. Use a positive integer to define the
#'   size of the median filtering. Larger values are effective at removing
#'   noise, but adversely affect edges.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param watershed If `TRUE` (default) performs watershed-based object
#'   detection. This will detect objects even when they are touching one other.
#'   If `FALSE`, all pixels for each connected set of foreground pixels are set
#'   to a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param threshold By default (`threshold = "Otsu"`), a threshold value based
#'   on Otsu's method is used to reduce the grayscale image to a binary image.
#'   If a numeric value is informed, this value will be used as a threshold.
#'   Inform any non-numeric value different than "Otsu" to iteratively chosen
#'   the threshold based on a raster plot showing pixel intensity of the index.
#' @param edge The number of pixels in the edge of the bounding rectangle.
#'   Defaults to `2`.
#' @param extension,tolerance,object_size Controls the watershed segmentation of
#'   objects in the image. See [analyze_objects()] for more details.
#' @param plot Shows the image with bounding rectangles? Defaults to
#'   `TRUE`.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 50% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param ...
#' * For `object_isolate()`, further arguments passed on to [object_coord()].
#' * For `object_id()`, further arguments passed on to [analyze_objects()].
#' @return
#' * `object_id()` An image of class `"Image"` containing the object's
#' identification.
#' * `object_coord()` A list with the coordinates for the bounding rectangles.
#' If `id = "all"` or a numeric vector, a list with a vector of coordinates is
#' returned.
#' * `object_isolate()` An image of class `"Image"` containing the isolated
#' object.
#' @export
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg")
#' # Get the object's (leaves) identification
#' object_id(img)
#'
#' # Get the coordinates and draw a bounding rectangle around leaves 1 and 3
#' object_coord(img, id = c(1, 3))
#'
#' # Isolate leaf 3
#' isolated <- object_isolate(img, id = 3)
#' plot(isolated)
#'
#' }
object_coord <- function(img,
                         id =  NULL,
                         index = "NB",
                         watershed = TRUE,
                         invert = FALSE,
                         filter = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         edge = 2,
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         parallel = FALSE,
                         workers = NULL,
                         plot = TRUE){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, img, object_coord, id, index, invert,
                fill_hull, threshold, edge, extension, tolerance,
                object_size, plot)
    } else{
      lapply(img, object_coord, id, index, invert, fill_hull, threshold,
             edge, extension, tolerance, object_size, plot)
    }
  } else{
    img2 <- help_binary(img,
                        index = index,
                        invert = invert,
                        filter = filter,
                        fill_hull = fill_hull,
                        threshold = threshold)
    if(is.null(id)){
      data_mask <- img2@.Data
      coord <- t(as.matrix(bounding_box(data_mask, edge)))
      colnames(coord) <- c("xleft", "xright", "ybottom", "ytop")
      if(plot == TRUE){
        plot(img)
        rect(xleft = coord[1],
             xright = coord[2],
             ybottom = coord[3],
             ytop = coord[4])
      }
    } else{
      if(isTRUE(watershed)){
        res <- length(img2)
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        parms2 <- parms[parms$object_size == object_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        nmask <- EBImage::watershed(EBImage::distmap(img2),
                                    tolerance = tol,
                                    ext = ext)
      } else{
        nmask <- EBImage::bwlabel(img2)
      }
      data_mask <- nmask@.Data
      ifelse(id == "all",
             ids <- 1:max(data_mask),
             ids <- id)
      list_mask <- list()
      for (i in ids) {
        temp <- data_mask
        temp[which(data_mask != i)] <- FALSE
        list_mask[[i]] <- temp
      }
      list_mask <- list_mask[ids]
      coord <- t(sapply(list_mask, bounding_box, edge))
      colnames(coord) <- c("xleft", "xright", "ybottom", "ytop")
      if(plot == TRUE){
        plot(img)
        rect(xleft = coord[,1],
             xright = coord[,2],
             ybottom = coord[,3],
             ytop = coord[,4])
      }
    }
    return(coord)
  }
}
#' @name utils_objects
#' @inheritParams analyze_objects
#' @export
#'

object_contour <- function(img,
                           pattern = NULL,
                           dir_original = NULL,
                           center =  FALSE,
                           index = "NB",
                           invert = FALSE,
                           filter = FALSE,
                           fill_hull = FALSE,
                           threshold = "Otsu",
                           watershed = TRUE,
                           extension = NULL,
                           tolerance = NULL,
                           object_size = "medium",
                           parallel = FALSE,
                           workers = NULL,
                           plot = TRUE,
                           verbose = TRUE){
  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }

  if(is.null(pattern) && inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, img, object_contour, pattern, dir_original, center, index, invert, filter, fill_hull, threshold,
                watershed, extension, tolerance, object_size, plot = plot)
    } else{
      lapply(img, object_contour, pattern, dir_original, center, index, invert, filter, fill_hull, threshold,
             watershed, extension, tolerance, object_size, plot = plot)
    }
  } else{
    if(is.null(pattern)){
      img2 <- help_binary(img,
                          index = index,
                          invert = invert,
                          filter = filter,
                          fill_hull = fill_hull,
                          threshold = threshold)
      if(isTRUE(watershed)){
        res <- length(img2)
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        parms2 <- parms[parms$object_size == object_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        nmask <- EBImage::watershed(EBImage::distmap(img2),
                                    tolerance = tol,
                                    ext = ext)
      } else{
        nmask <- EBImage::bwlabel(img2)
      }
      contour <- EBImage::ocontour(nmask)
      if(isTRUE(center)){
        contour <-
          lapply(contour, function(x){
            transform(x,
                      X1 = X1 - mean(X1),
                      X2 = X2 - mean(X2))
          })
      }
      dims <- sapply(contour, function(x){dim(x)[1]})
      contour <- contour[which(dims > mean(dims * 0.1))]
      if(isTRUE(plot)){
        if(isTRUE(center)){
          plot_polygon(contour)
        } else{
          plot(img)
          plot_contour(contour, col = "red")
        }
      }
      return(contour)
    } else{
      if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
        pattern <- "^[0-9].*$"
      }
      plants <- list.files(pattern = pattern, diretorio_original)
      extensions <- as.character(sapply(plants, file_extension))
      names_plant <- as.character(sapply(plants, file_name))
      if(length(grep(pattern, names_plant)) == 0){
        stop(paste("Pattern '", pattern, "' not found in '",
                   paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
             call. = FALSE)
      }
      if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
        stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
      }


      help_contour <- function(img){
        img <- image_import(img)
        img2 <- help_binary(img,
                            index = index,
                            invert = invert,
                            filter = filter,
                            fill_hull = fill_hull,
                            threshold = threshold)
        if(isTRUE(watershed)){
          res <- length(img2)
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
          parms2 <- parms[parms$object_size == object_size,]
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          nmask <- EBImage::watershed(EBImage::distmap(img2),
                                      tolerance = tol,
                                      ext = ext)
        } else{
          nmask <- EBImage::bwlabel(img2)
        }
        contour <- EBImage::ocontour(nmask)
        if(isTRUE(center)){
          contour <-
            lapply(contour, function(x){
              transform(x,
                        X1 = X1 - mean(X1),
                        X2 = X2 - mean(X2))
            })
        }
        dims <- sapply(contour, function(x){dim(x)[1]})
        contour[which(dims > mean(dims * 0.1))]
      }


      if(parallel == TRUE){
        init_time <- Sys.time()
        nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.5), workers)
        cl <- parallel::makePSOCKcluster(nworkers)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))

        if(verbose == TRUE){
          message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
        }
        ## declare alias for dopar command
        `%dopar%` <- foreach::`%dopar%`
        results <-
          foreach::foreach(i = seq_along(plants), .packages = c("pliman", "EBImage")) %dopar%{
            help_contour(plants[[i]])
          }

      } else{

        pb <- progress(max = length(plants), style = 4)
        foo <- function(plants, ...){
          run_progress(pb, ...)
          help_contour(plants)
        }
        results <-
          lapply(seq_along(plants), function(i){
            foo(plants[i],
                actual = i,
                text = paste("Processing image", names_plant[i]))
          })

      }
      names(results) <- plants
      return(results)
    }
  }
}



#' @name utils_objects
#' @export
object_isolate <- function(img,
                           id = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           ...){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, img, object_isolate, id, ...)
    } else{
      lapply(img, object_isolate, id, ...)
    }
  } else{
    coord <- object_coord(img,
                          id = id,
                          plot = FALSE,
                          ...)
    segmented <- img[coord[1]:coord[2],
                     coord[3]:coord[4],
                     1:3]
    return(segmented)
  }
}
#' @name utils_objects
#' @export
object_id <- function(img,
                      parallel = FALSE,
                      workers = NULL,
                      ...){
  if(inherits(img, "list")){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, img, object_id, ...)
    } else{
      lapply(img, object_id, ...)
    }
  } else{
    analyze_objects(img, verbose = FALSE, marker = "id", ...)
  }
}




#' Splits objects from an image into multiple images
#'
#' Using threshold-based segmentation, objects are first isolated from
#' background. Then, a new image is created for each single object. A list of
#' images is returned.
#'
#' @inheritParams analyze_objects
#' @param lower_size Plant images often contain dirt and dust. To prevent dust from
#'   affecting the image analysis, objects with lesser than 10% of the mean of all objects
#'   are removed. Set `lower_limit = 0` to keep all the objects.
#' @param edge The number of pixels to be added in the edge of the segmented
#'   object. Defaults to 5.
#' @param remove_bg If `TRUE`, the pixels that are not part of objects are
#'   converted to white.
#' @param ... Additional arguments passed on to [image_combine()]
#' @return A list of objects of class `Image`.
#' @export
#' @seealso [analyze_objects()], [image_binary()]
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg", plot = TRUE)
#' imgs <- object_split(img) # set to NULL to use 50% of the cores
#'
object_split <- function(img,
                         index = "NB",
                         lower_size = NULL,
                         watershed = TRUE,
                         invert = FALSE,
                         fill_hull = FALSE,
                         filter = 2,
                         threshold = "Otsu",
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         edge = 3,
                         remove_bg = FALSE,
                         plot = TRUE,
                         verbose = TRUE,
                         ...){

  img2 <- help_binary(img,
                      filter = filter,
                      index = index,
                      invert = invert,
                      fill_hull = fill_hull,
                      threshold = threshold)
  if(isTRUE(watershed)){
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
  } else{
    nmask <- EBImage::bwlabel(img2)
  }

  objcts <- get_area_mask(nmask)
  av_area <- mean(objcts)
  ifelse(!is.null(lower_size),
         cutsize <- lower_size,
         cutsize <-  av_area * 0.1)
  selected <- which(objcts > cutsize)

  split_objects <- function(img, nmask){
    objects <- help_isolate_object(img[,,1], img[,,2], img[,,3], nmask, remove_bg, edge)
    lapply(seq_along(objects), function(x){
      dimx <- dim(objects[[x]][[1]])
      EBImage::Image(array(c(objects[[x]][[1]], objects[[x]][[2]], objects[[x]][[3]]), dim = c(dimx, 3)), colormode = "Color")
    })
  }
  list_objects <- split_objects(img, nmask)
  names(list_objects) <- 1:length(list_objects)
  list_objects <- list_objects[selected]
  if(isTRUE(verbose)){
    cat("==============================\n")
    cat("Summary of the procedure\n")
    cat("==============================\n")
    cat("Number of objects:", length(objcts), "\n")
    cat("Average area     :", mean(objcts), "\n")
    cat("Minimum area     :", min(objcts), "\n")
    cat("Maximum area     :", max(objcts), "\n")
    cat("Objects created  :", length(list_objects), "\n")
    cat("==============================\n")
  }
  if(isTRUE(plot)){
    image_combine(list_objects, ...)
  }
  return(list_objects)
}

#' Export multiple objects from an image to multiple images
#'
#' Givin an image with multiple objects, `object_export()` will split the
#' objects into a list of objects using [object_split()] and then export them to
#' multiple images into the current working directory (or a subfolder). Batch
#' processing is performed by declaring a file name pattern that matches the
#' images within the working directory.
#'
#' @inheritParams object_split
#' @inheritParams utils_image
#' @inheritParams analyze_objects
#'
#' @param pattern A pattern of file name used to identify images to be
#'   processed. For example, if `pattern = "im"` all images in the current
#'   working directory that the name matches the pattern (e.g., img1.-,
#'   image1.-, im2.-) will be imported and processed. Providing any number as
#'   pattern (e.g., `pattern = "1"`) will select images that are named as 1.-,
#'   2.-, and so on. An error will be returned if the pattern matches any file
#'   that is not supported (e.g., img1.pdf).
#'@param dir_original The directory containing the original images. Defaults to
#'  `NULL`. It can be either a full path, e.g., `"C:/Desktop/imgs"`, or a
#'  subfolder within the current working directory, e.g., `"/imgs"`.
#' @param dir_processed Optional character string indicating a subfolder within the
#'   current working directory to save the image(s). If the folder doesn't
#'   exist, it will be created.
#' @param format The format of image to be exported.
#' @param squarize Squarizes the image before the exportation? If `TRUE`,
#'   [image_square()] will be called internally.
#' @return A `NULL` object.
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg")
#' object_export(img)
#'
#' }
object_export <- function(img,
                          pattern = NULL,
                          dir_original = NULL,
                          dir_processed = NULL,
                          format = ".jpg",
                          squarize = FALSE,
                          index = "NB",
                          lower_size = NULL,
                          watershed = TRUE,
                          invert = FALSE,
                          fill_hull = FALSE,
                          filter = 2,
                          threshold = "Otsu",
                          extension = NULL,
                          tolerance = NULL,
                          object_size = "medium",
                          edge = 20,
                          remove_bg = FALSE,
                          parallel = FALSE,
                          verbose = TRUE){
  if(is.null(pattern)){
    list_objects <- object_split(img = img,
                                 index = index,
                                 lower_size = lower_size,
                                 watershed = watershed,
                                 invert = invert,
                                 fill_hull = fill_hull,
                                 filter = filter,
                                 threshold = threshold,
                                 extension = extension,
                                 tolerance = tolerance,
                                 object_size = object_size,
                                 edge = edge,
                                 remove_bg = remove_bg,
                                 plot = FALSE,
                                 verbose = FALSE)

    a <- lapply(seq_along(list_objects), function(i){
      tmp <- list_objects[[i]]
      if(isTRUE(squarize)){
        tmp <- image_square(tmp,
                            plot = FALSE,
                            sample_left = 10,
                            sample_top = 10,
                            sample_right = 10,
                            sample_bottom = 10)
      }
      image_export(tmp,
                   name = paste0("obj_",
                                 leading_zeros(as.numeric(names(list_objects[i])),
                                               n = 4), format),
                   subfolder = dir_processed)
    })
  } else{

    if(is.null(dir_original)){
      diretorio_original <- paste0("./")
    } else{
      diretorio_original <-
        ifelse(grepl("[/\\]", dir_original),
               dir_original,
               paste0("./", dir_original))
    }
    if(is.null(dir_processed)){
      diretorio_processada <- paste0("./")
    } else{
      diretorio_processada <-
        ifelse(grepl("[/\\]", dir_processed),
               dir_processed,
               paste0("./", dir_processed))
    }

    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("Pattern '", pattern, "' not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }

    if(isTRUE(parallel)){

      init_time <- Sys.time()
      nworkers <- trunc(detectCores()*.3)
      cl <- parallel::makePSOCKcluster(nworkers)
      doParallel::registerDoParallel(cl)
      on.exit(stopCluster(cl))

      if(verbose == TRUE){
        message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
      }
      ## declare alias for dopar command
      `%dopar%` <- foreach::`%dopar%`

      results <-
        foreach::foreach(i = seq_along(plants), .packages = c("pliman")) %dopar%{

          tmpimg <- image_import(plants[[i]], path = diretorio_original)

          list_objects <- object_split(img = tmpimg,
                                       index = index,
                                       lower_size = lower_size,
                                       watershed = watershed,
                                       invert = invert,
                                       fill_hull = fill_hull,
                                       filter = filter,
                                       threshold = threshold,
                                       extension = extension,
                                       tolerance = tolerance,
                                       object_size = object_size,
                                       edge = edge,
                                       remove_bg = remove_bg,
                                       verbose = FALSE,
                                       plot = FALSE)

          a <- lapply(seq_along(list_objects), function(j){
            tmp <- list_objects[[j]]
            if(isTRUE(squarize)){
              tmp <- image_square(tmp,
                                  plot = FALSE,
                                  sample_left = 10,
                                  sample_top = 10,
                                  sample_right = 10,
                                  sample_bottom = 10)
            }
            image_export(tmp,
                         name = paste0(file_name(plants[[i]]), "_",
                                       leading_zeros(as.numeric(names(list_objects[j])),
                                                     n = 4), format),
                         subfolder = diretorio_processada)
          }
          )
        }

      message("Done!")
      message("Elapsed time: ", sec_to_hms(as.numeric(difftime(Sys.time(),  init_time, units = "secs"))))

    } else{

      for(i in seq_along(plants)){
        tmpimg <- image_import(plants[[i]], path = diretorio_original)

        list_objects <- object_split(img = tmpimg,
                                     index = index,
                                     lower_size = lower_size,
                                     watershed = watershed,
                                     invert = invert,
                                     fill_hull = fill_hull,
                                     filter = filter,
                                     threshold = threshold,
                                     extension = extension,
                                     tolerance = tolerance,
                                     object_size = object_size,
                                     edge = edge,
                                     remove_bg = remove_bg,
                                     verbose = FALSE,
                                     plot = FALSE)

        a <- lapply(seq_along(list_objects), function(j){
          tmp <- list_objects[[j]]
          if(isTRUE(squarize)){
            tmp <- image_square(tmp,
                                plot = FALSE,
                                sample_left = 10,
                                sample_top = 10,
                                sample_right = 10,
                                sample_bottom = 10)
          }
          image_export(tmp,
                       name = paste0(file_name(plants[[i]]), "_",
                                     leading_zeros(as.numeric(names(list_objects[j])),
                                                   n = 4), format),
                       subfolder = diretorio_processada)
        }
        )
      }
    }

  }
}



#' Extract red, green and blue values from objects
#'
#' Given an image and a matrix of labels that identify each object, the function
#' extracts the red, green, and blue values from each object.
#'
#' @param img An `Image` object
#' @param labels A mask containing the labels for each object. This can be
#'   obtained with [EBImage::bwlabel()] or [EBImage::watershed()]
#'
#' @return A data.frame with `n` rows (number of pixels for all the objects) and
#'   the following columns:
#'  * `id`: the object id;
#'  * `R`: the value for the red band;
#'  * `G`: the value for the blue band;
#'  * `B`: the value for the green band;
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' # segment the objects using the "B" (blue) band (default)
#'
#' labs <- object_label(img, watershed = TRUE)
#' rgb <- object_rgb(img, labs[[1]])
#' head(rgb)
object_rgb <- function(img, labels){
  dd <- help_get_rgb(img[,,1], img[,,2], img[,,3], labels)
  df2 <- data.frame(do.call(rbind,  lapply(dd, function(x){
    matrix(x, ncol = 4, byrow = TRUE)
  })))
  colnames(df2) <- c("id", "R", "G", "B")
  return(df2)
}



#' Apply color to image objects
#'
#' The function applies the color informed in the argument `color` to segmented
#' objects in the image. The segmentation is performed using image indexes. Use
#' [image_index()] to identify the better candidate index to segment objects.
#'
#' @inheritParams image_binary
#' @param color The color to apply in the image objects. Defaults to `"blue"`.
#' @param plot Plots the modified image? Defaults to `TRUE`.
#' @param ... Additional arguments passed on to [image_binary()].
#'
#' @return An object of class `Image`
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg")
#' img2 <- object_to_color(img, index = "G-R")
#' image_combine(img, img2)
#'
object_to_color <- function(img,
                            index = "NB",
                            color = "blue",
                            plot = TRUE,
                            ...){
  bin <- help_binary(img,
                     index = index,
                     ...)
  pix_ref <- which(bin == 1)
  colto <- col2rgb(color) / 255
  img@.Data[,,1][pix_ref] <- colto[1]
  img@.Data[,,2][pix_ref] <- colto[2]
  img@.Data[,,3][pix_ref] <- colto[3]
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}
