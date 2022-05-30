#' Utilities for working with image objects
#'
#' * `object_id()` get the object identification in an image.
#' * `object_coord()` get the object coordinates and (optionally) draw a
#' bounding rectangle around multiple objects in an image.
#' * `object_contour()` returns the coordinates (`x` and `y`) for the contours
#' of each object in the image.
#' * `object_isolate()` isolates an object from an image.
#' @name utils_objects
#' @param image An image of class `Image` or a list of `Image` objects.
#' @param id
#' * For `object_coord()`, a vector (or scalar) of object `id` to compute the
#' bounding rectangle. Object ids can be obtained with [object_id()]. Set `id =
#' "all"` to compute the coordinates for all objects in the image. If `id =
#' NULL` (default) a bounding rectangle is drawn including all the objects.
#' * For `object_isolate()`, a scalar that identifies the object to be extracted.
#' @param index The index to produce a binary image used to compute bounding
#'   rectangle coordinates. See [image_binary()] for more details.
#' @param invert Inverts the binary image, if desired. Defaults to `FALSE`.
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
#' @param show_image Shows the image with bounding rectangles? Defaults to
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
object_coord <- function(image,
                         id =  NULL,
                         index = "NB",
                         watershed = TRUE,
                         invert = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         edge = 2,
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         parallel = FALSE,
                         workers = NULL,
                         show_image = TRUE){
  check_ebi()
  if(inherits(image, "list")){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, image, object_coord, id, index, invert,
                fill_hull, threshold, edge, extension, tolerance,
                object_size, show_image)
    } else{
      lapply(image, object_coord, id, index, invert, fill_hull, threshold,
             edge, extension, tolerance, object_size, show_image)
    }
  } else{
    # helper function to get coordinates from a mask
    get_coordinates <- function(data_mask, edge){
      nrows <- nrow(data_mask)
      ncols <- ncol(data_mask)
      a <-
        apply(data_mask, 2, function(x){
          any(!is.na(x))
        })
      col_min <- min(which(a == TRUE))
      col_min <- col_min - edge
      col_min <- ifelse(col_min < 1, 1, col_min)
      col_max <- max(which(a == TRUE))
      col_max <- col_max + edge
      col_max <- ifelse(col_max > ncols, ncols, col_max)
      b <-
        apply(data_mask, 1, function(x){
          any(!is.na(x))
        })
      row_min <- min(which(b == TRUE))
      row_min <- row_min - edge
      row_min <- ifelse(row_min < 1, 1, row_min)
      row_max <- max(which(b == TRUE))
      row_max <- row_max + edge
      row_max <- ifelse(row_max > nrows, nrows, row_max)
      return(list(col_min = col_min,
                  col_max = col_max,
                  row_min = row_min,
                  row_max = row_max))
    }
    img2 <- image_binary(image,
                         index = index,
                         invert = invert,
                         fill_hull = fill_hull,
                         threshold = threshold,
                         show_image = FALSE,
                         resize = FALSE)[[1]]
    if(is.null(id)){
      data_mask <- img2@.Data
      data_mask[which(data_mask == FALSE)] <- NA
      coord <- get_coordinates(data_mask, edge)
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
        temp[which(data_mask != i)] <- NA
        list_mask[[i]] <- temp
      }
      list_mask <- list_mask[ids]
      coord <- sapply(list_mask, get_coordinates, edge)
      coord <- list(col_min = as.numeric(coord[1,]),
                    col_max = as.numeric(coord[2,]),
                    row_min = as.numeric(coord[3,]),
                    row_max = as.numeric(coord[4,]))
    }
    if(show_image == TRUE){
      plot(image)
      rect(xleft = coord$row_min,
           xright = coord$row_max,
           ybottom = coord$col_min,
           ytop = coord$col_max)
    }
    return(list(col_min = coord$col_min,
                col_max = coord$col_max,
                row_min = coord$row_min,
                row_max = coord$row_max))
  }
}
#' @name utils_objects
#' @export
#'
object_contour <- function(image,
                           index = "NB",
                           invert = FALSE,
                           fill_hull = FALSE,
                           threshold = "Otsu",
                           watershed = TRUE,
                           extension = NULL,
                           tolerance = NULL,
                           object_size = "medium",
                           parallel = FALSE,
                           workers = NULL,
                           show_image = TRUE){
  check_ebi()
  if(inherits(image, "list")){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, image, object_contour, index, invert, fill_hull, threshold,
                watershed, extension, tolerance, object_size, show_image)
    } else{
      lapply(image, object_contour, index, invert, fill_hull, threshold,
             watershed, extension, tolerance, object_size, show_image)
    }
  } else{
    img2 <- image_binary(image,
                         index = index,
                         invert = invert,
                         fill_hull = fill_hull,
                         threshold = threshold,
                         show_image = FALSE,
                         resize = FALSE)[[1]]
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
    dims <- sapply(contour, function(x){dim(x)[1]})
    contour <- contour[which(dims > mean(dims * 0.1))]
    if(isTRUE(show_image)){
      plot(image)
      plot_contour(contour, col = "red")
    }
    return(lapply(contour, function(x){data.frame(x)}))
  }
}
#' @name utils_objects
#' @export
object_isolate <- function(image,
                           id = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           ...){
  if(inherits(image, "list")){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, image, object_isolate, id, ...)
    } else{
      lapply(image, object_isolate, id, ...)
    }
  } else{
    coord <- object_coord(image,
                          id = id,
                          show_image = FALSE,
                          ...)
    segmented <- image[coord$row_min:coord$row_max,
                       coord$col_min:coord$col_max,
                       1:3]
    return(segmented)
  }
}
#' @name utils_objects
#' @export
object_id <- function(image,
                      parallel = FALSE,
                      workers = NULL,
                      ...){
  if(inherits(image, "list")){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      parLapply(clust, image, object_id, ...)
    } else{
      lapply(image, object_id, ...)
    }
  } else{
    analyze_objects(image, verbose = FALSE, marker = "id", ...)
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
#' @param keep_location A logical argument (defaults to `TRUE`). If `FALSE`, the
#'   new image is created with the object in the exactly position of the
#'   original image.
#' @param col_background The background color in the new image. Defaults to
#'   `white`. Use the built-in color names which `R` knows about (see
#'   ?[grDevices::colors()]) or a numeric vector with R, G, and B intensities
#'   (see ?[grDevices::rgb()]).
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param ... Additional arguments passed on to [image_combine()]
#' @return A list of objects of class `Image`.
#' @export
#' @seealso [analyze_objects()], [image_binary()]
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("la_leaves.jpg", plot = TRUE)
#' imgs <- object_split(img, workers = 2) # set to NULL to use 50% of the cores
#'
object_split <- function(img,
                         index = "NB",
                         lower_size = NULL,
                         my_index = NULL,
                         watershed = TRUE,
                         invert = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         extension = NULL,
                         tolerance = NULL,
                         object_size = "medium",
                         keep_location = FALSE,
                         col_background = "white",
                         show_image = TRUE,
                         verbose = TRUE,
                         workers = NULL,
                         ...){

  img2 <- image_binary(img,
                       index = index,
                       my_index = my_index,
                       invert = invert,
                       fill_hull = fill_hull,
                       threshold = threshold,
                       resize = FALSE,
                       show_image = FALSE)[[1]]
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

  objcts <-
    sapply(1:max(nmask), function(i){
      length(which(nmask == i))
    })
  av_area <- mean(objcts)
  ifelse(!is.null(lower_size),
         cutsize <- lower_size,
         cutsize <-  av_area * 0.1)
  ifelse(is.character(col_background),
         col_background <- col2rgb(col_background) / 255,
         col_background <- col_background / 255)
  selected <- which(objcts > cutsize)

  list_crop <- function(img, selected, keep_location, col_background){
    id = which(nmask != selected)
    img@.Data[, , 1][id] <- col_background[1]
    img@.Data[, , 2][id] <- col_background[2]
    img@.Data[, , 3][id] <- col_background[3]
    if(!isTRUE(keep_location)){
      image_autocrop(img, plot = FALSE, index = "R")
    }
  }
  nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
  clust <- makeCluster(nworkers)
  clusterExport(clust,
                varlist = c("img", "nmask", "list_crop", "image_autocrop", "selected", "col_background"),
                envir=environment())
  on.exit(stopCluster(clust))
  list_objects <-
    parLapply(clust, seq_along(selected),
              function(i){
                list_crop(img, selected[i], keep_location, col_background)
              })

  names(list_objects) <- 1:length(list_objects)
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
  if(isTRUE(show_image)){
    image_combine(list_objects, ...)
  }
  return(list_objects)
}
