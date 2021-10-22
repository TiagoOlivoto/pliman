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
                fill_hull, edge, extension, tolerance, object_size, show_image)
    } else{
      lapply(image, object_coord, id, index, invert,
             fill_hull, edge, extension, tolerance, object_size, show_image)
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
      parLapply(clust, image, object_contour, index, invert, fill_hull,
                watershed, extension, tolerance, object_size, show_image)
    } else{
      lapply(image, object_contour, index, invert, fill_hull,
             watershed, extension, tolerance, object_size, show_image)
    }
  } else{
    img2 <- image_binary(image,
                         index = index,
                         invert = invert,
                         fill_hull = fill_hull,
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

