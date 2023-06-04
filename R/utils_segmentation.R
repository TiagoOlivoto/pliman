#' Alternative watershed algorithm
#'
#' This is a basic watershed algorithm that can be used as a faster alternative
#' to [EBImage::watershed()]. I strongly suggest using this only with round
#' objects, since it doesn't consider both 'extension' and 'tolerance' arguments
#' of [EBImage::watershed()].
#'
#' @param binary A binary image
#' @param dist_thresh The distance threshold to create the
#' @param plot If `TRUE` (default) plots the labeled objects
#' @return The labelled version of `binary`.
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' binary <- image_binary(img, "B")[[1]]
#' wts <- watershed2(binary)
#' range(wts)
watershed2 <- function(binary,
                       dist_thresh = 0.75,
                       plot = TRUE){
  dt <- help_dist_transform(1 - binary)
  sure_fg <- dt > dist_thresh * max(dt)
  markers <- EBImage::bwlabel(sure_fg)
  wts <- EBImage::Image(help_watershed(binary, markers, dt))
  if(isTRUE(plot)){
    plot(EBImage::colorLabels(wts))
  }
  return(wts)
}


#' Distance map transform
#'
#' Computes the distance map transform of a binary image. The distance map is a
#' matrix which contains for each pixel the distance to its nearest background
#' pixel.
#'
#' @param binary A binary image
#'
#' @return An `Image` object or an array, with pixels containing the distances
#'   to the nearest background points
#' @export
#' @examples
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' binary <- image_binary(img, "B")[[1]]
#' wts <- dist_transform(binary)
#' range(wts)
dist_transform <- function(binary){
  help_dist_transform(1 - binary)
}


#' Labels objects
#'
#' All pixels for each connected set of foreground (non-zero) pixels in x are
#' set to an unique increasing integer, starting from 1. Hence, max(x) gives the
#' number of connected objects in x. This is a wrapper to [EBImage::bwlabel] or
#' [EBImage::watershed] (if `watershed = TRUE`).
#' @inheritParams image_binary
#' @inheritParams analyze_objects
#' @return A list with the same length of `img` containing the labeled objects.
#' @export
#'
#' @examples
#'img <- image_pliman("soybean_touch.jpg")
#'# segment the objects using the "B" (blue) band.
#'object_label(img, index = "B")
#'object_label(img, index = "B", watershed = TRUE)
object_label <- function(img,
                         index = "B",
                         invert = FALSE,
                         fill_hull = FALSE,
                         threshold = "Otsu",
                         k = 0.1,
                         windowsize = NULL,
                         filter = FALSE,
                         watershed = FALSE,
                         tolerance = NULL,
                         extension = NULL,
                         object_size = "medium",
                         plot = TRUE,
                         ncol = NULL,
                         nrow = NULL,
                         verbose = TRUE){
  img2 <- image_binary(img,
                       index = index,
                       invert = invert,
                       fill_hull = fill_hull,
                       threshold = threshold,
                       k = k,
                       windowsize = windowsize,
                       filter = filter,
                       resize = FALSE,
                       plot = FALSE)
  labels <- list()
  img2_len <- length(img2)
  for (i in 1:length(img2)){
    if(img2_len > 1){
      tmp <- img2[[i]][[1]]
    } else{
      tmp <- img2[[i]]
    }
    if(isTRUE(watershed)){
      parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
      res <- length(tmp)
      parms2 <- parms[parms$object_size == object_size,]
      rowid <-
        which(sapply(as.character(parms2$resolution), function(x) {
          eval(parse(text=x))}))
      ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
      tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
      labels[[i]] <- EBImage::watershed(EBImage::distmap(tmp),
                                        tolerance = tol,
                                        ext = ext)
    } else{
      labels[[i]] <- EBImage::bwlabel(tmp)
    }
  }
  if(plot == TRUE){
    num_plots <- length(labels)
    if (is.null(nrow) && is.null(ncol)){
      ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
      nrow <- ceiling(num_plots/ncol)
    }
    if (is.null(ncol)){
      ncol <- ceiling(num_plots/nrow)
    }
    if (is.null(nrow)){
      nrow <- ceiling(num_plots/ncol)
    }
    op <- par(mfrow = c(nrow, ncol))
    on.exit(par(op))
    index <- names(labels)
    for(i in 1:length(labels)){
      plot( EBImage::colorLabels(labels[[i]]))
      if(verbose == TRUE){
        dim <- image_dimension(labels[[i]], verbose = FALSE)
        text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
      }
    }
  }
  invisible(labels)
}

#' Calculate Otsu's threshold
#'
#' Given a numeric vector with the pixel's intensities, returns the threshold
#' value based on Otsu's method, which minimizes the combined intra-class
#' variance
#'
#' @param values A numeric vector with the pixel values.
#'
#' @return
#' A double (threshold value).
#'
#' @references Otsu, N. 1979. Threshold selection method from gray-level
#'   histograms. IEEE Trans Syst Man Cybern SMC-9(1): 62–66. doi:
#'   \doi{10.1109/tsmc.1979.4310076}

#' @export
#'
#' @examples
#' img <- image_pliman("soybean_touch.jpg")
#' thresh <- otsu(img@.Data[,,3])
#' plot(img[,,3] < thresh)
#'
otsu <- function(values){
  help_otsu(values)
}
