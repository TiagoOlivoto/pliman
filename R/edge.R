#' Object edges
#'
#' Applies the Sobel-Feldman Operator to detect edges. The operator is based on
#' convolving the image with a small, separable, and integer-valued filter in
#' the horizontal and vertical directions.
#'
#' @param img An image or a list of images of class `Image`.
#' @param sigma Gaussian kernel standard deviation used in the gaussian blur.
#' @param threshold The theshold method to be used.  If `threshold = "Otsu"`
#'   (default), a threshold value based on Otsu's method is used to reduce the
#'   grayscale image to a binary image. If any non-numeric value different than
#'   `"Otsu"` is used, an iterative section will allow you to choose the
#'   threshold based on a raster plot showing pixel intensity of the index.
#'   Alternatively, provide a numeric value to be used as the threshold value.
#' @param thinning Logical value indicating whether a thinning procedure should
#'   be applied to the detected edges. See [image_skeleton()]
#' @param plot Logical value indicating whether a plot should be created
#' @return A binary version of `image`.
#' @references Sobel, I., and G. Feldman. 1973. A 3×3 isotropic gradient
#'   operator for image processing. Pattern Classification and Scene Analysis:
#'   271–272.
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("sev_leaf_nb.jpg", plot = TRUE)
#' object_edge(img)
#'
object_edge <- function(img,
                        sigma = 1,
                        threshold = "Otsu",
                        thinning = FALSE,
                        plot = TRUE){
  gray <- image_index(img,
                      "GRAY",
                      plot = FALSE,
                      verbose = FALSE)[[1]]
  if(!isFALSE(sigma)){
    gray <- EBImage::gblur(gray, sigma = sigma)
  }
  edata <- sobel_help(gray)

  if (threshold == "Otsu") {
    threshold <- help_otsu(edata@.Data)
  }  else {
    if (is.numeric(threshold)) {
      threshold <- threshold
    }
    else {
      pixels <- raster::raster(t(edata@.Data))
      raster::plot(pixels, col = custom_palette(),  axes = FALSE, asp = NA)
      threshold <- readline("Selected threshold: ")
    }
  }
  edata <- EBImage::Image(edata > threshold)
  if(isTRUE(thinning)){
    edata <- image_thinning(edata, verbose = FALSE, plot = FALSE)
  }
  if(isTRUE(plot)){
    plot(edata)
  }
  return(edata)
}
