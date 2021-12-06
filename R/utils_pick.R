#' Picking up points in an image
#'
#' * `pick_count()` opens an interactive section where the user will be able to
#' click in the image to count objects (points) manually. In each mouse click, a
#' point is drawn and an upward counter is shown in the console. After `n`
#' counts or after the user press Esc, the interactive process is terminated and
#' a data.frame with the `x` and `y` coordinates for each point is returned.
#' * `pick_palette()`  creates an image palette by picking up color point(s)
#' from the image.
#' * `pick_rgb()` Picks up the RGB values from selected point(s) in the image.
#' @param image An `Image` object.
#' @param n The number of points of the `pick_*` function. Defaults to `Inf`.
#'   This means that picking will run until the user press Esc.
#' @param r The radius of neighborhood pixels. Defaults to `3`.
#' @param col,size The color and size for the marker point.
#' @param shape A character vector indicating the shape of the brush around the
#'   selected pixel. It  can be `"box"`, `"disc"`, `"diamond"`, `"Gaussian"` or
#'   `"line"`. Defaults to `"box"`. In this case, if `'r = 1'`, all the 8
#'   surrounding pixels are sampled. Setting to `"disc"` and increasing the
#'   radius (`r`) will select surrounding pixels towards the format of a sphere
#'   around the selected pixel.
#' @param random Randomize the selected pixels? Defaults to `TRUE`.
#' @param width,height The width and height of the generated palette. Defaults
#'   to `100` for both, i.e., a square image of 100 x 100.
#' @param verbose If `TRUE` (default) shows a counter in the console.
#' @param plot Call a new `plot(image)` before processing? Defaults to `TRUE`.
#' @param palette Plot the generated palette? Defaults to `TRUE`.
#' @return
#' * `pick_count()` returns `data.frame` with the `x` and `y` coordinates of the
#' selected point(x).
#' * `pick_rgb()` returns a `data.frame` with the R, G, and B values of the
#' selected point(s).
#' * `pick_palette()` returns an object of class `Image`.
#' @name utils_pick
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@gmail.com}
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#'
#' # start a counting process
#' pick_count(img)
#'
#' # get rgb from point(s)
#' pick_rgb(img)
#'
#' # create a palette from point(s)
#' pick_palette(img)
#' }
pick_count <- function(image,
                       n = Inf,
                       col = "red",
                       size = 0.8,
                       plot = TRUE,
                       verbose = TRUE){
  if (isTRUE(interactive())) {
    if (isTRUE(plot)) {
      plot(image)
    }
    on.exit(return(data.frame(x = x, y = y)))
    if(isTRUE(verbose)){
      message("Use the first mouse button to pick up points in the plot.\nPress Esc to exit.")
    }
    x <- y <- NULL
    i <- 1
    while (i <= n) {
      d <- locator(n = 1)
      if (is.null(d)) {
        break
      }
      x <- c(x, d$x)
      y <- c(y, d$y)
      points(x, y, type = "p", col = col, cex = size, pch = 19)
      cat("Number of objects:", i, "\r")
      i <- i + 1
    }
    if (i >= n) {
      warning("Maximum number of count achieved. Please, increase the argument `n`.", call. = FALSE)
    }
  }
}

#' @name utils_pick
#' @export
pick_rgb <- function(image,
                     n = Inf,
                     col = "red",
                     size = 0.8,
                     plot = TRUE,
                     verbose = TRUE){
  if (isTRUE(interactive())) {
    if (isTRUE(plot)) {
      plot(image)
    }
    on.exit(return(pixels))
    if(isTRUE(verbose)){
    message("Use the first mouse button to pick up points in the plot.\nPress Esc to exit.")
    }
    pixels <- NULL
    x <- y <- NULL
    i <- 1
    while (i <= n) {
      d <- locator(n = 1)
      if (is.null(d)) {
        break
      }
      x <- d$x
      y <- d$y
      pixels <- rbind(pixels, image@.Data[x, y, ])
      points(x, y, type = "p", col = col, cex = size, pch = 19)
      if(isTRUE(verbose)){
        cat("Number of objects:", i, "\r")
      }
      i <- i + 1
    }
    pixels <- data.frame(pixels)
    colnames(pixels) <- c("R", "G", "B")
    pixels$id <- 1:nrow(pixels)
    pixels <- pixels[, c("id", "R", "G", "B")]
    if (i >= n) {
      warning("Maximum number of count achieved. Please, increase the argument `n`.", call. = FALSE)
    }
  }
}

#' @name utils_pick
#' @export
pick_palette <- function(image,
                         n = Inf,
                         r = 3,
                         shape = "box",
                         random = TRUE,
                         width = 100,
                         height = 100,
                         col = "red",
                         size = 0.8,
                         plot = TRUE,
                         palette = TRUE,
                         verbose = TRUE){
  if (isTRUE(interactive())) {
    if (isTRUE(plot)) {
      plot(image)
    }
    if(isTRUE(verbose)){
      message("Use the first mouse button to pick up points in the plot.\nPress Esc to exit.")
    }
    bind <- NULL
    i <- 1
    while (i <= n) {
      d <- locator(n = 1)
      if (is.null(d)) {
        break
      }
      xrmin <- trunc(d$x) - r
      xrmax <- trunc(d$x) + r
      yrmin <- trunc(d$y) - r
      yrmax <- trunc(d$y) + r
      sqr <- xrmax - xrmin + 1
      kern <- as.logical(EBImage::makeBrush(sqr, shape = shape))
      R <- image@.Data[xrmin:xrmax, yrmin:yrmax, 1][kern]
      G <- image@.Data[xrmin:xrmax, yrmin:yrmax, 2][kern]
      B <- image@.Data[xrmin:xrmax, yrmin:yrmax, 3][kern]
      rect(xrmin, yrmin, xrmax, yrmax, border = col, lwd = size)
      bind <- rbind(bind, cbind(R, G, B))
      if(isTRUE(verbose)){
        cat("Number of points:", i, "\r")
      }
      i <- i + 1
    }
    if(i == 1){
      stop("Process interrupted", call. = FALSE)
    }
    dim_mat <- trunc(sqrt(nrow(bind)))
    if(isTRUE(random)){
      bind <- bind[sample(1:nrow(bind)), ][1:dim_mat^2, ]
    } else{
      bind <- bind[1:dim_mat^2, ]
    }
    pal <-
      EBImage::Image(c(bind[, 1], bind[, 2], bind[, 3]),
                     dim = c(dim_mat, dim_mat, 3),
                     colormode = "Color") %>%
      image_resize(width = width, height = height)
    if(isTRUE(palette)){
      plot(pal)
    }
    if(i > 1){
    on.exit(return(pal))
    }
  }
}
