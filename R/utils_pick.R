#' Utilities for picking up points in an image
#'
#' * `pick_count()` opens an interactive section where the user will be able to
#' click in the image to count objects (points) manually. In each mouse click, a
#' point is drawn and an upward counter is shown in the console. After `n`
#' counts or after the user press Esc, the interactive process is terminated and
#' the number of counts is returned.
#' * `pick_coord()` Picks coordinates from the image
#' * `pick_palette()`  creates an image palette by picking up color point(s)
#' from the image.
#' * `pick_rgb()` Picks up the RGB values from selected point(s) in the image.
#'
#' @param img An `Image` object.
#' @param n The number of points of the `pick_*` function. Defaults to `Inf`.
#'   This means that picking will run until the user press Esc.
#' @param r The radius of neighborhood pixels. Defaults to `1`.
#' @param col,size The color and size for the marker point.
#' @param shape A character vector indicating the shape of the brush around the
#'   selected pixel. It  can be `"box"`, `"disc"`, `"diamond"`, `"Gaussian"` or
#'   `"line"`. Defaults to `"box"`. In this case, if `'r = 1'`, all the 8
#'   surrounding pixels are sampled. Setting to `"disc"` and increasing the
#'   radius (`r`) will select surrounding pixels towards the format of a sphere
#'   around the selected pixel.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param title The title of the map view when `viewer`is used.
#' @param show How to plot in mapview viewer, either `'rgb` or `'index'`.
#' @param index The index to use for the index view. Defaults to 'B'.
#' @param random Randomize the selected pixels? Defaults to `TRUE`.
#' @param width,height The width and height of the generated palette. Defaults
#'   to `100` for both, i.e., a square image of 100 x 100.
#' @param verbose If `TRUE` (default) shows a counter in the console.
#' @param plot Call a new `plot(img)` before processing? Defaults to `TRUE`.
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
pick_count <- function(img,
                       n = Inf,
                       col = "red",
                       viewer = get_pliman_viewer(),
                       size = 0.8,
                       plot = TRUE,
                       verbose = TRUE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    if(vieweropt == "base"){
      if (isTRUE(plot)) {
        plot(img)
      }
      on.exit(invisible(length(x)))
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
      cat("\n")
      if (i >= n) {
        warning("Maximum number of count achieved. Please, increase the argument `n`.", call. = FALSE)
      }
    } else {
      points <- mv_points(img, title = "Use the 'Draw Marker' tool to pick up points in the plot")
      invisible(nrow(points))
    }
  }
}

#' @name utils_pick
#' @export
pick_coords <- function(img,
                        n = Inf,
                        col = "red",
                        viewer = get_pliman_viewer(),
                        size = 0.8,
                        verbose = TRUE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    pixels <- NULL
    if(vieweropt == "base"){
      plot(img)
      on.exit(invisible(data.frame(x = x, y = y)))
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
        x <- append(x, d$x)
        y <- append(y, d$y)
        points(x, y, type = "p", col = col, cex = size, pch = 19)
        i <- i + 1
      }
      if (i >= n) {
        warning("Maximum number of count achieved. Please, increase the argument `n`.", call. = FALSE)
      }
    } else{
      points <- mv_points(img, title = "Use the 'Draw Marker' tool to pick up points in the plot")
      invisible(points)
    }
  }
}

#' @name utils_pick
#' @export
pick_rgb <- function(img,
                     n = Inf,
                     col = "red",
                     viewer = get_pliman_viewer(),
                     size = 0.8,
                     plot = TRUE,
                     verbose = TRUE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    pixels <- NULL
    if(vieweropt == "base"){
      if (isTRUE(plot)) {
        plot(img)
      }
      on.exit(invisible(pixels))
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
        x <- d$x
        y <- d$y
        pixels <- rbind(pixels, img@.Data[x, y, ])
        points(x, y, type = "p", col = col, cex = size, pch = 19)
        if(isTRUE(verbose)){
          cat("Number of objects:", i, "\r")
        }
        i <- i + 1
      }
      pixels <- data.frame(pixels)

      if (i >= n) {
        warning("Maximum number of count achieved. Please, increase the argument `n`.", call. = FALSE)
      }
    } else{
      points <- mv_points(img, title = "Use the 'Draw Marker' tool to pick up points in the plot")
      pixels <-
        do.call(rbind,
                lapply(1:nrow(points), function(i){
                  img@.Data[points[i, 1], points[i, 2], ]
                })) |>
        as.data.frame()
    }
    colnames(pixels) <- c("R", "G", "B")
    pixels$id <- 1:nrow(pixels)
    pixels <- pixels[, c("id", "R", "G", "B")]
  }
  invisible(pixels)
}

#' @name utils_pick
#' @export
pick_palette <- function(img,
                         n = Inf,
                         r = 1,
                         shape = "box",
                         viewer = get_pliman_viewer(),
                         show = "rgb",
                         title = "Pick colors in the image",
                         index =  "B",
                         random = TRUE,
                         width = 100,
                         height = 100,
                         col = "red",
                         size = 0.8,
                         plot = TRUE,
                         palette = TRUE,
                         verbose = TRUE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    if(vieweropt == "base"){

      if (isTRUE(plot)) {
        plot(img)
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
        R <- img@.Data[xrmin:xrmax, yrmin:yrmax, 1][kern]
        G <- img@.Data[xrmin:xrmax, yrmin:yrmax, 2][kern]
        B <- img@.Data[xrmin:xrmax, yrmin:yrmax, 3][kern]
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
      if(i > 1){
        on.exit(invisible(pal))
      }
    } else{
      mvpoin <- mv_points(img, show = show, title = title, index = index)
      bind <- NULL
      for (i in 1:nrow(mvpoin)) {
        xrmin <- trunc(mvpoin[, 1][i]) - r
        xrmax <- trunc(mvpoin[, 1][i]) + r
        yrmin <- trunc(mvpoin[, 2][i]) - r
        yrmax <- trunc(mvpoin[, 2][i]) + r
        sqr <- xrmax - xrmin + 1
        kern <- as.logical(EBImage::makeBrush(sqr, shape = shape))
        R <- img@.Data[xrmin:xrmax, yrmin:yrmax, 1][kern]
        G <- img@.Data[xrmin:xrmax, yrmin:yrmax, 2][kern]
        B <- img@.Data[xrmin:xrmax, yrmin:yrmax, 3][kern]
        bind <- rbind(bind, cbind(R, G, B))
      }
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
    return(pal)
  }
}
