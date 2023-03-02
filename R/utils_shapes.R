#' Utilities for drawing coordinates of known shapes
#' @description The functions computes the coordinates of common shapes such as
#'   squares triangles, rectangles and circles.
#' * `draw_circle()` Draws a perfect circle with a desired radius.
#' * `draw_square()` Draws a square with a desired side.
#' * `draw_rectangle()` Draws a rectangle given two desired sides.
#' * `draw_trian_equi()` Draws an equilateral triangle with a desired side.
#' * `draw_trian_rect()` Draws a triangle rectangle given two cathetus.
#' * `draw_n_tagon()` Draws polygons with `n` sides
#' @param radius The radius of the circle. Defaults to `1`.
#' @param n The number of points
#' @param plot Plots the result? Defaults to `TRUE`.
#' @name utils_shapes
#' @return A data frame with the `x` and `y` coordinates
#' @export
#' @examples
#' ########## An example of a circle ##########
#' library(pliman)
#' radius <- 3
#' circ <- draw_circle(radius = radius)
#'
#' # area
#' pi * radius ^ 2
#' poly_area(circ)
#'
#' # perimeter
#' 2 * pi * radius
#' poly_perimeter(circ)
draw_circle <- function(radius = 1,
                        n = 1000,
                        plot = TRUE){
  theta <- seq(0, 2 * pi, length = n)
  circ <-
    data.frame(x = radius * cos(theta),
               y = radius * sin(theta))
  if (isTRUE(plot)) {
    plot(circ,
         type = "l",
         xlab = "",
         ylab = "",
         asp = 1)
    with(circ, polygon(x, y, col = "red"))
    points(x = 0, y = 0, pch = 16)
    arrows(x0 = 0,
           y0 = 0,
           x1 = 0,
           y1 = radius)
    text(x = 0.1 * radius,
         y = 0.5 * radius,
         label = radius)
  }
  return(as.matrix(circ))
}


#' @param side The side of the square/equilateral triangle. Defaults to `2`.
#' @name utils_shapes
#' @export
#'
#' @examples
#'
#' ############ An example of a square ############
#' side <- 2
#' (square <- draw_square(side = side))
#'
#' # area
#' side ^ 2
#' poly_area(square)
#'
#' # perimeter
#' side * 4
#' poly_perimeter(square)
draw_square <- function(side = 2,
                        plot = TRUE){
  x <- c(0, 0, side, side, 0)
  y <- c(0, side, side, 0, 0)
  df <- data.frame(x = x, y = y)
  df <- df[-1, ]
  if (isTRUE(plot)){
    plot(df, type = "n", xlab = "", ylab = "", asp = 1)
    with(df, polygon(x, y, col = "red"))
  }
  return(as.matrix(df))
}


#' @param side1,side2 The first and second sides of the rectangle. Defaults to
#'   `2` and `3`, respectively.
#' @name utils_shapes
#' @export
#'
#' @examples
#'
#' ############ An example of a rectangle ############
#' side1 <- 2
#' side2 <- 3
#' (rect <- draw_rectangle())
#'
#' # area
#' poly_area(rect)
#'
#' # perimeter
#' poly_perimeter(rect)
draw_rectangle <- function(side1 = 2,
                           side2 = 3,
                           plot = TRUE){
  x <- c(0, 0, side2, side2, 0)
  y <- c(0, side1, side1, 0, 0)
  df <- data.frame(x = x, y = y)
  df <- df[-1, ]
  if(isTRUE(plot)){
    plot(df, type = "n", xlab = "", ylab = "", asp = 1)
    with(df, polygon(x, y, col = "red"))
  }
  return(as.matrix(df))
}

#' @name utils_shapes
#' @export
#'
#' @examples
#' ###########  An example of an equilateral triangle #########
#' side <- 1 # defaults
#' (trig <- draw_trian_equi(side = side))
#'
#' ### area (b*h / 2)
#' # height of the triangle
#' (h <- (side * sqrt(3)) / 2)
#' side * h / 2
#'
#' poly_area(trig)
#'
#' ### perimeter (side * 3)
#' poly_perimeter(trig)
#'
draw_trian_equi <- function(side = 2,
                            plot = TRUE){
  angle <- pi/3
  x <- c(0, side, side * cos(angle),0)
  y <- c(0, 0, side * sin(angle), 0)

  df <- data.frame(x = x, y = y)
  df <- df[-1, ]
  if(isTRUE(plot)){
    plot(df, type = "n", xlab = "", ylab = "", asp = 1)
    with(df, polygon(x, y, col = "red"))
  }
  return(as.matrix(df))
}

#' @param cat1,cat2 The first and second cathetus of the right triangle.
#'   Defaults to `1`, and `2`, respectively.
#' @name utils_shapes
#' @export
#'
#' @examples
#' ########### An example of a rectangle triangle ##########
#' cat1 <- 2
#' cat2 <- 3
#' (df <- draw_trian_rect(cat1, cat2))
#' # area
#' (cat1 * cat2) / 2
#' poly_area(df)
#'
#' # perimeter
#' cat1 + cat2 + sqrt(cat1^2 + cat2^2)
#' poly_perimeter(df)
draw_trian_rect <- function(cat1 = 1,
                            cat2 = 2,
                            plot = TRUE){
  x <- c(0, cat2, 0)
  y <- c(0, 0, cat1)
  df <- data.frame(x = x, y = y)
  if(isTRUE(plot)){
    plot(df, type = "n", xlab = "", ylab = "", asp = 1)
    with(df, polygon(x, y, col = "red"))
  }
  return(as.matrix(rbind(df, df[1,])))
}

#' @param n The number of sides in the `n`-tagon.
#' @name utils_shapes
#' @export
#'
#' @examples
#' ############ An creating shapes with n sides ############
#' side <- 2
#' (square <- draw_square(side = side))
#'
#' # area
#' side ^ 2
#' poly_area(square)
#'
#' # perimeter
#' side * 4
#' poly_perimeter(square)
draw_n_tagon <- function(n, plot = TRUE){
  coords <- draw_circle(n = n + 1, plot = FALSE)
  coords <- coords[-1, ]
  if(isTRUE(plot)){
    plot(coords,
         type = "n",
         xlab = "",
         ylab = "",
         asp = 1)
    polygon(coords, col = "red")
  }
  return(as.matrix(coords))
}
