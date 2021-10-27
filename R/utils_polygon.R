#' Utilities for Polygons
#' @description
#'* `conv_hull()` Compute convex hull of a set of points.
#'* `poly_area()` Compute the area of a polygon given by the vertices in the
#'vectors `x` and `y`.
#'* `plot_contour()` Plot contour lines using [graphics::lines()].
#'
#' @details
#' `poly_area()` computes the area of a polygon given a set of x and y
#' coordinates using the Shoelace formula, as follows (Lee and Lim, 2017).
#' \deqn{A=\frac{1}{2}\left|\sum_{i=1}^{n}\left(x_{i} y_{i+1}-x_{i+1}
#' y_{i}\right)\right|}, where `x` and `y` are the coordinates which form the
#' corners of a polygon, and `n` is the number of coordinates.
#' @param x,y Coordinate vectors of points. This can be specified as two vectors
#'   (`x` and `y`), or a 2-column matrix `x`. If `x` is a list of vector
#'   coordinates the function will be applied to each element using
#'   [base::lapply()].
#' @param closed If `TRUE` (default) returns the vector of points of a closed
#'   polygon, i.e., the first point is replicated as the last one.
#' @param col,lwd The color and width of the lines, respectively.
#' @param id The object identification (numeric) to plot the contour. By default
#'   (`id = NULL`), the contour is plotted to all objects
#' @name utils_polygon
#' @return
#' * `conv_hull()` returns a matrix with `x` and `y` coordinates for the convex
#' hull in clockwise order. If `x` is a list, a list of points is returned.
#' * `poly_area()` returns a `double`, or a list if `x` is a list of vector
#' points.
#' * `plot_contour()` returns a `double`, or a list if `x` is a list of vector
#' points.
#' @importFrom grDevices chull
#' @export
#' @references Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area
#'   of a Polygon and the Vector Cross Product. The Mathematics Teacher, 110(8),
#'   631–636. doi:10.5951/MATHTEACHER.110.8.0631
#' @examples
#' \donttest{
#' library(pliman)
#' # A 2 x 2 square
#' x <- c(0, 0, 2, 2, 0)
#' y <- c(0, 2, 2, 0, 0)
#' df <- data.frame(x = x, y = y)
#' plot(df)
#' with(df, polygon(x, y, col = "red"))
#'
#' poly_area(x, y)
#' poly_area(df)
#'
#'# The convex hull will be the vertices of the square
#' (conv_square <- conv_hull(df))
#' poly_area(conv_square)
#'
#' with(conv_square,
#'      lines(x, y,
#'            col  = "blue",
#'            lwd = 4))
#'
#'# polygon
#' x <- c(0, 2, 2,  5, 2, -1, 0, 0)
#' y <- c(5, 6, 3,  1, 1,  0, 2, 5)
#' df_poly <- data.frame(x = x, y = y)
#'
#'# area of the polygon
#' poly_area(df_poly)
#'
#' plot(df_poly, pch = 19, col = "red")
#' with(df_poly, polygon(x, y, col = "red"))
#'
#'
#'# vertices of the convex hull
#' (conv_poly <- conv_hull(df_poly))
#'
#'# area of the convex hull
#' poly_area(conv_poly)
#'
#' with(conv_poly,
#'      polygon(x, y,
#'            col  = rgb(1, 0, 0, 0.2)))
#' }
conv_hull <- function(x, y = NULL, closed = TRUE){
  if(inherits(x, "list")){
    lapply(x, conv_hull,  y, closed)
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        vx <- x[, 1]
        vy <- x[, 2]
      }
    } else{
      vx <- x
      vy <- y
    }
    vcts <- data.frame(x = vx, y = vy)
    vec <- chull(vcts)
    if(isTRUE(closed)){
      vcts <- vcts[c(vec, vec[1]), ]
    } else{
      vcts <- vcts[vec, ]
    }
    return(vcts)
  }
}
#' @name utils_polygon
#' @export
# Compute the area of a polygon using shoelace formula
# adapted from https://www.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/
poly_area <- function(x, y = NULL){
  if(inherits(x, "list")){
    lapply(x, poly_area, y)
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        vx <- x[, 1]
        vy <- x[, 2]
      }
    } else{
      vx <- x
      vy <- y
    }
    nr <- length(vx)
    trunc(abs(sum(vx[-nr] * vy[-1] - vx[-1] * vy[-nr]) / 2))
  }
}
#' @name utils_polygon
#' @export
plot_contour <- function(x,
                         y = NULL,
                         id = NULL,
                         col = "white",
                         lwd = 1){
  if(inherits(x, "list")){
    if(is.null(id)){
      invisible(lapply(x, plot_contour, y, id, col, lwd))
    } else{
      invisible(lapply(x[id], plot_contour, y, id, col, lwd))
    }
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        vx <- x[, 1]
        vy <- x[, 2]
      }
    } else{
      vx <- x
      vy <- y
    }
    lines(vx, vy, col = col, lwd = lwd)
  }
}
