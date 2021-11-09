#' Utilities for Polygons
#' @description
#'* `conv_hull()` Compute convex hull of a set of points.
#'* `poly_area()` Compute the area of a polygon given by the vertices in the
#'vectors `x` and `y`.
#'* `poly_mass()` Compute the center of mass of a polygon given by the vertices
#'in the vectors `x` and `y`.
#'* `poly_spline()` Smooths a polygon contour.
#'* `plot_contour()` Plot contour lines.
#'* `plot_ellipse()` Plots an ellipse that fits the major and minor axis for
#'each object.
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
#' @param vertices The number of spline vertices to create.
#' @param k The number of points to wrap around the ends to obtain a smooth
#'   periodic spline.
#' @param object An object computed with [analyze_objects()].
#' @param closed If `TRUE` (default) returns the vector of points of a closed
#'   polygon, i.e., the first point is replicated as the last one.
#' @param col,lwd,cex The color, width of the lines, and size of point,
#'   respectively.
#' @param id The object identification (numeric) to plot the contour/ellipse. By
#'   default (`id = NULL`), the contour is plotted to all objects
#' @param arrow If `TRUE` (default) plots two arrows connecting the center of
#'   mass to the minimum and maximum radius.
#' @param ...
#' * For `plot_contour()` and `plot_ellipse()` further arguments passed on to
#' [graphics::lines()].
#' * For `plot_mass()`, further arguments passed on to [graphics::points()].
#' @name utils_polygon
#' @return
#' * `conv_hull()` and `poly_spline()` returns a matrix with `x` and `y`
#' coordinates for the convex hull/smooth line in clockwise order. If `x` is a
#' list, a list of points is returned.
#' * `poly_area()` returns a `double`, or a list if `x` is a list of vector
#' points.
#' * `poly_mass()` returns a `data.frame` containing the coordinates for the
#' center of mass, as well as for the maximum and minimum distance from contour
#' to the center of mass.
#' * `plot_contour()`, `plot_mass()`, and `plot_ellipse()` return a `NULL`
#' object.
#' @importFrom grDevices chull
#' @importFrom graphics arrows
#' @importFrom stats spline
#' @export
#' @references Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area
#'   of a Polygon and the Vector Cross Product. The Mathematics Teacher, 110(8),
#'   631–636. \doi{10.5951/mathteacher.110.8.0631}
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
#'# center of mass of the square
#' cm <- poly_mass(df)
#' plot_mass(cm)
#'
#'# The convex hull will be the vertices of the square
#' (conv_square <- conv_hull(df))
#' plot_contour(conv_square,
#'              col = "blue",
#'              lwd = 6)
#' poly_area(conv_square)
#'
#'############# Example with a polygon#############
#' x <- c(0, 1,   2, 3,  5, 2, -1, 0, 0)
#' y <- c(5, 6.5, 7, 3,  1, 1,  0, 2, 5)
#' df_poly <- data.frame(x = x, y = y)
#'
#'# area of the polygon
#' poly_area(df_poly)
#' plot(df_poly, pch = 19, col = "red")
#' with(df_poly, polygon(x, y, col = "red"))
#'
#'# center of mass of polygon
#'# arrows from center of mass to maximum and minimum radius
#' cm <- poly_mass(df_poly)
#' plot_mass(cm, arrow = TRUE, col = "blue")
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
    abs(sum(vx[-nr] * vy[-1] - vx[-1] * vy[-nr]) / 2)
  }
}
#' @name utils_polygon
#' @export
poly_mass <- function(x, y = NULL){
  # adapted from https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
  if(inherits(x, "list")){
    do.call(rbind, lapply(x, poly_mass, y))
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
    if(min(vx) < 0 | min(vy) < 0){
      nr <- length(vx)
      ar <- sum(vx[-nr] * vy[-1] - vx[-1] * vy[-nr])/2
      cx <- 1 / (6 * ar) * sum((vx[-nr] + vx[-1]) * (vx[-nr] * vy[-1] - vx[-1] * vy[-nr]))
      cy <- 1 / (6 * ar) * sum((vy[-nr] + vy[-1]) * (vx[-nr] * vy[-1] - vx[-1] * vy[-nr]))
      ct <- sqrt((vx - cx)^2 + (vy - cy)^2)
    } else{
      nr <- length(vx)
      ar <- abs(sum(vx[-nr] * vy[-1] - vx[-1] * vy[-nr])/2)
      cx <- 1 / (6 * ar) * abs(sum((vx[-nr] + vx[-1]) * (vx[-nr] * vy[-1] - vx[-1] * vy[-nr])))
      cy <- 1 / (6 * ar) * abs(sum((vy[-nr] + vy[-1]) * (vx[-nr] * vy[-1] - vx[-1] * vy[-nr])))
      ct <- sqrt((vx - cx)^2 + (vy - cy)^2)
    }

    return(data.frame(x = cx,
                      y = cy,
                      min_x = vx[which.min(ct)],
                      min_y = vy[which.min(ct)],
                      max_x = vx[which.max(ct)],
                      max_y = vy[which.max(ct)]))
  }
}
#' @name utils_polygon
#' @export
poly_spline <- function(x,
                        y = NULL,
                        vertices = 100,
                        k = 2,
                        ...){
  # adapted from https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map
  if(inherits(x, "list")){
    lapply(x, poly_spline, y, vertices, k, ...)
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        xy <- x
      }
    } else{
      xy <- cbind(x, y)
    }
    n <- dim(xy)[1]
    if (k >= 1) {
      data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
    } else {
      data <- xy
    }
    # Spline the x and y coordinates.
    data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
    sl <- cbind(x1, x2)[k < x & x <= n+k, ]
    # ensure a closed line
    rbind(sl, sl[1,])
  }
}
#' @name utils_polygon
#' @export
plot_contour <- function(x,
                         y = NULL,
                         id = NULL,
                         col = "black",
                         lwd = 1,
                         ...){
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
    lines(vx, vy, col = col, lwd = lwd, ...)
  }
}
#' @name utils_polygon
#' @export
plot_mass <- function(x,
                      y = NULL,
                      id = NULL,
                      arrow = TRUE,
                      col = "black",
                      cex = 1,
                      lwd = 1){
  if(inherits(x, "list")){
    if(is.null(id)){
      invisible(lapply(x, plot_mass, y, id, arrow, col, cex, lwd))
    } else{
      invisible(lapply(x[id], plot_mass, y, id, arrow, col, cex, lwd))
    }
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        vx <- x[, 1]
        vy <- x[, 2]
        if("min_x"  %in% colnames(x)){
          min_x <- x$min_x
          min_y <- x$min_y
          max_x <- x$max_x
          max_y <- x$max_y
        }
      }
    } else{
      vx <- x
      vy <- y
    }
    points(vx, vy, pch = 16, col = col, cex = cex)
    if(arrow == TRUE){
      if(!exists("min_x")){
        stop("Arrows can only be plotted using an object computed with 'poly_mass()`.")
      }
      arrows(vx, vy, min_x, min_y,
             angle = 25,
             col = col,
             lwd = lwd,
             length = 0.15)
      arrows(vx, vy, max_x, max_y,
             angle = 25,
             col = col,
             lwd = lwd,
             length = 0.15)
    }
  }
}
#' @name utils_polygon
#' @export
plot_ellipse <- function(object,
                         id = NULL,
                         col = "black",
                         lwd = 1){
  res <- object$results
  if(is.null(id)){
    ids <- res$id
  } else{
    ids <- id
  }
  points <-
    lapply(1:length(ids), function(i){
      xc <- res$x[i] # center x_c or h
      yc <- res$y[i] # y_c or k
      a <- res$major_axis[i]/2 # major axis length
      b <- res$minor_axis[i]/2 # minor axis length
      tt <- seq(0,  2 * pi, length = 2000)
      sa <- sin(res$theta[i])
      ca <- cos(res$theta[i])
      ct <- cos(tt)
      st <- sin(tt)
      x <- xc + a * ct * ca - b * st * sa
      y <- yc + a * ct * sa + b * st * ca
      lines(x, y, col = col, lwd = lwd, )
    })
}
