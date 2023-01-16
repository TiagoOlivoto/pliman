#' Utilities for Polygons
#'
#' @description
#' Several useful functions to analyze polygons. All of them are
#' based on a set of coordinate points describing the edge of the object(s).
#'
#' * Area measures
#'   - `conv_hull()` Computes convex hull of a set of points.
#'   - `poly_area()` Computes the area of a polygon given by the vertices in the
#'   vectors `x` and `y` using the Shoelace formula, as follows (Lee and Lim,
#'   2017). \deqn{A=\frac{1}{2}\left|\sum_{i=1}^{n}\left(x_{i} y_{i+1}-x_{i+1}
#'   y_{i}\right)\right|}, where `x` and `y` are the coordinates which form the
#'   corners of a polygon, and `n` is the number of coordinates.
#'
#'   - `poly_lw()` Returns the length and width of a polygon based on their
#'   alignment to the x-axis. The length is defined as the range along the
#'   x-axis and the width as the range on the y-axis.
#'
#'   - `poly_mass()` Computes the center of mass of a polygon given by the
#'   vertices in the vectors `x` and `y`.
#'
#'   - `poly_solidity()` Computes the solidity of a shape as the ratio of
#'   the shape area and the convex hull area.
#'
#' * Perimeter measures
#'   - `poly_slide()` Slides the coordinates of a polygon given by the vertices
#'   in the vectors `x` and `y` so that the id-th point become the first one.
#'
#'   - `poly_distpts()` Computes the euclidean distance between every points
#'   of a polygon given by the vertices in the vectors `x` and `y`.
#'
#'   - `poly_centdist()` Computes the euclidean distance between every point of
#'   the perimeter and the centroid of the object.
#'
#'   - `poly_perimeter()` Computes the perimeter of a polygon given by the
#'   vertices in the vectors `x` and `y`.
#'
#'   - `poly_caliper()` Computes the caliper (Also called the Feret's diameter)
#'   of a polygon given by the vertices in the vectors `x` and `y`.
#'
#'
#' * Circularity measures (Montero et al. 2009).
#'    - `poly_circularity()` computes the circularity (C), also called shape
#'    compactness, or roundness measure of an object. It is given by `C = P^2 /
#'    A`, where `P` is the perimeter and `A` is the area of the object.
#'
#'    - `poly_circularity_norm()` computes the normalized circularity (Cn), to
#'    be unity for a circle. This measure is invariant under translation,
#'    rotation, scaling transformations, and dimensionless. It is given by: `Cn
#'    = P^2 / 4*pi*A`.
#'
#'    - `poly_circularity_haralick()` computes the Haralick's circularity (CH).
#'    The method is based on the computation of all the Euclidean distances from
#'    the object centroid to each boundary pixel. With this set of distances,
#'    the media (m) and the standard deviation (sd) are computed. These
#'    statistical parameters are used on a ratio that calculates the
#'    circularity, CH, of a shape, as `CH =  m/sd`
#'
#'    - `poly_convexity()` Computes the convexity of a shape using a ratio
#'    between the perimeter of the convex hull and the perimeter of the polygon.
#'
#'    - `poly_eccentricity()` Computes the eccentricity of a shape using the
#'    ratio of the eigen values (inertia axes of coordinates).
#'
#'    - `poly_elongation()` Computes the elongation of a shape as `1 - width / length`.
#'
#' * Utilities for polygons
#'    - `poly_check()` Checks a set of coordinate points and return a matrix
#'    with `x` and `y` columns.
#'
#'    - `poly_is_closed()` Returns a logical value indicating if a polygon is
#'    closed.
#'
#'    - `poly_close()`, and `poly_unclose()` close and unclose a polygon,
#'    respectively.
#'
#'   - `poly_rotate()` Rotates the polygon coordinates by a `angle` (0-360
#'   degrees) in the trigonometric direction (anti-clockwise).
#'
#'   - `poly_flip_x()`,`poly_flip_y()` flips shapes along the x and y axis,
#'   respectively.
#'
#'   - `poly_align()` Aligns the coordinates along their longer axis using
#'   var-cov matrix and eigen values.
#'
#'   - `poly_sample()` Samples `n` coordinates among existing points. Defaults
#'   to `50`
#'
#'   - `poly_sample_prop()` Samples a proportion of coordinates among existing
#'   points. Defaults to `0.1`
#'
#'   - `poly_spline()` Interpolates a polygon contour.
#'
#'   - `poly_smooth()` Smooths a polygon contour using a simple moving average.
#'
#'   - `poly_jitter()` Add a small amount of noise to a set of point
#'   coordinates. See [base::jitter()] for more details.
#'
#' * `poly_measures()` Is a wrapper around the `poly_*()` functions.
#'
#' @param x,y Coordinate vectors of points. This can be specified as two vectors
#'   (`x` and `y`), or a 2-column matrix `x`. If `x` is a list of vector
#'   coordinates the function will be applied to each element using
#'   [base::lapply()].
#' @param fp The id of the point that will become the new first point. Defaults
#'   to `1`.
#' @param angle The angle (0-360) to rotate the object.
#' @param plot Plots the object? Defaults to `TRUE`
#' @param noise_x,noise_y A numeric factor to define the noise added to `x` and
#'   `y` axis, respectively. See [base::jitter()] for more details.
#' @param vertices The number of spline vertices to create.
#' @param k The number of points to wrap around the ends to obtain a smooth
#'   periodic spline.
#' @param n,prop The number and proportion of coordinates to sample from the
#'   perimeter coordinates. In `poly_smooth()` these arguments can be used to
#'   sample points from the object's perimeter before smoothing.
#' @param niter An integer indicating the number of smoothing iterations.
#' @name utils_polygon
#' @return
#'  * `conv_hull()` and `poly_spline()` returns a matrix with `x` and `y`
#'  coordinates for the convex hull/smooth line in clockwise order. If `x` is a
#'  list, a list of points is returned.
#'
#'  * `poly_area()` returns a `double`, or a numeric vector if `x` is a list of
#'  vector points.
#'
#'  * `poly_mass()` returns a `data.frame` containing the coordinates for the
#'  center of mass, as well as for the maximum and minimum distance from contour
#'  to the center of mass.
#'
#'  * `poly_slides()`, `poly_distpts()`, `poly_spline()`, `poly_close()`,
#'  `poly_unclose()`, `poly_rotate()`, `poly_jitter()`, `poly_sample()`,
#'  `poly_sample_prop()`, and `poly_measures` returns a `data.frame`.
#'
#'  * `poly_perimeter()`, `poly_lw()`, `poly_eccentricity()`,
#'  `poly_convexity()`, `poly_caliper()`, `poly_elongation()`,
#'  `poly_circularity_norm()`, `poly_circularity_haralick()` returns a `double`.
#'
#'
#' @importFrom grDevices chull
#' @importFrom graphics arrows
#' @importFrom stats spline
#' @export
#' @references
#'  Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area of a
#'  Polygon and the Vector Cross Product. The Mathematics Teacher, 110(8),
#'  631–636. \doi{10.5951/mathteacher.110.8.0631}
#'
#'  Montero, R. S., Bribiesca, E., Santiago, R., & Bribiesca, E. (2009). State
#'  of the Art of Compactness and Circularity Measures. International
#'  Mathematical Forum, 4(27), 1305–1335.
#'
#'  Chen, C.H., and P.S.P. Wang. 2005. Handbook of Pattern Recognition and
#'  Computer Vision. 3rd ed. World Scientific.

#'
#' @examples
#' \donttest{
#' library(pliman)
#' # A 2 x 2 square
#' df <- draw_square(side = 2)
#'
#' # square area
#' poly_area(df)
#'
#' # polygon perimeter
#' poly_perimeter(df)
#'
#'# center of mass of the square
#' cm <- poly_mass(df)
#' plot_mass(cm)
#'
#'# The convex hull will be the vertices of the square
#' (conv_square <- conv_hull(df) |> poly_close())
#' plot_contour(conv_square,
#'              col = "blue",
#'              lwd = 6)
#' poly_area(conv_square)
#'
#'
#'################### Example with a polygon ##################
#' x <- c(0, 1,   2, 3,  5, 2, -1, 0, 0)
#' y <- c(5, 6.5, 7, 3,  1, 1,  0, 2, 5)
#' df_poly <- data.frame(x = x, y = y)
#'
#' # area of the polygon
#' plot_polygon(df_poly, fill = "red")
#' poly_area(df_poly)
#'
#' # perimeter of the polygon
#' poly_perimeter(df_poly)
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
#' plot_polygon(conv_poly,
#'              fill = "red",
#'               alpha = 0.2,
#'                add = TRUE)
#'
#'
#' ############ example of circularity measures ################
#' tri <- draw_circle(n = 200, plot = FALSE)
#' plot_polygon(tri, aspect_ratio = 1)
#' poly_circularity_norm(tri)
#'
#' set.seed(1)
#' tri2 <-
#'   draw_circle(n = 200, plot = FALSE) |>
#'   poly_jitter(noise_x = 100, noise_y = 100, plot = FALSE)
#'
#' plot_polygon(tri2, aspect_ratio = 1)
#' poly_circularity_norm(tri2)
#' }

## check a coordinate set of points
poly_check <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_check, y)
  } else {
    if (is.null(y)) {
      if (is.data.frame(x) | is.matrix(x)) {
        vx <- x[, 1]
        vy <- x[, 2]
      }
      if (is.integer(x) &  length(x) == 2) {
        vx <- x[1]
        vy <- x[2]
      }
    } else{
      vx <- x
      vy <- y
    }
    coord <- data.frame(x = vx, y = vy)
    if (!all(sapply(coord, class) %in% c("numeric", 'integer'))){
      stop("Columns must be numeric.")
    }
    colnames(coord) <- c("x", "y")
    return(coord)
  }
}

#' @name utils_polygon
#' @export
# check if a polygon is closed
poly_is_closed <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_is_closed, y)
  } else {
    coord <- poly_check(x, y) |> as.matrix()
    return(identical(coord[1, ], coord[nrow(coord), ]))
  }
}

#' @name utils_polygon
#' @export
# close a polygon, if unclosed
poly_close <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_close, y)
  } else {
    coord <- poly_check(x, y)
    ifelse(poly_is_closed(coord),
           return(coord),
           return(rbind(coord, coord[1, ])))
  }
}
#' @name utils_polygon
#' @export
# unclose a polygon, if closed
poly_unclose <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_unclose, y)
  } else {
    coord <- poly_check(x, y)
    ifelse(poly_is_closed(coord),
           return(coord[-nrow(coord), ]),
           return(coord))
  }
}

#' @name utils_polygon
#' @export
# unclose a polygon, if closed
poly_limits <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_limits, y)
  } else {
    coord <- poly_check(x, y)
    lims <- data.frame(apply(coord, 2, range))
    colnames(lims) <- c("x", "y")
    rownames(lims) <- c("min", "max")
    return(lims)
  }
}

#' @name utils_polygon
#' @export
# computes the convex hull
conv_hull <- function(x, y = NULL){
  if (inherits(x, "list")) {
    lapply(x, conv_hull,  y)
  } else{
    vcts <- poly_check(x, y)
    vec <- chull(vcts)
    return(vcts[vec, ])
  }
}

#' @name utils_polygon
#' @export
# Compute the area of a polygon using shoelace formula
# adapted from https://www.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/
poly_area <- function(x, y = NULL){
  if (inherits(x, "list")) {
    sapply(x, poly_area, y)
  } else{
    coord <- poly_check(x, y) |> poly_close()
    vx <- coord[, 1]
    vy <- coord[, 2]
    nr <- length(vx)
    abs(sum(vx[-nr] * vy[-1] - vx[-1] * vy[-nr]) / 2)
  }
}

#' @name utils_polygon
#' @export
# Slides the coordinates so that the id-th point become the first one.
# adapted from https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_slide <- function(x, y = NULL, fp = 1) {
  if (inherits(x, "list")) {
    lapply(x, poly_slide, y, fp)
  } else{
    coord <- poly_check(x, y)
    n <- nrow(coord)
    slided.rows <- unique(c(fp:n, 1:(fp - 1)))
    return(coord[slided.rows, ])
  }
}

#' @name utils_polygon
#' @export
# calculates the euclidean distance between every points of a shape
# adapted from https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_distpts <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_distpts, y)
  } else{
    coord <- poly_check(x, y)
    n <- nrow(coord)
    d <- sqrt(apply((coord - poly_slide(coord, fp = n))^2, 1, sum))[-1]
    return(d)
  }
}

#' @name utils_polygon
#' @export
# Computes the euclidean distance between every points and the centroid
poly_centdist <- function(x, y = NULL) {
  ed <- function(x, y) {
    return(sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2))
  }
  if (inherits(x, "list")) {
    sapply(x, poly_centdist, y)
  } else{
    coord <- poly_check(x, y)
    return(apply(coord, 1, function(x){
      ed(apply(coord, 2, mean), x)
    }))
  }
}

#' @name utils_polygon
#' @export
# calculates the perimeter
# adapted from https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_perimeter <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_perimeter, y)
  } else{
    coord <- poly_check(x, y)
    return(sum(poly_distpts(coord)))
  }
}

#' @name utils_polygon
#' @export
# rotate
# adapted from https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_rotate <- function(x,
                        y = NULL,
                        angle = 0,
                        plot = TRUE) {
  if (inherits(x, "list")) {
    rot <- lapply(x, poly_rotate, y, angle, plot =  FALSE)
    plot_polygon(rot)
    return(rot)

  } else{
    theta <- angle * pi / 180
    coord <- poly_check(x, y)
    rmat <- matrix(c(cos(-theta), sin(-theta), -sin(-theta), cos(-theta)), nrow = 2)
    rotated <- as.matrix(coord) %*% rmat
    colnames(rotated) <- c("x", "y")
    if(isTRUE(plot)){
      plot_polygon(rotated)
    }
    return(rotated)
  }
}

#' @name utils_polygon
#' @export
# align a polygon
# adapted from https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_align <- function(x,
                       y = NULL,
                       plot = TRUE) {
  if (inherits(x, "list")) {
    coords <- lapply(x, poly_align, y, plot = FALSE)
    if(isTRUE(plot)){
      plot_polygon(coords, merge = FALSE, aspect_ratio = 1)
    }
    return(coords)
  } else{
    coord <- poly_check(x, y)
    aligned <- as.matrix(coord) %*% svd(var(coord))$u
    if(isTRUE(plot)){
      plot_polygon(aligned)
    }
    return(aligned)
  }
}


#' @name utils_polygon
#' @export
# computes width and length
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_lw <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lw <- sapply(x, poly_lw, y) |> t()
    colnames(lw) <- c("length", "width")
    return(lw)
  } else{
    coord <- poly_check(x, y)
    d <- apply(poly_align(coord, plot = FALSE), 2, range)
    lw <- abs(d[2, ] - d[1, ]) |> t()
    colnames(lw) <- c("length", "width")
    return(lw)
  }
}

#' @name utils_polygon
#' @export
# computes the eccentricity
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_eccentricity <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_eccentricity, y)
  } else{
    coord <- poly_check(x, y)
    eig <- eigen(cov(coord))$values
    return(eig[2]/eig[1])
  }
}

#' @name utils_polygon
#' @export
# computes the convexity
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_convexity <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_convexity, y)
  } else{
    coord <- poly_check(x, y)
    return(poly_perimeter(conv_hull(coord))/poly_perimeter(coord))
  }
}

#' @name utils_polygon
#' @export
# computes caliper of a polygon
# https://github.com/MomX/Momocs/blob/356a9b2a59bc48c3a1182021d1240e8f165829f4/R/coo-utilities.R
poly_caliper <- function(x, y = NULL) {
  if (inherits(x, "list")){
    sapply(x, poly_caliper, y)
  } else{
    coord <- poly_check(x, y)
    d <- dist(coord, method = "euclidean")
    return(max(d))
  }
}

#' @name utils_polygon
#' @export
# computes the elongation
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_elongation <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_elongation, y)
  } else{
    coord <- poly_check(x, y)
    lw <- poly_lw(coord)
    return(1 - lw[2]/lw[1])
  }
}

#' @name utils_polygon
#' @export
# computes the solidity
poly_solidity <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_solidity, y)
  } else{
    coord <- poly_check(x, y)
    return(poly_area(coord) / poly_area(conv_hull(coord)))
  }
}

#' @name utils_polygon
#' @export
# flips shapes in along the y axis
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_flip_y <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_flip_y, y)
  } else{
    coord <- poly_check(x, y)
    m <- matrix(c(1, 0, 0, -1), nrow = 2)
    return(as.matrix(coord) %*% m)
  }
}

#' @name utils_polygon
#' @export
# flips shapes in along the x axis
# adapted from https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R#L19
poly_flip_x <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    lapply(x, poly_flip_x, y)
  } else{
    coord <- poly_check(x, y)
    m <- matrix(c(-1, 0, 0, 1), nrow = 2)
    return(as.matrix(coord) %*% m)
  }
}


#' @name utils_polygon
#' @export
# sample n points
poly_sample <- function(x, y = NULL, n = 50) {
  if (inherits(x, "list")) {
    lapply(x, poly_sample, y, n)
  } else{
    coord <- poly_check(x, y)
    if (nrow(coord) < n){
      stop("Less coordinates than n, try coo_interpolate", call. = FALSE)
    }
    sampled <- round(seq(1, nrow(coord), len = n + 1)[-(n + 1)])
    return(coord[sampled, ])
  }
}

#' @name utils_polygon
#' @export
# sample a proportion of points
poly_sample_prop <- function(x, y = NULL, prop = 0.1) {
  if (inherits(x, "list")) {
    lapply(x, poly_sample_prop, y, prop)
  } else{
    coord <- poly_check(x, y)
    return(poly_sample(coord, n = nrow(coord) * prop))
  }
}


#' @name utils_polygon
#' @export
# add an small amount of noise in a set of coordinates
poly_jitter <- function(x,
                        y = NULL,
                        noise_x = 1,
                        noise_y = 1,
                        plot = TRUE) {
  if (inherits(x, "list")) {
    lapply(x, poly_jitter, y, noise_x, noise_y, plot)
  } else{
    coord <- poly_check(x, y)
    coord$x <- jitter(coord$x, factor = noise_x)
    coord$y <- jitter(coord$y, factor = noise_y)
    if(isTRUE(plot)){
      plot_polygon(coord)
    }
    return(coord)
  }
}

#' @name utils_polygon
#' @export
#  calculates the 'circularity measure'. Also called 'compactness' and 'shape factor'
# https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R
poly_circularity <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_circularity, y)
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        coord <- x
      }
    } else{
      coord <- cbind(x = x, y = y)
    }
    return(poly_perimeter(coord)^2 / poly_area(coord))
  }
}

#' @name utils_polygon
#' @export
# calculates calculates Haralick's circularity
# https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R
poly_circularity_haralick <- function(x, y = NULL) {
  ed <- function(x, y) {
    return(sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2))
  }
  coo_centdist <- function(coo) {
    return(
      apply(coo, 1, function(x){
        ed(apply(coo, 2, mean), x)
      })
    )
  }
  if (inherits(x, "list")) {
    sapply(x, poly_circularity_haralick, y)
  } else{
    if(is.null(y)){
      if(is.data.frame(x) | is.matrix(x)){
        coord <- x
      }
    } else{
      coord <- cbind(x = x, y = y)
    }
    cd <- coo_centdist(coord)
    return(mean(cd) / sd(cd))
  }
}

#' @name utils_polygon
#' @export
# calculates calculates Haralick's circularity
# https://github.com/MomX/Momocs/blob/66bd81530ff0c55a46f506208a94ee325c402e0f/R/coo-shapedescriptors.R
poly_circularity_norm <- function(x, y = NULL) {
  if (inherits(x, "list")) {
    sapply(x, poly_circularity_norm, y)
  } else{
    coord <- poly_check(x, y)
    return(poly_perimeter(coord)^2/(poly_area(coord) * 4 * pi))
  }
}


#' @name utils_polygon
#' @export
poly_measures <- function(x, y = NULL){
  if (inherits(x, "list")) {
    valid <- which(sapply(x, function(x){length(as.matrix(x))}) > 2)
    coord <- x[valid]
    res <- do.call(rbind, lapply(coord, poly_measures))
    res$id <- 1:nrow(res)
    shape <- res[, c(ncol(res), 1:ncol(res) -1) ]
    return(shape)
  } else{
    coord <- poly_check(x, y)
    mass <- poly_mass(coord)
    ch <- conv_hull(coord)
    area_ch <- poly_area(ch)
    lw <- poly_lw(coord)
    centdist <- poly_centdist(coord)
    radius_mean = mean_list(centdist)
    radius_min = min_list(centdist)
    radius_max = max_list(centdist)
    radius_sd = sd_list(centdist)
    area <- poly_area(coord)
    shape <- data.frame(x = mass$x,
                        y = mass$y,
                        area = area,
                        area_ch = area_ch,
                        perimeter = poly_perimeter(coord),
                        radius_mean = radius_mean,
                        radius_min = radius_min,
                        radius_max = radius_max,
                        radius_sd = radius_sd,
                        radius_ratio = radius_max / radius_min,
                        diam_mean = radius_mean * 2,
                        diam_min = radius_min * 2,
                        diam_max = radius_max * 2,
                        caliper = poly_caliper(coord),
                        length = lw[1],
                        width = lw[2],
                        solidity = area / area_ch,
                        convexity = poly_convexity(coord),
                        elongation = poly_elongation(coord),
                        circularity = poly_circularity(coord),
                        circularity_haralick = radius_min / radius_sd,
                        circularity_norm = poly_circularity_norm(coord),
                        eccentricity = poly_eccentricity(coord))
  }
  return(shape)
}

#' @name utils_polygon
#' @export
poly_mass <- function(x, y = NULL){
  # adapted from https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
  if (inherits(x, "list") | inherits(x, "iefourier_lst")) {
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
                        k = 2){
  # adapted from https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map
  if (inherits(x, "list")) {
    lapply(x, poly_spline, y, vertices, k)
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
    data.spline <- spline(1:(n+2*k), data[,1], n=vertices)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n+2*k), data[,2], n=vertices)$y
    sl <- cbind(x1, x2)[k < x & x <= n+k, ]
    # ensure a closed line
    rbind(sl, sl[1,])
  }
}


#' @name utils_polygon
#' @export
poly_smooth <- function(x,
                        y = NULL,
                        niter = 10,
                        n = NULL,
                        prop = NULL,
                        plot = TRUE) {
  if (inherits(x, "list")) {
    coords <- lapply(x, poly_smooth, y, niter, n, prop, plot)
    if(isTRUE(plot)){
      plot_polygon(coords, merge = FALSE, aspect_ratio = 1)
    }
    return(coords)
  } else{
    coords <- poly_check(x, y)
    if(!is.null(n) ){
      coords <-
        poly_sample(coords,
                    n = n)
    }
    if(!is.null(prop) ){
      coords <-
        poly_sample_prop(coords,
                         prop = prop)
    }
    p <- nrow(coords)
    a <- 0
    while (a < niter) {
      a <- a + 1
      coo_i <- rbind(coords[-1, ], coords[1, ])
      coo_s <- rbind(coords[p, ], coords[-p, ])
      coords <- coords/2 + coo_i/4 + coo_s/4
    }
    if(isTRUE(plot)){
      plot_polygon(coords, aspect_ratio = 1)
    }
    return(coords)
  }
}

#' Utilities for plotting polygons
#'
#' @description
#'   - `plot_contour()` Plot contour lines.
#'   - `plot_polygon()` Plots a polygon describing the objects.
#'   - `plot_mass()` Plots the center of mass along with maximum and minimum
#'   radius.
#'   - `plot_ellipse()` Plots an ellipse that fits the major and minor axis for
#'   each object.
#'
#' @param x,y Coordinate vectors of points. This can be specified as two vectors
#'   (`x` and `y`), or a 2-column matrix `x`. If `x` is a list of vector
#'   coordinates the function will be applied to each element using
#'   [base::lapply()].
#' @param id The object identification (numeric) to plot the contour/ellipse. By
#'   default (`id = NULL`), the contour is plotted to all objects.
#' @param col,lwd,cex The color, width of the lines, and size of point,
#'   respectively.
#' @param fill,border,alpha The color to fill the polygon, the color of the
#'   polygon's border, and the alpha transparency (1 opaque, 0 transparent).
#' @param random_fill Fill multiple objects with random colors? Defaults to
#'   `TRUE`.
#' @param points Plot the points? Defaults to `FALSE`.
#' @param merge Merge multiple objects into a single plot? Defaults to `TRUE`.
#'   If `FALSE`, a single call `plot()` will be used for each objects. Use
#'   `nrow` and `ncol` to control the number of rows and columns of the window.
#' @param add Add the current plot to a previous one? Defaults to `FALSE`.
#' @param nrow,ncol The number of rows and columns to use in the composite
#'   image. Defaults to `NULL`, i.e., a square grid is produced.
#' @param aspect_ratio The x/y aspect ratio. Defaults to `1`. This will set up
#'   the window so that one data unit in the y direction is equal to one data
#'   unit in the x direction. Set `aspect_ratio = NULL` to fit the object to the
#'   window size.
#' @param show_id Shows the object id? Defaults to `TRUE`.
#' @param xlim,ylim A numeric vector of length 2 (min; max) indicating the range
#'   of `x` and `y`-axes.
#' @param arrow If `TRUE` (default) plots two arrows connecting the center of
#'   mass to the minimum and maximum radius.
#' @param object An object computed with [analyze_objects()].
#' @param ...
#'  * For `plot_contour()` and `plot_ellipse()` further arguments passed on to
#'  [graphics::lines()].
#'  * For `plot_mass()`, further arguments passed on to [graphics::points()].
#'  * For `plot_polygon()`, further arguments passed on to
#'  [graphics::polygon()].
#'
#' @return a `NULL` object.
#' @name utils_polygon_plot
#' @export
#'
#' @examples
#' plot_polygon(contours)
#' plot_contour(contours[[1]], id = 6, col = "red", lwd = 3)
plot_contour <- function(x,
                         y = NULL,
                         id = NULL,
                         col = "black",
                         lwd = 1,
                         ...){
  if (inherits(x, "list") | inherits(x, "iefourier_lst")) {
    if(is.null(id)){
      invisible(lapply(x, plot_contour, y, id, col, lwd))
    } else{
      invisible(lapply(x[id], plot_contour, y, id, col, lwd))
    }
  } else{
    coord <- poly_check(x, y)
    lines(coord,
          col = col,
          xlab = "",
          ylab = "",
          lwd = lwd,
          asp = 1)
  }
}

#' @name utils_polygon_plot
#' @export
plot_polygon <- function(x,
                         y = NULL,
                         fill = "gray",
                         random_fill = TRUE,
                         points = FALSE,
                         merge = TRUE,
                         border = "black",
                         alpha = 1,
                         add = FALSE,
                         nrow = NULL,
                         ncol = NULL,
                         aspect_ratio = 1,
                         show_id = TRUE,
                         xlim = NULL,
                         ylim = NULL,
                         ...){
  if (inherits(x, "list") | inherits(x, "iefourier_lst")) {
    if(isTRUE(merge)){
      if(inherits(x[[1]][[1]], "landmarks_regradi")){
        x <- lapply(x, function(x){
          x[[1]]$coords
        })
      }
      binds <- do.call(rbind,
                       lapply(x, function(x){
                         x
                       }))
      lims <- apply(binds, 2, range)
      xlim <- c(lims[1], lims[2])
      ylim <- c(lims[3], lims[4])

      plot(binds,
           axes = FALSE,
           type = "n",
           xlab = "",
           ylab = "",
           xlim = xlim,
           ylim = ylim,
           asp = 1)
      axis(1)
      axis(2)

      ifelse(random_fill == TRUE,
             colf <- random_color(n = length(x)),
             colf <- rep(fill, length(x)))
      if(isTRUE(show_id)){
        mass <- poly_mass(x)
        xx <- mass[,1]
        xy <- mass[,2]
        labs <- names(x)
      }

      for (i in seq_along(x)) {
        plot_polygon(x[[i]],
                     add = TRUE,
                     fill = colf[[i]],
                     border = border,
                     alpha = alpha,
                     aspect_ratio = aspect_ratio,
                     ...)
        if(isTRUE(show_id)){
          text(xx[i], xy[i], label = labs[i])
        }
      }

    } else{
      if(inherits(x[[1]][[1]], "landmarks_regradi")){
        x <- lapply(x, function(x){
          x[[1]]$coords
        })
      }
      num_plots <- length(x)
      if (is.null(nrow) && is.null(ncol)){
        ncol <- ceiling(sqrt(num_plots))
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
      invisible(lapply(x, plot_polygon, y, fill, random_fill, points, merge, border, alpha, add, aspect_ratio = aspect_ratio, ...))
    }
  } else{
    if(inherits(x, "landmarks_regradi")){
      coords <- x$coords
    } else{
      coords <- poly_check(x, y)
    }
    if(isFALSE(add)){
      plot(coords,
           axes = FALSE,
           type = "n",
           xlab = "",
           ylab = "",
           xlim = xlim,
           ylim = ylim,
           asp = aspect_ratio)
      axis(1)
      axis(2)
    }
    if (!isFALSE(fill)) {
      polygon(coords,
              col = rgb(t(col2rgb(fill) / 255),  alpha = alpha),
              border = border,
              ...)
    } else{
      polygon(coords,
              border = border,
              ...)
    }
    if(isTRUE(points)){
      points(coords, pch = 16)
    }
  }
}


#' @name utils_polygon_plot
#' @export
plot_mass <- function(x,
                      y = NULL,
                      id = NULL,
                      arrow = TRUE,
                      col = "black",
                      cex = 1,
                      lwd = 1){
  if (inherits(x, "list")) {
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

#' @name utils_polygon_plot
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

