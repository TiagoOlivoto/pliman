#' Utilities for Polygons
#'
#' @description
#' Several useful functions for analyzing polygons. All of them are based on a
#' set of coordinate points that describe the edge of the object(s). If a list
#' of polygons is provided, it loops through the list and computes what is
#' needed for each element of the list.
#'
#'* Polygon measures
#'
#'    - `conv_hull()` Computes the convex hull of a set of points.
#'
#'    - `conv_hull_unified()` Computes the convex hull of a set of points. Compared
#' to `conv_hull()`, `conv_hull_unified()` binds (unifies) the coordinates when
#' x is a list of coordinates.
#'
#'    - `poly_area()` Computes the area of a polygon given by the vertices in the
#' vectors x and y using the Shoelace formula, as follows (Lee and Lim,
#' 2017): \deqn{A=\frac{1}{2}\left|\sum_{i=1}^{n}\left(x_{i} y_{i+1}-x_{i+1}
#' y_{i}\right)\right|} where `x` and `y` are the coordinates that form the
#' corners of a polygon, and `n` is the number of coordinates.
#'
#'    - `poly_angles()` Calculates the internal angles of the polygon using the
#' law of cosines.
#'
#'
#'    - `poly_lw()` Returns the length and width of a polygon based on its
#' alignment to the y-axis (with poly_align()). The length is defined as the
#' range along the x-axis, and the width is defined as the range on the y-axis.
#'
#'    - `poly_mass()` Computes the center of mass of a polygon given by the
#' vertices in the vectors in `x`.
#'
#'    - `poly_solidity()` Computes the solidity of a shape as the ratio of
#' the shape area and the convex hull area.
#'
#' * Perimeter measures
#'    - `poly_slide()` Slides the coordinates of a polygon given by the vertices
#' in the vectors x and y so that the id-th point becomes the first one.
#'
#'    - `poly_distpts()` Computes the Euclidean distance between every point
#' of a polygon given by the vertices in the vectors x and y.
#'
#'    - `poly_centdist()` Computes the Euclidean distance between every point on
#' the perimeter and the centroid of the object.
#'
#'    - `poly_perimeter()` Computes the perimeter of a polygon given by the
#' vertices in the vectors x and y.
#'
#'    - `poly_caliper()` Computes the caliper (also called the Feret's diameter)
#' of a polygon given by the vertices in the vectors x and y.
#'
#'
#'  * Circularity measures (Montero et al. 2009).
#'
#'    - `poly_circularity()` computes the circularity (C), also called shape
#' compactness or roundness measure, of an object. It is given by `C = P^2 / A`,
#'  where P is the perimeter and A is the area of the object.
#'
#'    - `poly_circularity_norm()` computes the normalized circularity (Cn), which
#' is unity for a circle. This measure is invariant under translation,
#' rotation, scaling transformations, and is dimensionless. It is given by:
#' `Cn = P^2 / 4*pi*A`.
#'
#'    - `poly_circularity_haralick()` computes Haralick's circularity (CH). The
#' method is based on computing all Euclidean distances from the object
#' centroid to each boundary pixel. With this set of distances, the mean (`m`)
#' and the standard deviation (`sd`) are computed. These statistical parameters
#' are used to calculate the circularity, CH, of a shape as `CH = m/sd`.
#'
#'    - `poly_convexity()` computes the convexity of a shape using the ratio
#' between the perimeter of the convex hull and the perimeter of the polygon.
#'
#'    - `poly_eccentricity()` computes the eccentricity of a shape using the
#' ratio of the eigenvalues (inertia axes of coordinates).
#'
#'    - `poly_elongation()` computes the elongation of a shape as `1 - width / length`.
#'
#'
#' * Utilities for polygons
#'
#'    - `poly_check()` Checks a set of coordinate points and returns a matrix
#' with x and y columns.
#'
#'    - `poly_is_closed()` Returns a logical value indicating if a polygon is
#' closed.
#'
#'    - `poly_close()` and `poly_unclose()` close and unclose a polygon,
#' respectively.
#'
#'    - `poly_rotate() `Rotates the polygon coordinates by an angle (0-360
#' degrees) in the counterclockwise direction.
#'
#'    - `poly_flip_x()`, `poly_flip_y()` flip shapes along the x-axis and y-axis,
#' respectively.
#'
#'    - `poly_align()` Aligns the coordinates along their longer axis using the
#' var-cov matrix and eigen values.
#'
#'    - `poly_center()` Centers the coordinates on the origin.
#'
#'    - `poly_sample()` Samples n coordinates from existing points. Defaults
#' to 50.
#'
#'    - `poly_sample_prop()` Samples a proportion of coordinates from existing
#' points. Defaults to 0.1.
#'
#'    - `poly_spline()` Interpolates the polygon contour.
#'
#'    - `poly_smooth()` Smooths the polygon contour using a simple moving average.
#'
#'    - `poly_jitter()` Adds a small amount of noise to a set of point
#' coordinates. See [base::jitter()] for more details.
#'
#'
#' * `poly_measures()` Is a wrapper around the `poly_*()` functions.
#'
#' @param x A 2-column matrix with the `x` and `y` coordinates. If `x` is a list
#' of vector coordinates, the function will be applied to each element using
#' [base::lapply()] or [base::sapply()].
#' @param fp The ID of the point that will become the new first point. Defaults
#' to 1.
#' @param angle The angle (0-360) to rotate the object.
#' @param plot Should the object be plotted? Defaults to `TRUE`.
#' @param noise_x,noise_y A numeric factor to define the noise added to the `x`
#' and `y` axes, respectively. See [base::jitter()] for more details.
#' @param vertices The number of spline vertices to create.
#' @param k The number of points to wrap around the ends to obtain a smooth
#' periodic spline.
#' @param n,prop The number and proportion of coordinates to sample from the
#' perimeter coordinates. In` poly_smooth()`, these arguments can be used to
#' sample points from the object's perimeter before smoothing.
#' @param niter An integer indicating the number of smoothing iterations.
#'
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
#' df_poly <- cbind(x, y)
#'
#' # area of the polygon
#' plot_polygon(df_poly, fill = "red")
#' poly_area(df_poly)
#'
#' # perimeter of the polygon
#' poly_perimeter(df_poly)
#'
#'# center of mass of polygon
#' cm <- poly_mass(df_poly)
#' plot_mass(cm,  col = "blue")
#'
#'# vertices of the convex hull
#' (conv_poly <- conv_hull(df_poly))
#'
#'# area of the convex hull
#' poly_area(conv_poly)
#'
#' plot_polygon(conv_poly,
#'              fill = "red",
#'              alpha = 0.2,
#'              add = TRUE)
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
poly_check <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, poly_check)
  } else{
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    return(x)
  }
}

#' @name utils_polygon
#' @export
# check if a polygon is closed
poly_is_closed <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, poly_is_closed)
  } else {
    coord <- poly_check(x)
    return(identical(coord[1, ], coord[nrow(coord), ]))
  }
}

#' @name utils_polygon
#' @export
# close a polygon, if unclosed
poly_close <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, poly_close)
  } else {
    coord <- poly_check(x)
    ifelse(poly_is_closed(coord),
           return(coord),
           return(rbind(coord, coord[1, ])))
  }
}

#' @name utils_polygon
#' @export
# unclose a polygon, if closed
poly_unclose <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, poly_unclose)
  } else {
    coord <- poly_check(x)
    ifelse(poly_is_closed(coord),
           return(coord[-nrow(coord), ]),
           return(coord))
  }
}

#' @name utils_polygon
#' @export
# compute the internal angles of a polygon
poly_angles <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, function(x){
      help_poly_angles(poly_unclose(x))
    })
  } else {
    coord <- poly_unclose(x)
    help_poly_angles(coord)
  }
}

#' @name utils_polygon
#' @export
poly_limits <- function(x) {
  if (inherits(x, "list")) {
    res <- do.call(rbind, lapply(x, help_limits))
    colnames(res) <- c("xmin", "xmax", "ymin", "ymax")
    return(res)
  } else {
    res <- t(as.data.frame(help_limits(x)))
    colnames(res) <- c("xmin", "xmax", "ymin", "ymax")
    rownames(res) <- NULL
    return(res)
  }
}

#' @name utils_polygon
#' @export
# computes the convex hull
conv_hull <- function(x){
  if (inherits(x, "list")) {
    lapply(x, conv_hull)
  } else{
    vcts <- poly_check(x)
    vec <- chull(vcts)
    return(vcts[vec, ])
  }
}

#' @name utils_polygon
#' @export
conv_hull_unified <- function(x){
  if(inherits(x, "list")){
    unified <- do.call(rbind, lapply(x, function(x){x}))
  } else{
    unified <- x
  }
  ch <- unified[chull(unified), ]
  ch <- rbind(ch,  ch[1, ])
  return(ch)
}


#' @name utils_polygon
#' @export
# Compute the area of a polygon using shoelace formula
# adapted from https://www.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/
poly_area <- function(x){
  help_area(x)
}


#' @name utils_polygon
#' @export
# Slides the coordinates so that the id-th point become the first one.
poly_slide <- function(x, fp = 1){
  if (inherits(x, "list")) {
    sapply(x, help_slide, fp)
  } else{
    help_slide(x, fp = fp)
  }
}




#' @name utils_polygon
#' @export
# calculates the euclidean distance between every points of a shape
poly_distpts <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, help_distpts)
  } else{
    help_distpts(x)
  }
}

#' @name utils_polygon
#' @export
# Computes the euclidean distance between every points and the centroid
poly_centdist <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, help_centdist)
  } else{
    help_centdist(x)
  }
}



#' @name utils_polygon
#' @export
# calculates the perimeter
poly_perimeter <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){sum(help_distpts(x))})
  } else{
    return(sum(help_distpts(x)))
  }
}


#' @name utils_polygon
#' @export
# rotate
poly_rotate <- function(x, angle, plot = TRUE) {
  if (inherits(x, "list")) {
    rotated <- lapply(x, help_rotate, angle)
    if(isTRUE(plot)){
      plot_polygon(rotated)
    }
    return(rotated)
  } else{
    rotated <- help_rotate(x, angle)
    if(isTRUE(plot)){
      plot_polygon(rotated)
    }
    return(rotated)
  }
}

#' @name utils_polygon
#' @export
poly_align <- function(x, plot = TRUE) {
  if (inherits(x, "list")) {
    aligned <- lapply(x, help_align)
    if(isTRUE(plot)){
      plot_polygon(aligned, merge = FALSE, aspect_ratio = 1)
    }
    return(aligned)
  } else{
    aligned <- help_align(x)
    if(isTRUE(plot)){
      plot_polygon(aligned, merge = FALSE, aspect_ratio = 1)
    }
    return(aligned)
  }
}

#' @name utils_polygon
#' @export
poly_center <- function(x, plot = TRUE) {
  if (inherits(x, "list")) {
    centered <- lapply(x, poly_center)
    if(isTRUE(plot)){
      plot_polygon(centered, merge = FALSE, aspect_ratio = 1)
    }
    return(centered)
  } else{
    centered <-
      data.frame(X1 = x[,1] - mean(x[,1]),
                 X2 = x[,2] - mean(x[,2]))
    if(isTRUE(plot)){
      plot_polygon(centered, merge = FALSE, aspect_ratio = 1)
    }
    return(centered)
  }
}

#' @name utils_polygon
#' @export
# computes width and length
poly_lw <- function(x) {
  lw <- help_lw(x)
  colnames(lw) <- c("length", "width")
  return(lw)
}

#' @name utils_polygon
#' @export
poly_eccentricity <- function(x) {
  help_eigen_ratio(x)
}


#' @name utils_polygon
#' @export
# computes the convexity
poly_convexity <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){sum(help_distpts(x[chull(x), ])) / sum(help_distpts(x))})
  } else{
    return(sum(help_distpts(x[chull(x), ])) / sum(help_distpts(x)))
  }
}

#' @name utils_polygon
#' @export
poly_caliper <- function(x) {
  if (inherits(x, "list")){
    sapply(x, poly_caliper)
  } else{
    if(nrow(x) > 250){
      x <- landmarks_regradi(x, n = 250, plot = FALSE)$coord
    }
    return(max(dist(x)))
  }
}



#' @name utils_polygon
#' @export
# computes the elongation
poly_elongation <- function(x) {
  help_elongation(x)
}


#' @name utils_polygon
#' @export
# computes the solidity
poly_solidity <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){
      help_area(x) / help_area(x[chull(x), ])
    })
  } else{
    help_area(x) / help_area(x[chull(x), ])
  }
}


#' @name utils_polygon
#' @export
# flips shapes in along the y axis
poly_flip_y <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, help_flip_y)
  } else{
    help_flip_y(x)
  }
}



#' @name utils_polygon
#' @export
# flips shapes in along the x axis
poly_flip_x <- function(x) {
  if (inherits(x, "list")) {
    lapply(x, help_flip_x)
  } else{
    help_flip_x(x)
  }
}


#' @name utils_polygon
#' @export
# sample n points
poly_sample <- function(x, n = 50) {
  if (inherits(x, "list")) {
    lapply(x, poly_sample, n)
  } else{
    coord <- poly_check(x)
    if (nrow(coord) < n){
      stop("Less coordinates than n", call. = FALSE)
    }
    sampled <- round(seq(1, nrow(coord), len = n + 1)[-(n + 1)])
    return(coord[sampled, ])
  }
}

#' @name utils_polygon
#' @export
# sample a proportion of points
poly_sample_prop <- function(x, prop = 0.1) {
  if (inherits(x, "list")) {
    lapply(x, poly_sample_prop, prop)
  } else{
    coord <- poly_check(x)
    return(poly_sample(coord, n = nrow(coord) * prop))
  }
}


#' @name utils_polygon
#' @export
# add an small amount of noise in a set of coordinates
poly_jitter <- function(x,
                        noise_x = 1,
                        noise_y = 1,
                        plot = TRUE) {
  if (inherits(x, "list")) {
    lapply(x, poly_jitter, noise_x, noise_y, plot)
  } else{
    coord <- poly_check(x)
    coord[,1] <- jitter(coord[,1], factor = noise_x)
    coord[,2] <- jitter(coord[,2], factor = noise_y)
    if(isTRUE(plot)){
      plot_polygon(coord)
    }
    return(coord)
  }
}

#' @name utils_polygon
#' @export
#  calculates the 'circularity measure'. Also called 'compactness' and 'shape factor'
poly_circularity <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){
      sum(help_distpts(x))^2 / help_area(x)
    })
  } else{
    sum(help_distpts(x))^2 / help_area(x)
  }
}

#' @name utils_polygon
#' @export
# calculates calculates Haralick's circularity
poly_circularity_norm <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){
      (help_area(x) * 4 * pi) / sum(help_distpts(x))^2
    })
  } else{
    (help_area(x) * 4 * pi) / sum(help_distpts(x))^2
  }
}


#' @name utils_polygon
#' @export
# calculates calculates Haralick's circularity
poly_circularity_haralick <- function(x) {
  if (inherits(x, "list")) {
    sapply(x, function(x){
      cd <- help_centdist(x)
      mean(cd) / sd(cd)
    })
  } else{
    cd <- help_centdist(x)
    mean(cd) / sd(cd)
  }
}



#' @name utils_polygon
#' @export
poly_mass <- function(x){
  if (inherits(x, "list") | inherits(x, "iefourier_lst")) {
    cm <-
      do.call(rbind,
              lapply(x, function(x){
                cm <- t(help_mc(x))
              })
      )
    colnames(cm) <- c("x", "y")
    return(cm)
  } else{
    cm <- t(help_mc(x))
    colnames(cm) <- c("x", "y")
    return(cm)
  }
}


#' @name utils_polygon
#' @export
poly_spline <- function(x,
                        vertices = 100,
                        k = 2){
  # adapted from https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map
  if (inherits(x, "list")) {
    lapply(x, poly_spline, vertices, k)
  } else{
    if(is.data.frame(x) | is.matrix(x)){
      xy <- x
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
                        niter = 10,
                        n = NULL,
                        prop = NULL,
                        plot = TRUE) {
  if (inherits(x, "list")) {
    coords <- lapply(x, poly_smooth,niter, n, prop, plot)
    if(isTRUE(plot)){
      plot_polygon(coords, merge = FALSE, aspect_ratio = 1)
    }
    return(coords)
  } else{
    coords <- poly_check(x)
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
    coords <- help_smoth(coords, niter)
    if(isTRUE(plot)){
      plot_polygon(coords, aspect_ratio = 1)
    }
    return(coords)
  }
}

#' @name utils_polygon
#' @export
poly_measures <- function(x){
  if (inherits(x, "list")) {
    valid <- which(sapply(x, function(x){length(as.matrix(x))}) > 2)
    coord <- x[valid]
    res <- do.call(rbind, lapply(coord, poly_measures))
    res$id <- 1:nrow(res)
    shape <- res[, c(ncol(res), 1:ncol(res) -1) ]
    return(shape)
  } else{
    coord <- poly_check(x)
    mass <- poly_mass(coord)
    ch <- coord[chull(coord), ]
    area_ch <- help_area(ch)
    lw <- help_lw(coord)
    centdist <- help_centdist(coord)
    radius_mean = mean(centdist)
    radius_min = min(centdist)
    radius_max = max(centdist)
    radius_sd = sd(centdist)
    area <- help_area(coord)
    shape <- data.frame(x = mass[[1]],
                        y = mass[[2]],
                        area = area,
                        area_ch = area_ch,
                        perimeter = sum(help_distpts(coord)),
                        radius_mean = radius_mean,
                        radius_min = radius_min,
                        radius_max = radius_max,
                        radius_sd = radius_sd,
                        radius_ratio = radius_max / radius_min,
                        diam_mean = radius_mean * 2,
                        diam_min = radius_min * 2,
                        diam_max = radius_max * 2,
                        caliper = lw[1],
                        length = lw[1],
                        width = lw[2],
                        solidity = area / area_ch,
                        convexity = poly_convexity(coord),
                        elongation = poly_elongation(coord),
                        circularity = poly_circularity(coord),
                        circularity_haralick = poly_circularity_haralick(coord),
                        circularity_norm = poly_circularity_norm(coord),
                        eccentricity = poly_eccentricity(coord),
                        pcv = poly_pcv(coord))
  }
  return(shape)
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
#' @param x A 2-column matrix with the `x` and `y` coordinates.
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
                         id = NULL,
                         col = "black",
                         lwd = 1,
                         ...){
  if (inherits(x, "list") | inherits(x, "iefourier_lst")) {
    if(is.null(id)){
      invisible(lapply(x, plot_contour, id, col, lwd))
    } else{
      invisible(lapply(x[id], plot_contour, id, col, lwd))
    }
  } else{
    coord <- poly_check(x)
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
      return(lapply(x, plot_polygon,  fill, random_fill, points, merge, border, alpha, add, aspect_ratio = aspect_ratio, ...))
    }
  } else{
    if(inherits(x, "landmarks_regradi")){
      coords <- x$coords
    } else{
      coords <- poly_check(x)
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
                      id = NULL,
                      col = "black",
                      cex = 1,
                      lwd = 1){
  if (inherits(x, "list")) {
    if(is.null(id)){
      return(lapply(x, plot_mass, y, id, col, cex, lwd))
    } else{
      return(lapply(x[id], plot_mass, y, id, col, cex, lwd))
    }
  } else{
    if(is.data.frame(x) | is.matrix(x)){
      vx <- x[, 1]
      vy <- x[, 2]
    }
    points(vx, vy, pch = 16, col = col, cex = cex)
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



#' Width at a given height
#'
#' The function computes the polygonal convex hull of the points in x and then
#' returns the number of points that lie below a specified set of heights along
#' the vertical axis of the convex hull.
#'
#' @details
#'
#' The convex hull computed from x is aligned along the major axis and then
#' converted to a binary image. For each height in the at vector, the function
#' computes the corresponding row number in the binary image (i.e., the row
#' number that corresponds to the specified height along the vertical axis of
#' the convex hull) and sums the values in that row to obtain the number of
#' points that lie below the specified height. If the convex hull contains
#' multiple polygons and unify = FALSE, the function loops over each polygon
#' and returns a list of the number of points below the specified heights for
#' each polygon. If the convex hull contains only one polygon or multiple
#' polygons and unify = TRUE, the function returns a vector of the number of
#' points below the specified heights for that single polygon.
#'
#'
#' @param x A vector containing two-dimensional data points (often produced with
#'   [object_contour]).
#' @param at A vector of heights along the vertical axis of the convex hull at
#'   which to count the number of points below. The default value is `c(0.05,
#'   0.25, 0.5, 0.75, 0.95)`, which means the function will return the number of
#'   points below the 5th, 25th, 50th, 75th, and 95th percentiles of the convex
#'   hull. If `at = "heights"` is used, the function returns the width for each
#'   point of the object length.
#' @param unify A logical value indicating whether to use the unified convex
#' hull calculation method. If unify = TRUE, coordinates in x will be
#' first bound before computing the convex hull.
#' @param plot A logical value that specifies whether the widths should be
#'   plotted.
#' @return A vector with the widths of the convex hull at the specified heights
#'   or a list of vectors with the widths of each component.
#'
#' @export
#'
#' @examples
#' cont <- contours[[2]]
#' plot_polygon(cont |> conv_hull() |> poly_align())

#' #  width below 5th, 25th, 50th, 75th, and 95th percentiles of the length
#' wd <- poly_width_at(cont)
#' wd
#'
#' # width along the height
#' poly_width_at(cont, at = "height", plot = TRUE)
#'
#'
poly_width_at <- function(x,
                          at = c(0.05, 0.25, 0.5, 0.75, 0.95),
                          unify = FALSE,
                          plot = FALSE){
  if(!is.numeric(at) && any(at != "height")){
    warning("`at` must be one of 'height' or a numeric vector in the range 0-1.")
  }
  if(isTRUE(unify)){
    chu <- conv_hull_unified(x) |> poly_align(plot = FALSE)
  } else{
    chu <- conv_hull(x) |> poly_align(plot = FALSE)
  }
  if(inherits(chu, "list")){
    res <- list()
    for(i in seq_along(chu)){
      bin <- EBImage::Image(polygon_to_binary(chu[[i]]))
      if(is.character(at)){
        wdts <- sum_true_cols(bin)
        wdts <- wdts[2:(length(wdts) - 1)]
        res[[i]] <- wdts
      } else{
        heights <- ceiling(dim(bin)[2] * at)
        res[[i]] <-
          sapply(seq_along(heights), function(i){
            sum(bin[, heights[i]:heights[i]]@.Data == TRUE)
          })
      }
    }
    names(res) <- paste0(1:length(res))
    return(res)
  } else{
    bin <- EBImage::Image(polygon_to_binary(chu))
    if(is.character(at)){
      wdts <- sum_true_cols(bin)
      wd <- wdts[2:(length(wdts) - 1)]
    } else{
      heights <- ceiling(dim(bin)[2] * at)
      wd <-
        sapply(seq_along(heights), function(i){
          sum(bin[, heights[i]:heights[i]]@.Data == TRUE)
        })
    }
    if(isTRUE(plot)){
      plot(wd,
           1:length(wd),
           type = "l",
           axes = FALSE,
           ylab = "Height",
           xlab = "Width (pixels)")
      axis(1)
      axis(2)
    }
    return(wd)
  }
}


#' Get the pixel indices for a given row of a binary image
#'
#' This function finds the first row in the bin matrix that has a value greater
#' than 0 (TRUE). It then calculates the minimum, median, and maximum values for
#' the pixels in that row and creates an array containing the row index, the
#' minimum pixel index, the median pixel index, and the maximum pixel index.
#'
#' @param bin A logical matrix representing a binary image
#' @param row An optional row index. If not provided, the function selects the
#' first non-zero row.
#' @param direction The direction for row selection when row is not provided.
#' If set to `"updown"`, the function starts scanning from the top of the image
#' towards the bottom. If set to `"downup"`, the function starts scanning from
#' the bottom towards the top.
#'
#' @return A numeric vector containing the row index, the minimum pixel index,
#' the median pixel index, and the maximum pixel index.
#' @importFrom stats median
#' @examples
#' library(pliman)
#' leaf <- image_pliman("sev_leaf.jpg")
#' bin <- image_binary(leaf, "NB")[[1]]
#'
#' # first row with leaf (17)
#' pixel_index(bin)
#'
#' # index at the row 100
#' pixel_index(bin, row = 100)
#'
#' plot(leaf)
#' points(x = 248, y = 17, pch = 16, col = "red", cex = 2)
#' points(x = 163, y = 100, pch = 16, col = "red", cex = 2)
#' points(x = 333, y = 100, pch = 16, col = "red", cex = 2)
#'
#'
#' @export
pixel_index <- function(bin,
                        row = NULL,
                        direction = "updown"){
  if(is.null(row)){
    nonzcols <- which(colSums(bin) > 0)
    if(direction == "updown"){
      i <- nonzcols[1]
    } else{
      i <- nonzcols[length(nonzcols)]
    }
    pixels <- which(bin[, i] == TRUE)
    indexes <- c(i, min(pixels), ceiling(median(pixels)), max(pixels))
  }else{
    if(length(row) > 1){
      warning("'row' must be an escalar. The fist element will be used")
      row <- row[1]
    }
    pixels <- which(bin[, row] == TRUE)
    indexes <- c(row, min(pixels), ceiling(median(pixels)), max(pixels))
  }
  return(indexes)
}


#' Calculate the apex and base angles of an object
#'
#' This function calculates the apex and base angles of an object. It takes as
#' input a matrix of coordinates and returns the apex angle, base angle, and the
#' coordinates of the apex and base as a list. The angles are computed after the
#' object is aligned in the vertical axis with `poly_align()`.
#'
#' @param x A matrix of coordinates representing the contour of the object,
#'   often obtained with [object_contour()].
#' @param percentiles A numeric vector of two percentiles between 0 and 1
#' indicating the height of the points from the top to the bottom. The function
#' calculates the apex angle between the two percentiles and the base angle
#' between the lowest point and the highest point.
#' @param invert If `TRUE`, aligns the object along the horizontal axis.
#' @param plot Plots the polygon with the points? Defaults to `TRUE`.
#' @return A list containing the apex angle, base angle, apex coordinates, and
#' base coordinates.
#' @export
#'
#' @examples
#' library(pliman)
#' # a matrix of coordinates
#' angls <- poly_apex_base_angle(contours[[2]])
#' angls
#'
#' # or a list of coordinates
#' poly_apex_base_angle(contours)
poly_apex_base_angle <- function(x,
                                 percentiles = c(0.25, 0.75),
                                 invert = FALSE,
                                 plot = TRUE){
  if(inherits(x, "list")){
    res <- lapply(x, poly_apex_base_angle, percentiles, invert, plot = FALSE)
    angles <- data.frame(do.call(rbind, lapply(res, function(x){c(x$apex_angle, x$base_angle)})))
    angles$id <- rownames(angles)
    colnames(angles) <- c("apex_angle", "base_angle", "id")
    return(angles[, c(3, 1, 2)])
  } else{
    x <- poly_align(x, plot = FALSE)
    if(isTRUE(invert)){
      x <- poly_rotate(x, angle = 90)
    }
    bp <- x[which.min(x[, 2]), ]
    tp <- x[which.max(x[, 2]), ]

    # base points
    yqant <- ((tp[2] - bp[2]) * percentiles[1])  + bp[2]
    yc1 <- which.min(abs(x[, 2] - yqant))
    yc2 <- which.min(abs(x[-seq(yc1 - 2, yc1 + 2), 2] - yqant))
    coord_bottom <- rbind(bp, x[yc1, ], x[yc2, ])
    rownames(coord_bottom) <- c("base", "left", "right")

    # top points
    yqant <- ((tp[2] - bp[2]) * percentiles[2])  + bp[2]
    yc1 <- which.min(abs(x[,2] - yqant))
    yc2 <- which.min(abs(x[-seq(yc1 - 2, yc1 + 2), 2] - yqant))
    coord_top <- rbind(tp, x[yc1, ], x[yc2, ])
    rownames(coord_top) <- c("apex", "left", "right")

    if(isTRUE(plot)){
      plot_polygon(x)
      points(coord_top, col = "red", pch = 16)
      points(coord_bottom, col = "blue", pch = 16)
    }
    return(list(apex_angle = help_poly_angles(coord_top)[1],
                base_angle =  help_poly_angles(coord_bottom)[1],
                apex_coords = coord_top,
                base_coords = coord_bottom))
  }
}



#' Compute Perimeter Complexity Value (PCV)
#'
#' This function calculates the Perimeter Complexity Value (PCV) for a given set
#' of coordinates representing a contour. The PCV measures the variation of
#' distances between the original coordinates and the smoothed coordinates
#' relative to the perimeter length of the original contour. See more in
#' `details` section.
#'
#' @param x A matrix or a list of matrices representing the coordinates of the
#'   contour(s).
#' @param niter An integer specifying the number of iterations for smoothing the
#'   contour. See [poly_smooth()] for more details. Defaults to 100.
#' @return The PCV value(s) computed for the contour(s).
#'
#' @details
#'
#' The PCV is computed using the following formula: \deqn{PCV =
#' \frac{sum(dists) \times sd(dists)}{perim}} where \eqn{dists}
#' represents the distances between corresponding points in the original and
#' smoothed coordinates, and \eqn{perim} is the perimeter length of the smoothed
#' contour.
#'
#' The PCV is computed by first smoothing the input contour using a specified
#' number of iterations. The smoothed contour is then used to compute the
#' distances between corresponding points in the original and smoothed
#' coordinates. These distances reflect the variations in the contour shape
#' after smoothing. The sum of these distances represents the overall magnitude
#' of the variations. Next, the sum of distances is multiplied by the standard
#' deviation of the distances to capture the dispersion or spread of the
#' variations. Finally, this value is divided by the perimeter length of the
#' original contour to provide a relative measure of complexity. Therefore, the
#' PCV provides a relative measure of complexity by considering both the
#' magnitude and spread of the variations in the contour shape after smoothing.
#'
#' @param x A matrix or a list of matrices representing the coordinates of the
#'   polygon(s).
#' @param niter An integer specifying the number of smoothing iterations. See
#'   [poly_smooth()] for more details.
#'
#' @return If `x` is a matrix, returns the complexity value of the polygon's
#'   perimeter. If `x` is a list of matrices, returns a numeric vector of
#'   complexity values for each polygon.
#'
#' @examples
#' library(pliman)
#' set.seed(20)
#' shp <- efourier_shape(npoints = 1000)
#' poly_pcv(shp)
#'
#' # increase the complexity of the outline
#' shp2 <- poly_jitter(shp, noise_x = 20, noise_y = 250, plot = TRUE)
#'
#' smo <- poly_smooth(shp2, niter = 100, plot = FALSE)
#' plot_contour(smo, col = "red")
#' poly_pcv(shp2)
#'
#' @export
#'
poly_pcv <- function(x, niter = 100){
  complexity <- function(x, niter){
    x_smoth <- poly_smooth(x, niter = niter, plot = FALSE)
    dists <- sqrt((x[, 1] - x_smoth[, 1]) ^ 2 + (x[, 2] - x_smoth[, 2]) ^ 2)
    perim <- poly_perimeter(x)
    pcv <- (sum(dists) * sd(dists) ) / perim
    return(pcv)
  }
  if (inherits(x, "list")) {
    sapply(x, complexity, niter)
  } else{
    return(complexity(x, niter))
  }
}
