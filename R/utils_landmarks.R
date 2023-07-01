help_angle <- function(v1, v2){
  v1 <- complex(1, v1[1], v1[2])
  v2 <- complex(1, v2[1], v2[2])
  (pi + Arg(v1) - Arg(v2)) %% (2*pi) - pi
}

create_mat <- function(vect, nr){
  m1 <- matrix(NA, nr, nr)
  m1[lower.tri(m1, diag=FALSE)] <- vect
  mat <- t(m1)
  mat[lower.tri(mat, diag=FALSE)] <- vect
  diag(mat) <- 0
  rownames(mat) <- paste0("L", 1:ncol(mat))
  colnames(mat) <- paste0("L", 1:ncol(mat))
  return(mat)
}


#' Create image landmarks
#'
#' An interactive section where the user will be able to click on the image to
#' select landmarks manually is open. With each mouse click, a point is drawn
#' and an upward counter is shown in the console. After `n` counts or after the
#' user press Esc, the interactive process is interrupted and a `data.frame`
#' with the `x` and `y` coordinates for the landmarks is returned.
#'
#' @param img An `Image` object.
#' @param n The number of landmarks to produce. Defaults to `Inf`. In this case,
#'   landmarks are chosen up to the user press Esc.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param scale A known scale of the coordinate values. If `NULL` (default)
#'   `scale = 1` is used.
#' @param calibrate A logical argument indicating whether a calibration step
#'   must be performed before picking up the landmarks. If so, [calibrate()] is
#'   called internally. Users must then select two points and indicate a known
#'   distance. A scale value will internally be computed and used in the
#'   correction of the coordinates (from pixels to the unit of the known
#'   distance).
#' @return A `data.frame` with the `x` and `y`-coordinates from the landmarks.
#' @importFrom utils combn
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#'
#' @examples
#' if(isTRUE(interactive())){
#' library(pliman)
#' img <- image_pliman("potato_leaves.jpg")
#' x <- landmarks(img)
#' }
landmarks <- function(img,
                      n = Inf,
                      viewer = get_pliman_viewer(),
                      scale = NULL,
                      calibrate = FALSE){
  if (isTRUE(interactive())) {
    if(isTRUE(calibrate)){
      scale <- calibrate(img)
      message(paste0("scale: ", scale))
    } else{
      if(is.null(scale)){
        scale <- 1
      } else{
        scale <- scale
      }
    }
    viewopt <- c("base", "mapview")
    viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]
    if(viewopt == "base"){
      plot(img)
      on.exit(return(coords))
      message("Use the first mouse button to select landmarks in the plot.\nPress Esc to exit.")
      i <- 1
      coords <- data.frame(x = NA, y = NA)
      while (i <= n) {
        d <- unlist(locator(n = 1))
        if (is.null(d)) {
          break
        }
        coords[i, 1] <- d[[1]] / scale
        coords[i, 2] <- d[[2]] / scale
        rownames(coords)[[i]] <- paste0("L", i)
        points(d[[1]], d[[2]], type = "p", col = "red", cex = 1, pch = 19)
        text(d[[1]], d[[2]], pos=2, labels=paste0("L", i))
        cat("Number of landmarks:", i, "\r")
        i <- i + 1
      }
      if (i >= n) {
        warning("Number of landmarks achieved.", call. = FALSE)
      }
    } else{
      coords <- mv_points(img, title = "Use the first mouse button to select landmarks in the plot. Press 'Done' to exit.")
      rownames(coords) <- paste0("L", 1:nrow(coords))
      return(coords)
    }
  }
}

#' Pseudolandmarks with equally spaced angles
#'
#' Select `n` landmarks that are spaced with a regular sequence of angles taken
#' between the outline coordinates and the centroid.
#'
#' @param x A `matrix`, a `data.frame` a `list` of perimeter coordinates, often
#'   produced with [object_contour()].
#' @param n Number of points to be sampled. Defaults to 50.
#' @param close Return a closed polygon? Defaults to `TRUE`.
#' @param plot Create a plot? Defaults to `TRUE`.
#' @param ncol,nrow The number of rows or columns in the plot grid when a `list`
#'   is used in `x`. Defaults to `NULL`, i.e., a square grid is produced.
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @return A list with the following objects:
#' * `pixindices`: Vector of radius indices.
#' * `radii`: Vector of sampled radii lengths.
#' * `Xc`: The centroid coordinate of `x` axis.
#' * `Yc`: The centroid coordinate of `y` axis.
#' * `coords`: Coordinates of sampled points arranged in a two-column matrix.
#'
#' If `x` is a list, a list of objects described above is returned.
#' @export
#' @note Borrowed from Claude (2008), pp. 53
#'
#' @examples
#' library(pliman)
#' plot_polygon(contours[[1]])
#' ldm <- landmarks_regradi(contours)
#'
landmarks_regradi <- function(x,
                              n = 50,
                              close = TRUE,
                              plot = TRUE,
                              ncol = NULL,
                              nrow = NULL){
  if (inherits(x, "list")) {
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
    coords <- lapply(x, landmarks_regradi, n, close, plot)
    return(coords)
  } else{
    Rx <- x[,1]
    Ry <- x[,2]
    le <- length(Rx)
    M <- matrix(c(Rx, Ry), le, 2)
    M1 <- matrix(c(Rx - mean(Rx), Ry - mean(Ry)), le, 2)
    V1 <- complex(real=M1[,1], imaginary=M1[,2])
    M2 <- matrix(c(Arg(V1), Mod(V1)), le, 2)
    V2 <- NA
    for (i in 0:(n-1)){
      V2[i + 1] <- which.max((cos(M2[,1]-2 * i * pi / n)))
    }
    V2 <- sort(V2)
    coords <- M1[V2,] |> as.data.frame()
    colnames(coords) <- c("x", "y")
    coords$x <- coords[,1] + mean(Rx)
    coords$y <- coords[,2] + mean(Ry)
    if(isTRUE(close)){
      coords <- poly_close(coords)
    }
    Xc <- mean(coords[,1])
    Yc <- mean(coords[,2])
    if(isTRUE(plot)){
      plot(coords[,1],
           coords[,2],
           type = "l",
           lwd = 1.5,
           asp = 1,
           xlab = NA,
           ylab = NA,
           axes = FALSE)
      axis(1)
      axis(2)
      points(coords[,1], coords[,2], pch = 16)
      segments(mean(Rx),
               mean(Ry),
               coords[,1],
               coords[,2],
               col = "gray")
    }
    return(
      structure(
        list(pixindices = V2,
             radii = M2[V2,2],
             Xc = Xc,
             Yc = Yc,
             coords = coords),
        class = "landmarks_regradi"
      )
    )
  }
}



#' Artificially inflates the number of landmarks
#'
#' Interpolates supplementary landmarks that correspond to the mean coordinates
#' of two adjacent landmarks.
#'
#' @param x A `matrix`, a `data.frame` a `list` of perimeter coordinates, often
#'   produced with [object_contour()], [landmarks()], or [landmarks_regradi()].
#' @param n The number of iterations. Defaults to 3.
#' @param smooth_iter The number of smoothing iterations to perform. This will
#'   smooth the perimeter of the interpolated landmarks using [poly_smooth()].
#' @param plot Creates a plot? Defaults to `TRUE`.
#' @param ncol,nrow The number of rows or columns in the plot grid when a `list`
#'   is used in `x`. Defaults to `NULL`, i.e., a square grid is produced.
#' @return A Matrix of interpolated coordinates.
#' @export
#'
#' @examples
#' library(pliman)
#'
#' # equally spaced landmarks
#' plot_polygon(contours[[4]])
#' ldm <- landmarks_regradi(contours[[4]], plot = FALSE)
#' points(ldm$coords, pch = 16)
#' segments(mean(ldm$coords[,1]),
#'          mean(ldm$coords[,2]),
#'          ldm$coords[,1],
#'          ldm$coords[,2])
#'
#' ldm_add <- landmarks_add(ldm, plot = FALSE)
#' points(ldm_add, col = "red")
#' points(ldm$coords, pch = 16)
#'
#' # smoothed version
#' ldm_add_smo <- landmarks_add(ldm, plot = FALSE, smooth_iter = 10)
#' lines(ldm_add_smo, col = "blue", lwd = 3)
landmarks_add <- function(x,
                          n = 3,
                          smooth_iter = 0,
                          plot = TRUE,
                          nrow = NULL,
                          ncol = NULL){
  if (inherits(x, "list")) {
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
    coords <- lapply(x, landmarks_add, n, smooth_iter, plot)
    return(coords)
  } else{
    if(inherits(x, "landmarks_regradi")){
      M <- as.matrix(x$coords)
    } else{
      M <- as.matrix(x)
    }
    a <- 0
    while(a <= n) {
      p <- dim(M)[1]
      k <- dim(M)[2]
      N <- matrix(NA,2*p,k)
      N[((1:p)*2)-1,] <- M
      N[(1:p)*2,]<-(M+(rbind(M[-1,],M[1,])))/2
      M <- N
      a <- a + 1
    }
    if(smooth_iter != 0){
      M <- poly_smooth(M, niter = smooth_iter, plot = FALSE) |> as.matrix()
    }
    if(isTRUE(plot)){
      plot_polygon(M)
    }
    return(M)
  }
}


#' Calibrates distances of landmarks
#'
#' Calibrating the actual size is possible if any interlandmark distance on the
#' image is known. calibrate() can be used to determine the size of a known
#' distance (cm) on the graph. I invite users to photograph the object together
#' with a scale (e.g., ruler, micrometer...).
#'
#' @param img An `Image` object
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @return A numeric (`double`) scalar value indicating the scale (in pixels per
#'   unit of known distance).
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#'
#' @examples
#' if(isTRUE(interactive())){
#' library(pliman)
#' #### compute scale (dots per unit of known distance) ####
#' # only works in an interactive section
#' # objects_300dpi.jpg has a known resolution of 300 dpi
#' img <- image_pliman("objects_300dpi.jpg")
#' # Larger square: 10 x 10 cm
#' # 1) Run the function calibrate()
#' # 2) Use the left mouse button to create a line in the larger square
#' # 3) Declare a known distance (10 cm)
#' # 4) See the computed scale (pixels per cm)
#' calibrate(img)
#'
#' # scale ~118
#' # 118 * 2.54 ~300 DPI
#' }
calibrate <- function(img, viewer = get_pliman_viewer()){
  viewopt <- c("base", "mapview")
  viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]
  if(viewopt == "base"){
    message("Use the first mouse button to create a line in the plot.")
    plot(img)
    a <- locator(2,type="o",pch=8,lwd=2,col="red",lty="11")
    scale <- sqrt(sum(diff(a$x)^2 + diff(a$y)^2))
  } else{
    mv <- mv_two_points(img,
                        title = "Use the 'Draw Polyline' tool to create a line in the plot")
    a <- abs(mv[[1]] - mv[[3]])
    b <- abs(mv[[2]] - mv[[4]])
    scale <- sqrt(a ^ 2 + b ^ 2)
  }
  known <- as.numeric(readline("known distance (cm): "))
  scale / known
}


#' Distances between landmarks
#'
#' Computes the distance between two landmarks as the square root of the sum of
#' the squared differences between each coordinate (Claude, 2008).
#'
#' @param x An object computed with [landmarks()].
#'
#' @return A matrix with the distances for each landmark combination.
#' @note Borrowed from Claude (2008), pp. 49
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#'
#' @examples
#' if(isTRUE(interactive())){
#' library(pliman)
#' img <- image_pliman("potato_leaves.jpg")
#' x <- landmarks(img)
#' landmarks_dist(x)
#' }
#'

landmarks_dist <- function(x){
  if(inherits(x, "landmarks_regradi")){
    x <- x$coords
  }
  n <- nrow(x)
  create_mat(dist(x), n)
}


#' Angles between landmarks
#'
#' Computes the angle from two interlandmark vectors using the difference of
#' their arguments using complex vectors (Claude, 2008).
#'
#' @param x An object computed with [landmarks()].
#' @param unit The unit of the angle. Defaults to radian (rad). Use `unit =
#'   "deg"` to return the angles in degrees.
#'
#' @return A matrix with the angles for each landmark combination.
#' @note Borrowed from Claude (2008), pp. 50
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#'
#' @examples
#' if(isTRUE(interactive())){
#' library(pliman)
#' img <- image_pliman("potato_leaves.jpg")
#' x <- landmarks(img)
#' landmarks_angle(x)
#' }
landmarks_angle <- function(x, unit = c("rad", "deg")){
  unit <- unit[[1]]
  if(!unit %in% c("rad", "deg")){
    warning("`unit` must be 'rad' or 'deg'. Defaulting to 'rad'")
    unit <- "rad"
  }
  if(inherits(x, "landmarks_regradi")){
    x <- x$coords
  }
  angles <- combn(nrow(x), 2, FUN = function(y){
    help_angle(x[y[1], 1:2], x[y[2], 1:2])
  })
  if(unit == "deg"){
    angles <- angles * 180 / pi
  }
  n <- nrow(x)
  create_mat(angles, n)
}



#' Confidence ellipse
#'
#' Produces a confidence ellipse that is an iso-contour of the Gaussian
#' distribution, allowing to visualize a 2D confidence interval.
#'
#' @param x A `matrix`, a `data.frame` or a `list` of perimeter coordinates, often
#'   produced with [object_contour()].
#' @param conf The confidence level. Defaults to `0.95`
#' @param np Number of sampled points on the ellipse.
#' @param plot Create a plot? Defaults to `TRUE`.
#' @param fill The color to fill the ellipse. Defaults to `"green"`.
#' @param alpha The alpha value to define the opacity of ellipse. Defaults to
#'   `0.3`
#' @param random_fill Fill multiple ellipses with random colors? Defaults to
#'   `TRUE`.
#' @return A matrix with coordinates of points sampled on the ellipse.
#' @note Borrowed from Claude (2008), pp. 85
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @importFrom graphics segments
#' @importFrom stats cor qnorm
#' @export
#'
#' @examples
#' library(pliman)
#' ellipse(contours)
ellipse <- function(x,
                    conf = 0.95,
                    np = 100,
                    plot = TRUE,
                    fill = "green",
                    alpha = 0.3,
                    random_fill = TRUE){
  if (inherits(x, "list")) {
    d <- lapply(x, ellipse, conf, np, plot =  FALSE)
    if(isTRUE(plot)){
      ifelse(random_fill == TRUE,
             colf <- random_color(n = length(x)),
             colf <- rep(fill, length(x)))
      plot_polygon(x, merge = TRUE, random_fill = TRUE)
      for (i in seq_along(x)) {
        plot_polygon(d[[i]],
                     add = TRUE,
                     fill = colf[[i]],
                     border = NA,
                     alpha = alpha,
                     aspect_ratio = 1)
      }
    }
    invisible(d)
  } else{
    y <- x[, 2]
    x <- x[, 1]
    centroid <- apply(cbind(x, y), 2, mean)
    ang <- seq(0, 2 * pi,length = np)
    z <- cbind(cos(ang) ,sin(ang))
    rcoef <- qnorm((1 - conf)/2, lower.tail=F)
    vcvxy <- var(cbind(x,y))
    r <- cor(x, y)
    M1 <- matrix(c(1,1,-1,1),2,2)
    M2 <- matrix(c(var(x), var(y)),2,2)
    M3 <- matrix(c(1+r,1-r),2,2, byrow=T)
    ellpar <- M1 * sqrt(M2 * M3/2)
    coords <- t(centroid + rcoef * ellpar %*%t(z))
    xlim <- c(min(c(min(x), min(coords[,1]))), max(c(max(x), max(coords[,1]))))
    ylim <- c(min(c(min(y), min(coords[,2]))), max(c(max(y), max(coords[,2]))))
    if(isTRUE(plot)){
      plot_polygon(coords,
                   add = FALSE,
                   fill = fill,
                   alpha = alpha,
                   xlim = xlim,
                   ylim = ylim,
                   merge = TRUE)
      points(x, y,
             pch = 16)

    }
    invisible(coords)
  }
}

