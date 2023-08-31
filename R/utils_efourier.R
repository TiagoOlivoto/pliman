#' Elliptical Fourier Analysis
#'
#' Computes Elliptical Fourier Analysis of closed outlines based on `x` and
#' `y`-coordinates coordinates.
#'
#' @param x A `matrix`, a `data.frame` a `list` of perimeter coordinates, often
#'   produced with [object_contour()] or a vector of landmarks produced with
#'   [landmarks()] or [landmarks_regradi()].
#' @param nharm An integer indicating the number of harmonics to use. Defaults
#'   to 10.
#' @param align Align the objects before computing Fourier analysis? Defaults to
#'   `FALSE`. If `TRUE`, the object is first aligned along the major caliper
#'   with [poly_align()].
#' @param center Center the objects on the origin before computing Fourier
#'   analysis? Defaults to `FALSE`. If `TRUE`, the object is first centered on
#'   the origin with [poly_center()].
#' @param smooth_iter The number of smoothing iterations to perform. This will
#'   smooth the perimeter of the objects using [poly_smooth()].
#' @return A list of class `efourier` with:
#' * the harmonic coefficients (`an`, `bn`, `cn` and `dn`)
#' * the estimates of the coordinates of the centroid of the configuration
#'   (`a0` and `c0`).
#' * The number of rows (points) of the perimeter outline (`nr`).
#' * The number of harmonics used (`nharm`).
#' * The original coordinates (`coords`).
#'
#' If `x` is a list of perimeter coordinates, a list of `efourier` objects will
#' be returned as an object of class `iefourier_lst`.
#' @details Adapted from Claude (2008). pp. 222-223.
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#'
#' Kuhl, F. P., and Giardina, C. R. (1982). Elliptic Fourier features of a
#' closed contour. Computer Graphics and Image Processing 18, 236–258. doi:
#' \doi{10.1016/0146-664X(82)90034-X}
#'
#' @export
#' @examples
#' library(pliman)
#' leaf1 <- contours[[4]]
#' plot_polygon(leaf1)
#'
#' #### default options
#' # 10 harmonics (default)
#' # without alignment
#'
#' ef <- efourier(leaf1)
#' efourier_coefs(ef)
#'
#' # object is aligned along the major caliper with `poly_align()`
#' # object is centered on the origin with `poly_center()`
#' # using a list of object coordinates
#' ef2 <- efourier(contours, align = TRUE, center = TRUE)
#' efourier_coefs(ef2)
#'
#' # reconstruct the perimeter of the object
#' # Use only the first one for simplicity
#' plot_polygon(contours[[1]] |> poly_align() |> poly_center())
#' efourier_inv(ef2[[1]]) |> plot_contour(col = "red", lwd = 4)
efourier <- function(x,
                     nharm = 10,
                     align = FALSE,
                     center = FALSE,
                     smooth_iter = 0) {
  if (inherits(x, "list")) {
    if(inherits(x[[1]][[1]], "landmarks_regradi")){
      x <- lapply(x, function(x){
        x[[1]]$coords
      })
    }
    if(inherits(x[[1]], "landmarks_regradi")){
      x <- lapply(x, function(x){
        x$coords
      })
    }
    coords <- lapply(x, efourier, nharm, align, center, smooth_iter)
    return(structure(coords, class = "efourier_lst"))
  } else{
    if(inherits(x, "landmarks_regradi")){
      coord <-x$coords
    } else{
      coord <- poly_check(x) |> poly_unclose()
    }
    if(isTRUE(align)){
      coord <- poly_align(coord, plot = FALSE)
    }
    if(isTRUE(center)){
      coord <- poly_center(coord, plot = FALSE)
    }
    nr <- nrow(coord)
    if (nharm * 2 > nr) {
      nharm <- floor(nr/2)
      message("'nharm' must be lower than half the number of points. It has been set to ", nharm, " harmonics")
    }
    if (nharm == -1) {
      nharm = floor(nr/2)
      message("the number of harmonics used has been set to: ", nharm)
    }
    if (smooth_iter != 0) {
      coord <- poly_smooth(coord, niter = smooth_iter)
    }
    Dx <- coord[, 1] - coord[, 1][c(nr, (1:(nr - 1)))]
    Dy <- coord[, 2] - coord[, 2][c(nr, (1:(nr - 1)))]
    Dt <- sqrt(Dx^2 + Dy^2)
    Dt[Dt < 1e-10] <- 1e-10  # to avoid Nan
    t1 <- cumsum(Dt)
    t1m1 <- c(0, t1[-nr])
    t <- sum(Dt)
    an <- bn <- cn <- dn <- numeric(nharm)
    for (i in 1:nharm) {
      ti <- (t/(2 * pi^2 * i^2))
      r <- 2 * i * pi
      an[i] <- ti * sum((Dx/Dt) * (cos(r * t1/t) - cos(r * t1m1/t)))
      bn[i] <- ti * sum((Dx/Dt) * (sin(r * t1/t) - sin(r * t1m1/t)))
      cn[i] <- ti * sum((Dy/Dt) * (cos(r * t1/t) - cos(r * t1m1/t)))
      dn[i] <- ti * sum((Dy/Dt) * (sin(r * t1/t) - sin(r * t1m1/t)))
    }
    a0 <- 2 * sum(coord[, 1] * Dt/t)
    c0 <- 2 * sum(coord[, 2] * Dt/t)
    coefs <- list(an = an,
                  bn = bn,
                  cn = cn,
                  dn = dn,
                  a0 = a0,
                  c0 = c0,
                  nr = nr,
                  nharm = nharm,
                  coords = coord)
    return(structure(coefs, class = "efourier"))
  }
}


#' Inverse Elliptical Fourier Analysis
#'
#' Performs an inverse elliptical Fourier transformation to construct a shape,
#' given a list with Fourier coefficients computed with [efourier()].
#' @param x An object of class `efourier` or `efourier_lst` computed with
#'   [efourier()].
#' @param nharm An integer indicating the number of harmonics to use. If not
#'   specified the number of harmonics used in `x` is used.
#' @param a0,c0 the estimates of the coordinates of the centroid of the
#'   configuration. If `NULL` (default), the generated coordinates will be
#'   centered on the position of the original shape given by [efourier()].
#' @param npoints The number of interpolated points on the constructed outline.
#'   Defaults to 500.
#' @details Adapted from Claude (2008). pp. 223.
#' @references Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#' @examples
#' library(pliman)
#' plot_polygon(contours, aspect_ratio = 1)
#' # without alignment
#' ef <- efourier(contours, nharm = 10, align = FALSE)
#' ief <- efourier_inv(ef)
#' plot_contour(ief, col = "red", lwd = 2)
efourier_inv <- function(x,
                         nharm = NULL,
                         a0 = NULL,
                         c0 = NULL,
                         npoints = 500) {
  if (inherits(x, "efourier_lst")) {
    coords <- lapply(x, efourier_inv, nharm, a0, c0, npoints)
    return(structure(coords, class = "iefourier_lst"))
  } else{
    if (is.null(x$a0)){
      x$a0 <- 0
    }
    if (is.null(x$c0)){
      x$c0 <- 0
    }
    an <- x$an
    bn <- x$bn
    cn <- x$cn
    dn <- x$dn
    a0 <- ifelse(is.null(a0), x$a0, a0)
    c0 <- ifelse(is.null(c0), x$c0, c0)

    if (is.null(nharm)) {
      nharm <- length(an)
    }
    if(nharm > length(an)){
      warning("Number of harmonics used must be less than or equal to the number used in `efourier()` function\nSetting `nharm` to ", length(an), call. = FALSE)
      nharm <- length(an)
    }
    theta <- seq(0, 2 * pi, length = npoints + 1)[-(npoints + 1)]
    hx <- matrix(NA, nharm, npoints)
    hy <- matrix(NA, nharm, npoints)
    for (i in 1:nharm) {
      hx[i, ] <- an[i] * cos(i * theta) + bn[i] * sin(i * theta)
      hy[i, ] <- cn[i] * cos(i * theta) + dn[i] * sin(i * theta)
    }
    x <- (a0/2) + apply(hx, 2, sum)
    y <- (c0/2) + apply(hy, 2, sum)
    coord <- cbind(x, y)
    colnames(coord) <- c("x", "y")
    return(structure(coord, class = "iefourier"))
  }
}

#' Erros between the original and reconstructed outline
#'
#' Computes the sum of squared distances between the original data and
#' reconstructed outline. It allows examining reconstructed outlines with the
#' addition of successive contributing harmonics indicated in the argument
#' `nharm`.
#'
#' @param x An object computed with [efourier()].
#' @param nharm An integer or vector of integers indicating the number of
#'   harmonics to use. If not specified the number of harmonics used in `x` is
#'   used.
#' @param type The type of plot to produce. By default, a line plot with the sum
#'   of squared distances (y-axis) and the number of harmonics (x-axis) is
#'   produced. If `type = "outline"` is used, a plot with the original polygon
#'   and the constructed outline is produced. If `type = "deviations"` is used,
#'   a plot with the deviations from the original outline and reconstructed
#'   outline (y-axis) and points along the outline (x-axis) is produced.
#' @param plot A logical to inform if a plot should be produced. Defaults to
#'   `TRUE`.
#' @param ncol,nrow The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @return A list with the objects:
#' * `dev_points` A list with the deviations (distances) from original and
#' predicted outline for each pixel of the outline.
#' *  `data.frame` object with the minimum, maximum and average
#' deviations (based on the outline points).
#'
#' If `x` is an object of class `efourier_lst`, a list will be returned.
#' @importFrom graphics legend
#' @export
#' @examples
#' library(pliman)
#' ef <-
#'   contours[[1]] |>
#'   efourier(nharm = 30)
#'
#' efourier_error(ef)
#'
#' efourier_error(ef,
#'                nharm = 30,
#'                type = "outline")
#'
#' efourier_error(ef,
#'                nharm = c(1, 4, 20),
#'                type = "deviations")
efourier_error <- function(x,
                           nharm = NULL,
                           type = c("error", "outline", "deviations"),
                           plot = TRUE,
                           ncol = NULL,
                           nrow = NULL){
  type <- type[[1]]
  if(!type %in% c("error", "outline", "deviations")){
    warning("`type` argument must be one of 'error', 'outline', or 'deviations'. Defaulting to 'error'", call. = FALSE)
    type <- "error"
  }
  if(inherits(x, "efourier_lst")){
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
    error <- lapply(x, efourier_error, nharm, type, plot)
    err <-
      do.call(rbind,
              lapply(error, function(x){x$stats}))

    err$object <- sapply(strsplit(rownames(err), "\\."), "[", 1)
    err <- err[, c(4, 1, 2, 3)]
    rownames(err) <- NULL

    dev_points <- lapply(error, function(x){x$dev_points})
    return(list(
      stats = err,
      dev_points = dev_points
    ))
  } else{
    ef <- x[["coords"]]
    if(is.null(nharm)){
      nh <- 1:x$nharm
    } else{
      nh <- nharm
    }
    averagedev <- maxdev <- mindev <- numeric(length(nh))
    if(isTRUE(plot) & type == "outline"){
      num_plots <- length(nh)
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
    }
    dev_points <- list()
    for (i in seq_along(nh)) {
      index <- nh[[i]]
      ief <- efourier_inv(x, nharm = index, npoints = nrow(ef)) |> unclass()
      deviation <- sqrt(apply((ef-ief)^2, 1, sum))
      dev_points[[paste0(nh[[i]])]] <- deviation
      averagedev[i] <- mean(deviation)
      maxdev[i] <- max(deviation)
      mindev[i] <- min(deviation)
      if(isTRUE(plot)){
        if(type == "outline"){
          plot(ef,
               axes = FALSE,
               type = "n",
               xlab = "",
               ylab = "",
               asp = 1,
               main=paste("Harmonics 0 to",index))
          axis(1)
          axis(2)
          polygon(ef,
                  col = rgb(t(col2rgb("gray") / 255),  alpha = 1),
                  border = NA)
          plot_contour(ief, col = "red", lwd = 2)
        }
      }
    }
    if(isTRUE(plot)){
      if(type == "error"){
        plot(nh, averagedev, pch=16,
             ylab = "Average deviation from original outline",
             xlab = "Number of harmonics")
        lines(nh, averagedev)
      }
      if(type == "deviations"){
        max <- sapply(dev_points, max)
        xl <- length(dev_points[[1]])
        yl <- length(dev_points[[1]])
        plot(x = NA,
             xlim = c(1, xl),
             ylim = c(0, max(max)),
             xlab = "Points along the outline",
             ylab = "Deviation from the original outline (pixels)")
        cols <- ggplot_color(n = length(dev_points))
        legends <- character()
        for (i in 1:length(dev_points)) {
          lines(x = 1:xl,
                y = dev_points[[i]],
                col = cols[[i]],
                lwd = 2)
          legends <- append(legends, paste0("H", names(dev_points[i])))
        }
        legend(x = "topright",
               legend = legends,
               col = cols,
               lwd = 2,
               lty = 1,
               pch = c(NA,NA),
               cex = 0.5)
      }
    }
    invisible(list(dev_points = dev_points,
                   stats = data.frame(average_dev = averagedev,
                                      max_dev = maxdev,
                                      min_dev = mindev)))
  }
}



#' Normalized Fourier coefficients
#'
#' The first harmonic defines an ellipse that best fits the outlines. One can
#' use the parameters of the first harmonic to “normalize” the data so that they
#' can be invariant to size, rotation, and starting position of the outline
#' trace. This approach is referred to in the literature as the normalized
#' elliptic Fourier. [efourier_norm()] calculates a new set of Fourier
#' coefficients `An`, `Bn`, `Cn`, `Dn` that one can use for further multivariate
#' analyses (Claude, 2008).
#'
#' @param x An object computed with [efourier()].
#' @param start Logical value telling whether the position of the starting point
#'   has to be preserved or not.
#'
#' @return A list with the following components:
#' * `A`, `B`, `C`, `D` for harmonic coefficients.
#' * `size` the magnitude of the semi-major axis of the first fitting ellipse.
#' * `theta` angle, in radians, between the starting and the semi-major axis of
#' the first fitting ellipse.
#' * `psi` orientation of the first fitting ellipse
#' * `a0` and `c0`, harmonic coefficients.
#' * `lnef` the concatenation of coefficients.
#' * `nharm` the number of harmonics used.
#' @details Adapted from Claude (2008). pp. 226.
#' @references Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#' @examples
#' library(pliman)
#' leaf1 <- contours[[4]]
#' plot_polygon(leaf1)
#'
#' # compute the Fourier coefficients
#' ef <- efourier(leaf1)
#' efourier_coefs(ef)
#'
#' # Normalized Fourier coefficients
#'
#' efn <- efourier_norm(ef)
#' efourier_coefs(efn)

efourier_norm <- function(x, start = FALSE) {
  if (inherits(x, "efourier_lst")) {
    coords <- lapply(x, efourier_norm, start)
    return(structure(coords, class = "nefourier_lst"))
  } else{
    A1 <- x$an[1]
    B1 <- x$bn[1]
    C1 <- x$cn[1]
    D1 <- x$dn[1]
    nharm <- length(x$an)
    theta <- 0.5 * atan(2 * (A1 * B1 + C1 * D1)/(A1^2 + C1^2 - B1^2 - D1^2))%%pi
    phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta),
                           cos(theta)), 2, 2)
    M2 <- matrix(c(A1, C1, B1, D1), 2, 2) %*% phaseshift
    v <- apply(M2 ^ 2, 2, sum)
    if (v[1] < v[2]) {
      theta <- theta + pi/2
    }
    theta <- (theta + pi/2) %% pi - pi/2
    Aa <- A1 * cos(theta) + B1 * sin(theta)
    Cc <- C1 * cos(theta) + D1 * sin(theta)
    scale <- sqrt(Aa^2 + Cc^2)
    psi <- atan(Cc / Aa) %% pi
    if (Aa < 0) {
      psi <- psi + pi
    }
    size <- 1/scale
    rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)),
                       2, 2)
    A <- B <- C <- D <- numeric(nharm)
    if (start) {
      theta <- 0
    }
    for (i in 1:nharm) {
      mat <- size * rotation %*%
        matrix(c(x$an[i], x$cn[i],
                 x$bn[i], x$dn[i]), 2, 2) %*%
        matrix(c(cos(i * theta), sin(i * theta),
                 -sin(i * theta), cos(i * theta)), 2, 2)
      A[i] <- mat[1, 1]
      B[i] <- mat[1, 2]
      C[i] <- mat[2, 1]
      D[i] <- mat[2, 2]
      lnef <- c(A[i], B[i], C[i], D[i])
    }
    coefs <-
      list(A = A,
           B = B,
           C = C,
           D = D,
           size = scale,
           theta = theta,
           psi = psi,
           a0 = x$a0,
           c0 = x$c0,
           lnef = lnef,
           nharm = nharm)
    return(structure(coefs, class = "nefourier"))
  }
}

#' Get Fourier coefficients
#'
#' Extracts the Fourier coefficients from objects computed with [efourier()] and
#' [efourier_norm()] returning a 'ready-to-analyze' data frame.
#'
#' @param x An object computed with [efourier()] or [efourier_norm()].
#'
#' @return A `data.frame` object
#' @export
#'
#' @examples
#' library(pliman)
#'
#' # a list of objects
#' efourier(contours) |> efourier_coefs()
#'
#' # one object, normalized coefficients
#' efourier(contours[[4]]) |>
#'   efourier_norm() |>
#'   efourier_coefs()
efourier_coefs <- function(x){
  if(!inherits(x, c("nefourier_lst", "efourier_lst", "efourier", "nefourier"))){
    stop("Object is not valid. Please, use an object computed with `efourier()` or `efourier_norm()`", call. = FALSE)
  }
  if(inherits(x, "efourier_lst") | inherits(x, "nefourier_lst")){
    if(inherits(x[[1]], "nefourier_lst")){
      x <- structure(lapply(x, function(x){x[[1]]}), class = "nefourier_lst")
    }
    if(inherits(x[[1]], "efourier_lst")){
      x <- structure(lapply(x, function(x){x[[1]]}), class = "efourier_lst")
    }
    nharm <- x[[1]]$nharm
    names_objs <- names(x)
    if(inherits(x, "nefourier_lst")){
      an <- do.call(rbind, lapply(x, function(x){x[["A"]]}))
      bn <- do.call(rbind, lapply(x, function(x){x[["B"]]}))
      cn <- do.call(rbind, lapply(x, function(x){x[["C"]]}))
      dn <- do.call(rbind, lapply(x, function(x){x[["D"]]}))
      colnames(an) <- paste0("A", 1:nharm)
      colnames(bn) <- paste0("B", 1:nharm)
      colnames(cn) <- paste0("C", 1:nharm)
      colnames(dn) <- paste0("D", 1:nharm)
      coefs <- cbind(an, bn, cn, dn) |> as.data.frame()
      coefs$object <- names_objs
      coefs <- coefs[, c(ncol(coefs), 1:(ncol(coefs) - 1))]
    } else{
      if(inherits(x, "efourier_lst")){
        an <- do.call(rbind, lapply(x, function(x){x[["an"]]}))
        bn <- do.call(rbind, lapply(x, function(x){x[["bn"]]}))
        cn <- do.call(rbind, lapply(x, function(x){x[["cn"]]}))
        dn <- do.call(rbind, lapply(x, function(x){x[["dn"]]}))
        colnames(an) <- paste0("A", 1:nharm)
        colnames(bn) <- paste0("B", 1:nharm)
        colnames(cn) <- paste0("C", 1:nharm)
        colnames(dn) <- paste0("D", 1:nharm)
        coefs <- cbind(an, bn, cn, dn) |> as.data.frame()
        coefs$object <- names_objs
        coefs <- coefs[, c(ncol(coefs), 1:(ncol(coefs) - 1))]
      }
    }
  } else{
    if(inherits(x, "efourier")){
      an <- x[["an"]]
      bn <- x[["bn"]]
      cn <- x[["cn"]]
      dn <- x[["dn"]]
      coefs <- cbind(an, bn, cn, dn)
    } else{
      an <- x[["A"]]
      bn <- x[["B"]]
      cn <- x[["C"]]
      dn <- x[["D"]]
      coefs <- cbind(an, bn, cn, dn)
    }
  }
  return(data.frame(coefs))
}



#' Power in Fourier Analysis
#'
#' Computes an spectrum of harmonic Fourier power. The power is proportional to
#' the harmonic amplitude and can be considered as a measure of shape
#' information. As the rank of harmonic increases, the power decreases and adds
#' less and less information. We can evaluate the number of harmonics that we
#' must select, so their cumulative power gathers 99% of the total cumulative
#' power (Claude, 2008).
#'
#' @param x An object of class `efourier`computed with [efourier()].
#' @param first Logical argument indicating whether to include the first
#'   harmonic for computing the power. See `Details`.
#' @param thresh A numeric vector indicating the threshold power. The number of
#'   harmonics needed for such thresholds will then be computed.
#' @param plot Logical argument indicating whether to produce a plot.
#' @param ncol,nrow The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @details Most of the shape "information" is contained in the first harmonic.
#'   This is not surprising because this is the harmonic that best fits the
#'   outline, and the size of ellipses decreases as for explaining successive
#'   residual variation. However, one may think that the first ellipse does not
#'   contain relevant shape information, especially when differences one wants
#'   to investigate concern complex outlines. By using `first = FALSE` it is
#'   possible to remove the first harmonic for this computation. When working on
#'   a set of outlines, high-rank-harmonics can contain information that may
#'   allow groups to be distinguished (Claude, 2008).
#' @return A list with the objects:
#' * `cum_power`, a `data.frame` object with the accumulated power depending on
#'  the number of harmonics
#' *
#' @details Adapted from Claude (2008). pp. 229.
#' @references Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @export
#'
#' @examples
#' library(pliman)
#' pw <- efourier(contours) |> efourier_power()
#'
efourier_power <- function(x,
                           first = TRUE,
                           thresh = c(0.8, 0.85, 0.9, 0.95, 0.99, 0.999),
                           plot = TRUE,
                           ncol = NULL,
                           nrow = NULL){
  if(inherits(x, "efourier_lst")){
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
    power <- lapply(x, efourier_power, first, thresh, plot)
    pwer <-
      do.call(rbind,
              lapply(power, function(x){x$cum_power}))
    pwer$object <- sapply(strsplit(rownames(pwer), "\\."), "[", 1)
    pwer <- pwer[, c(3, 1, 2)]
    rownames(pwer) <- NULL

    min_harm <- do.call(rbind,
                        lapply(power, function(x){x$min_harm}))
    min_harm$object <- names(power)
    min_harm <- min_harm[, c(ncol(min_harm), 1:ncol(min_harm)-1)]
    rownames(min_harm) <- NULL
    return(list(
      cum_power = pwer,
      min_harm = min_harm
    ))
  } else{
    a <- x$an
    b <- x$bn
    c <- x$cn
    d <- x$dn
    Power <- (x$an^2 + x$bn^2 + x$cn^2 + x$dn^2) / 2
    if(isTRUE(first)){
      power <- Power
      cump <- cumsum(power)/sum(power)
      harms <- 1:length(cump)
      if(isTRUE(plot)){
        plot(harms,
             cump,
             pch = 16,
             ylab = "Power",
             xlab = "Number of harmonics")
        lines(harms, cump)
      }
    } else{
      cump <- cumsum(Power[-1])/sum(Power[-1])
      harms <- 2:length(Power)
      if(isTRUE(plot)){
        plot(harms,
             cump,
             pch=16,
             ylab = "Power",
             xlab = "Number of harmonics")
        lines(harms, cump)
      }
    }

    cum_power <-
      data.frame(nharm = harms,
                 cum_power = cump)

    minh <- numeric(length(thresh))
    names(minh) <- paste0(thresh)
    for (i in seq_along(thresh)){
      wi <- which(cump > thresh[i])
      minh[i] <- ifelse(length(wi)==0, NA, min(wi))
    }
    minh <- data.frame(t(minh))
    colnames(minh) <- paste0("p", thresh)
    return(list(
      cum_power = cum_power,
      min_harm = minh
    ))
  }
}

#' Draw shapes based on Fourier coefficients
#'
#' Calculates a 'Fourier elliptical shape' given Fourier coefficients
#'
#' @details
#' \code{efourier_shape} can be used by specifying \code{nharm} and
#' \code{alpha}. The coefficients are then sampled in an uniform distribution
#' \eqn{(-\pi ; \pi)} and this amplitude is then divided by \eqn{harmonicrank ^
#' alpha}. If \code{alpha} is lower than 1, consecutive coefficients will thus
#' increase. See Claude (2008) pp.223 for the maths behind inverse ellipitical
#' Fourier
#'
#' @param an The \eqn{a_n} Fourier coefficients on which to calculate a shape.
#' @param bn The \eqn{b_n} Fourier coefficients on which to calculate a shape.
#' @param cn The \eqn{c_n} Fourier coefficients on which to calculate a shape.
#' @param dn The \eqn{d_n} Fourier coefficients on which to calculate a shape.
#' @param nharm The number of harmonics to use. It must be less than or equal to
#'   the length of `*_n` coefficients.
#' @param n The number of shapes to generate. Defaults to 1. If more than one
#'   shape is used, a list of coordinates is returned.
#' @param npoints The number of points to calculate.
#' @param alpha The power coefficient associated with the (usually decreasing)
#'   amplitude of the Fourier coefficients.
#' @param plot Logical indicating Whether to plot the shape. Defaults to ´TRUE`
#' @return A list with components:
#'  * `x` vector of x-coordrdinates
#'  * `y` vector of y-coordrdinates.
#' @details Adapted from Claude (2008). pp. 223.
#' @references Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' library(pliman)
#' # approximation of the third leaf's perimeter
#' # 4 harmonics
#' image_pliman("potato_leaves.jpg", plot = TRUE)
#'
#' efourier_shape(an = c(-7.34,  1.81,  -1.32, 0.50),
#'                bn = c(-113.88, 21.90, -0.31, -6.14),
#'                cn = c(-147.51, -20.89, 0.66, -14.06),
#'                dn = c(-0.48, 2.36, -4.36, 3.03))
#'
efourier_shape <- function(an = NULL,
                           bn = NULL,
                           cn = NULL,
                           dn = NULL,
                           n = 1,
                           nharm = NULL,
                           npoints = 150,
                           alpha = 4,
                           plot = TRUE) {
  shapes <- list()
  if(n == 1){
    if (is.null(nharm) & is.null(an)){
      nharm <- 10
    }
    if (is.null(nharm) & !is.null(an)){
      nharm <- length(an)
    }
    if (is.null(an)){
      an <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
    }
    if (is.null(bn)){
      bn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
    }
    if (is.null(cn)) {
      cn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
    }
    if (is.null(dn)) {
      dn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
    }
    ef <- list(an = an, bn = bn, cn = cn, dn = dn, a0 = rnorm(1, 10, 30), c0 = rnorm(1, 10, 30))
    shapes <- efourier_inv(ef, nharm = nharm, npoints = npoints) |> unclass()
  } else{
    for(i in 1:n){
      if (is.null(nharm) & is.null(an)){
        nharm <- 10
      }
      if (is.null(nharm) & !is.null(an)){
        nharm <- length(an)
      }
      an <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
      bn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
      cn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
      dn <- runif(nharm, -pi, pi)/(1:nharm) ^ alpha
      ef <- list(an = an, bn = bn, cn = cn, dn = dn, a0 = rnorm(1, 10, 30), c0 = rnorm(1, 10, 30))
      shapes[[i]] <- efourier_inv(ef, nharm = nharm, npoints = npoints) |> unclass()
    }
  }
  if (isTRUE(plot)){
    plot_polygon(shapes)
  }
  return(shapes)
}


