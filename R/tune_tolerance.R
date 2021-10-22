#' Tune tolerance parameter
#' @description
#' Provides options for tunning `tolerance` parameter utilized in
#' `[analyze_objects()]` in two ways:
#'  * Declaring the `actual` argument, an iterative algorithm will compute the
#'  first analysis and sequentially increase the parameter `tolerance` if the
#'  computed number of objects is greater than `actual` or reduce the parameter
#'  `tolerance` if the computed number of objects is less than `actual`. If the
#'  algorithm did not converge up to `maxiter` is reached, users can change the
#'  default `extension` value.
#'
#'  * The second way is to create a `grid` with `tolerance` and `extension`
#'  values. When `grid` is informed, all combinations (made by
#'  [base::expand.grid()]) are tested and the residual from `actual` value is
#'  plotted. Users can than find a better combination of parameters to use in
#'  `analyze_objects()`.
#'
#'
#' @param img The image to be analyzed.
#' @param actual The actual number of objects.
#' @param start_tol An starting value for tolerance. Defaults to `1`.
#' @param extension The extension value. Defaults to `1`.
#' @param grid A list with a numeric sequence for `tolerance` and `extension`
#'   values. When grid is informed, all combinations are tested and the residual
#'   from `actual` value is plotted.
#' @param maxiter The maximum number of iterations. Default to 200.
#' @param index,my_index A character value specifying the target mode for
#'   conversion to binary image when `foreground` and `background` are not
#'   declared. Defaults to `"NB"` (normalized blue). See [image_index()] for
#'   more details.
#' @param plot Logical. If `TRUE` (default) generates a plot showing the
#'   results.
#' @param fill_hull Fill holes in the binary image? Defaults to `FALSE`. This is
#'   useful to fill holes in objects that have portions with a color similar to
#'   the background. IMPORTANT: Objects touching each other can be combined into
#'   one single object, which may underestimate the number of objects in an
#'   image.
#' @param filter Performs median filtering after image processing? defaults to
#'   `FALSE`. See more at [image_filter()].
#' @param invert Inverts the binary image, if desired. This is useful to process
#'   images with black background. Defaults to `FALSE`.
#' @param workers The number of multiple sections to be used in the computation.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @export
#' @importFrom grDevices terrain.colors
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' library(pliman)
#' if(FALSE){
#'  img <- image_pliman("soybean_touch.jpg")
#'
#'  # starts with tolerance = 1
#'  tune_tolerance(img, actual = 30)
#'
#'  # Using a grid of tolerance and extension values
#'  tune2 <-
#'    tune_tolerance(img,
#'                   actual = 30,
#'                   grid = list(tolerance = seq(0, 5, by = 0.1),
#'                               extension = 1:4), plot = TRUE)
#'}
tune_tolerance <- function(img,
                           actual,
                           start_tol = NULL,
                           extension = NULL,
                           grid = NULL,
                           maxiter = 200,
                           index = "NB",
                           my_index = NULL,
                           plot = TRUE,
                           fill_hull = FALSE,
                           filter = FALSE,
                           invert = FALSE,
                           workers = NULL,
                           verbose = TRUE){
  check_ebi()
  get_num <- function(tolerance, extension){
    bin <-
      image_binary(img,
                   index = index,
                   my_index = my_index,
                   invert = invert,
                   fill_hull = fill_hull,
                   resize = FALSE,
                   show_image = FALSE)[[1]]
    nmask <- EBImage::watershed(EBImage::distmap(bin),
                                tolerance = tolerance,
                                ext = extension)
    shape <- EBImage::computeFeatures.shape(nmask)
    n <- nrow(shape[shape[,1] > mean(shape[,1]) * 0.1, ])
    return(n)
  }
  if(is.null(grid)){
    npix <- npixels(img)
    extension <- ifelse(is.null(extension),
                        ifelse(npix < 560000, 1, 2),
                        extension)
    tolerance <- ifelse(is.null(start_tol), 1, start_tol)
    n <- get_num(extension = extension, tolerance = tolerance)
    dif <- n - actual
    inic <- 0
    while(dif != 0){
      if(dif < 0){
        val_iter <- abs(tolerance * log(exp(abs(dif))*npix/10000, 1e20))
        tolerance <- ifelse(tolerance < val_iter, tolerance * 0.5, tolerance - tolerance * val_iter)
      } else{
        val_iter <- abs(tolerance * log(exp(abs(dif))*npix/10000, 1e20))
        tolerance <- ifelse(tolerance < val_iter, tolerance * 0.5, tolerance + tolerance * val_iter)
      }
      n <- get_num(extension = extension, tolerance = tolerance)
      dif <- n - actual
      inic <- inic + 1
      if(isTRUE(verbose)){
        cat("Iteration: ", inic,
            "| Difference:", dif,"\r")
      }
      if(inic >= maxiter){
        warning("\nMaximum number of iterations reached before convergence.\n", call. = FALSE)
        return(list(tolerance = tolerance,
                    extension = extension))
      }
    }
    if(isTRUE(verbose)){
      message("\nConvergence achieved!\n")
      cat("Tolerance: ", round(tolerance, 4),"\n")
      cat("Extension: ", round(extension, 4),"\n")
    }
    invisible(list(tolerance = tolerance,
                   extension = extension))
  } else{
    comb <- expand.grid(grid[["tolerance"]], grid[["extension"]])
    names(comb) <- c("tolerance", "extension")
    nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
    clust <- makeCluster(nworkers)
    clusterExport(clust,
                  varlist = c("img", "index",
                              "get_num", "comb",
                              "invert",
                              "image_binary"),
                  envir=environment())
    on.exit(stopCluster(clust))
    message("Tunning parameters using ",nworkers, " sections. Please wait.")
    n <-
      as.numeric(
        parLapply(clust, 1:nrow(comb),
                  function(i){
                    get_num(tol = comb[i, 1], ext = comb[i, 2])
                  })
      )
    comb$n <- as.numeric(n)
    comb <- transform(comb, diff = n - actual)
    if(actual %in% n){
      candidate <- comb[comb$n == actual,]
    } else{
      candidate <- comb[which.min(abs(comb$n - actual)),]
    }
    if(isTRUE(plot)){
      gg <-
        levelplot(diff ~ tolerance * extension,
                           data = comb,
                           col.regions = terrain.colors(300),
                           colorkey = list(title = "Difference",
                                           interpolate = TRUE,
                                           raster = TRUE))
      plot(gg)
    } else{
      gg <- NULL
    }
    return(list(results = comb,
                candidate = candidate,
                plot =  gg))
  }
}
