#' Utilities for object measures
#'
#'* `get_measures()` computes object measures (area, perimeter, radius) by using
#'either a known resolution (dpi) or an object with known measurements.
#' * `plot_measures()` draws the object measures given in an object to the
#' current plot. The object identification (`"id"`) is drawn by default.
#'
#' @name utils_measures
#' @param object An object computed with [count_objects()] or [leaf_area()].
#' @param dpi A known resolution of the image in DPI (dots per inch).
#' @param id An object in the image to indicate a known value.
#' @param measure For `plot_measures()`, a character string; for
#'   `get_measures()`, a two-sided formula, e.g., `measure = area ~ 100`
#'   indicating the known value of object `id`. The right-hand side is the known
#'   value and the left-hand side can be one of the following.
#' * `area` The known area of the object.
#' * `perimeter` The known perimeter of the object.
#' * `radius_mean` The known radius of the object.
#' * `radius_min` The known minimum radius of the object. If the object is a
#' square, then the `radius_min` of such object will be `L/2` where `L` is the
#' length of the square side.
#' * `radius_max` The known maximum radius of the object. If the object is a
#' square, then the `radius_max` of such object according to the Pythagorean
#' theorem will be `L x sqrt(2) / 2` where `L` is the length of the square side.
#' @param digits The number of significant figures. Defaults to `2.`
#' @param size The size of the text. Defaults to `0.9`.
#' @param col The color of the text. Defaults to `"white"`.
#' @param verbose If `FALSE`, runs the code silently.
#' @param ... Further arguments passed on to [graphics::text()].
#' @return
#' * `get_measures()` returns a data frame with the object `id` and the
#' measures. If `measure` is informed, the pixel values will be corrected by the
#' value of the known object, given in the unit of the right-hand side of
#' `measure`. If `dpi` is informed, then all the measures will be adjusted to
#' the known `dpi`.
#' * `plot_measures()` returns a `NULL` object, drawing the text according to
#' the x and y coordinates of the objects in `object`.
#' @export
#' @importFrom stats as.formula
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(image_pliman("objects_300dpi.jpg"))
#' image_show(img)
#' # Image with four objects with a known resolution of 300 dpi
#' # Higher square: 10 x 10 cm
#' # Lower square: 5 x 5 cm
#' # Rectangle: 4 x 2 cm
#' # Circle: 3 cm in diameter
#'
#' # Count the objects using the blue band to segment the image
#' results <-
#'    count_objects(img,
#'                  index = "B")
#' plot_measures(results, measure = "id")
#'
#' # Get object measures by declaring the known resolution in dots per inch
#' (measures <- get_measures(results, dpi = 300))
#'
#' # Calculated diagonal of the object 1
#' # 10 * sqrt(2) = 14.14
#'
#' # Observed diagonal of the object 1
#' measures[1, "radius_max"] * 2
#'
#'
#' # Get object measures by declaring the known area of object 1
#' get_measures(results,
#'              id = 1,
#'              area ~ 100)
#'}
get_measures <- function(object,
                         id = NULL,
                         measure = NULL,
                         dpi = NULL,
                         verbose = TRUE,
                         digits = 3){
  if(is.data.frame(object)){
    if(any(c("area", "perimeter", "radius_mean") %in% colnames(object) == FALSE)){
      stop("Object informed seems to be not an object computed with pliman.")
    }
    res <- object
  }
  if(any(class(object) == "plm_count")){
    res <- object$results
  }
  if(any(class(object) == "plm_la")){
    obj <- as.numeric(class(object)[[3]])
    res <- object[-obj,-ncol(object)]
    if(verbose == TRUE){
      cat("-----------------------------------------\n")
      cat(paste0("Total leaf area  : ", round(sum(object$area), 3)), "\n")
      cat(paste0("Average leaf area: ", round(mean(object$area), 3)), "\n")
      cat("-----------------------------------------\n")
    }
  }
  if(any(class(object) == "objects_rgb")){
    res <- object[["objects"]]
  }
  if(!is.null(id) & !is.null(dpi)){
    stop("Only one of 'dpi' or 'id' can be used.", call. = FALSE)
  }
  if(!is.null(id) & is.null(measure) ){
    stop("'measure' must be informed.", call. = FALSE)
  }
  if(!is.null(id)){
    if(class(measure) != "formula"){
      stop("'measure' must be a two-sided formula, e.g., 'area ~ 25'.")
    }
    terms <- as.formula(measure)
    var <- as.character(terms[[2]])
    value <- as.numeric(terms[[3]])
    measures <- c("area", "perimeter", "radius_mean", "radius_min",  "radius_max", "radius_ratio")
    if(!var %in% measures){
      stop("The left-hand side of 'measure' must be one of ", paste(measures, collapse = ", "), call. = FALSE)
    }
    if(var == "area"){
      id_val <- res[which(res$id == id), var]
      px_side <- sqrt(value / id_val)
      values <- res[, var]
      corrected <- values * value / id_val
      res$area <- corrected
      res$perimeter <- res$perimeter * px_side
      res$radius_mean <- res$radius_mean * px_side
      res$radius_min <- res$radius_min * px_side
      res$radius_max <- res$radius_max * px_side
    }
    if(var != "area"){
      id_val <- res[which(res$id == id), var]
      px_side <- value / id_val
      res$area <- res$area * px_side^2
      res$perimeter <- res$perimeter * px_side
      res$radius_mean <- res$radius_mean * px_side
      res$radius_min <- res$radius_min * px_side
      res$radius_max <- res$radius_max * px_side
    }
    if(verbose == TRUE){
      cat("-----------------------------------------\n")
      cat(paste0("measures corrected with:\nobject id: ", id, "\n", var,
                 "     : ",  value, "\n"))
      cat("-----------------------------------------\n")
      cat(paste0("Total    : ", round(sum(res[, var]), 3)), "\n")
      cat(paste0("Average  : ", round(mean(res[, var]), 3)), "\n")
      cat("-----------------------------------------\n")
    }
  }
  if(!is.null(dpi)){
    dpc <- dpi * 1 / 2.54
    res$area <- res$area * 1/dpc^2
    res$perimeter <- res$perimeter / dpc
    res$radius_mean <- res$radius_mean / dpc
    res$radius_min <- res$radius_min / dpc
    res$radius_max <- res$radius_max / dpc
  }
  res[,1:10] <- apply(res[,1:10], 2, round, digits)
  class(res) <- c("data.frame", "plm_measures")
  return(res)
}

#' @name utils_measures
#' @export
plot_measures <- function(object,
                          measure = "id",
                          digits = 2,
                          size = 0.9,
                          col = "white",
                          ...){
  if("plm_measures"  %in% class(object)){
    object <- object
  } else if(class(object) == "plm_count"){
    object <- object$results
  } else if(class(object) == "objects_rgb"){
    object <- cbind(object[["objects"]], index = object$indexes$index)
  } else{
    stop("Object of ivalid class.")
  }
  measures <- c("id", "area", "perimeter", "radius_mean", "radius_min", "radius_max", "radius_ratio", "index")
  if(!measure %in% measures){
    stop("'measure' must be one of the", paste(measures, collapse = ", "), call. = FALSE)
  }
  text(x = object[,2],
       y = object[,3],
       labels = round(object[, which(colnames(object) == measure)], digits),
       col = col,
       cex = size,
       ...)
}

