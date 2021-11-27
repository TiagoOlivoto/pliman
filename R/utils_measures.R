#' Utilities for object measures
#'
#' * `get_measures()` computes object measures (area, perimeter, radius) by using
#'either a known resolution (dpi) or an object with known measurements.
#' * `plot_measures()` draws the object measures given in an object to the
#' current plot. The object identification (`"id"`) is drawn by default.
#'
#' @name utils_measures
#' @param object An object computed with [analyze_objects()].
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
#' @param sep Regular expression to manage file names. The function combines in
#'   the `merge` object the object measures (sum of area and mean of all the
#'   other measures) of all images that share the same filename prefix, defined
#'   as the part of the filename preceding the first hyphen (-) or underscore
#'   (_) (no hyphen or underscore is required). For example, the measures of
#'   images named L1-1.jpeg, L1-2.jpeg, and L1-3.jpeg would be combined into a
#'   single image information (L1). This feature allows the user to treat
#'   multiple images as belonging to a single sample, if desired. Defaults to
#'   `sep = "\\_|-"`.
#' @param hjust,vjust A numeric value to adjust the labels horizontally and
#'   vertically. Positive values will move labels to right (hjust) and top
#'   (vjust). Negative values will move the labels to left and bottom,
#'   respectively.
#' @param digits The number of significant figures. Defaults to `2.`
#' @param size The size of the text. Defaults to `0.9`.
#' @param col The color of the text. Defaults to `"white"`.
#' @param verbose If `FALSE`, runs the code silently.
#' @param ... Further arguments passed on to [graphics::text()].
#' @return
#' * For `get_measures()`, if `measure` is informed, the pixel values will be
#' corrected by the value of the known object, given in the unit of the
#' right-hand side of `measure`. If `dpi` is informed, then all the measures
#' will be adjusted to the known `dpi`.
#'
#'    -  If applied to an object of class `anal_obj`, returns a data frame with the
#' object `id` and the (corrected) measures.
#'    - If applied to an object of class `anal_obj_ls`, returns a list of class
#'    `measures_ls`, with two objects: (i) `results`, a data frame containing
#'    the identification of each image (img) and object within each image (id);
#'    and (ii) `summary` a data frame containing the values for each image. If
#'    more than one object is detected in a given image, the number of objects
#'    (`n`), total area (`area_sum`), mean area (`area_mean`) and the standard
#'    deviation of the area (`area_sd`) will be computed. For the other measures
#'    (perimeter and radius), the mean values are presented.
#' * `plot_measures()` returns a `NULL` object, drawing the text according to
#' the x and y coordinates of the objects in `object`.
#' @export
#' @importFrom stats as.formula
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("objects_300dpi.jpg")
#' plot(img)
#' # Image with four objects with a known resolution of 300 dpi
#' # Higher square: 10 x 10 cm
#' # Lower square: 5 x 5 cm
#' # Rectangle: 4 x 2 cm
#' # Circle: 3 cm in diameter
#'
#' # Count the objects using the blue band to segment the image
#' results <-
#'    analyze_objects(img,
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
                         sep = "\\_|-",
                         verbose = TRUE,
                         digits = 3){
  if(is.data.frame(object)){
    if(any(c("area", "perimeter", "radius_mean") %in% colnames(object) == FALSE)){
      stop("Object informed seems to be not an object computed with pliman.")
    }
    res <- object
  }
  if(any(class(object) %in% c("anal_obj", "anal_obj_ls"))){
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
  if(class(object) == "plm_disease"){
    res <- object$shape
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
      res$area_ch <- res$area_ch * px_side^2
      res$perimeter <- res$perimeter * px_side
      res$radius_mean <- res$radius_mean * px_side
      res$radius_min <- res$radius_min * px_side
      res$radius_max <- res$radius_max * px_side
      res$diam_mean <- res$diam_mean * px_side
      res$diam_min <- res$diam_min * px_side
      res$diam_max <- res$diam_max * px_side
      res$major_axis <- res$major_axis * px_side
      res$minor_axis <- res$minor_axis * px_side
    }
    if(var != "area"){
      id_val <- res[which(res$id == id), var]
      px_side <- value / id_val
      res$area <- res$area * px_side^2
      res$area_ch <- res$area_ch * px_side^2
      res$perimeter <- res$perimeter * px_side
      res$radius_mean <- res$radius_mean * px_side
      res$radius_min <- res$radius_min * px_side
      res$radius_max <- res$radius_max * px_side
      res$diam_mean <- res$diam_mean * px_side
      res$diam_min <- res$diam_min * px_side
      res$diam_max <- res$diam_max * px_side
      res$major_axis <- res$major_axis * px_side
      res$minor_axis <- res$minor_axis * px_side
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
    res$area_ch <- res$area_ch * 1/dpc^2
    res$perimeter <- res$perimeter / dpc
    res$radius_mean <- res$radius_mean / dpc
    res$radius_min <- res$radius_min / dpc
    res$radius_max <- res$radius_max / dpc
    res$diam_mean <- res$diam_mean / dpc
    res$diam_min <- res$diam_min / dpc
    res$diam_max <- res$diam_max / dpc
    res$major_axis <- res$major_axis / dpc
    res$minor_axis <- res$minor_axis / dpc
  }
  if("img" %in% names(res)){
    smr <-
      do.call(cbind,
              lapply(5:ncol(res), function(i){
                if(i == 5){
                  n <- aggregate(res[[i]] ~ img, res, length)[[2]]
                  s <- aggregate(res[[i]] ~ img, res, sum, na.rm = TRUE)[2]
                  a <- aggregate(res[[i]] ~ img, res, mean, na.rm = TRUE)[2]
                  d <- aggregate(res[[i]] ~ img, res, sd, na.rm = TRUE)[2]
                  cbind(n, s, a, d)
                } else{
                  aggregate(res[[i]] ~ img, res, mean, na.rm = TRUE)[2]
                }
              })
      )
    names(smr) <- c("n", "area_sum", "area_mean", "area_sd",  names(res[6:ncol(res)]))
    smr$img <- unique(res$img)
    smr <- smr[,c(ncol(smr), 1:ncol(smr)-1)]
    smr$area_sd[is.na(smr$area_sd)] <- 0
    merg <- smr
    merg$img = sapply(strsplit(as.character(merg$img), sep), "[", 1)
    mergt <-
      do.call(cbind,
              lapply(2:ncol(merg), function(i){
                if(i %in% 2:3){
                  aggregate(merg[[i]] ~ img, merg, sum, na.rm = TRUE)[2]
                } else{
                  aggregate(merg[[i]] ~ img, merg, mean, na.rm = TRUE)[2]
                }
              })
      )
    mergt$img <- unique(merg$img)
    mergt <- mergt[,c(ncol(mergt), 1:ncol(mergt)-1)]
    names(mergt) <- names(smr)
    smr[,3:ncol(smr)] <- apply(smr[,3:ncol(smr)], 2, round, digits)
    res[,3:ncol(res)] <- apply(res[,3:ncol(res)], 2, round, digits)
    rownames(res) <- NULL
    mergt[,3:ncol(mergt)] <- apply(mergt[,3:ncol(mergt)], 2, round, digits)
    out <-
      list(results = res,
           summary = smr,
           merge = mergt)
    class(out) <- c("measures_ls")
    return(out)
  } else{
    res[,2:ncol(res)] <- apply(res[,2:ncol(res)], 2, round, digits)
    class(res) <- c("data.frame", "measures")
    return(res)
  }
}

#' @name utils_measures
#' @export
plot_measures <- function(object,
                          id = NULL,
                          measure = "id",
                          hjust = NULL,
                          vjust = NULL,
                          digits = 2,
                          size = 0.9,
                          col = "white",
                          ...){
  if("measures"  %in% class(object)){
    object <- object
  } else if(class(object) == "anal_obj"){
    index <- object$object_index
    object <- object$results
  } else if(class(object) == "objects_rgb"){
    object <- object$objects
  } else if(class(object) == "plm_disease"){
    object <- object$shape
  } else{
    stop("Object of ivalid class.")
  }
  if(is.null(id)){
    id <- object$id
  } else{
    id <- id
  }
  object <- object[which(object$id %in% id), ]
  if(measure %in% colnames(object)){
    hjust <- ifelse(is.null(hjust), 0, hjust)
    vjust <- ifelse(is.null(vjust), 0, vjust)
    text(x = object[,2] + hjust,
         y = object[,3] - vjust,
         labels = round(object[, which(colnames(object) == measure)], digits),
         col = col,
         cex = size,
         ...)
  } else{
    if(!is.null(index)){
      measures <- colnames(index)
      if(!measure %in% measures){
        stop("'measure' must be one of {", paste(c(colnames(object), measures), collapse = ", "),"}.", call. = FALSE)
      }
      text(x = object[,2],
           y = object[,3],
           labels = round(index[object$id , which(colnames(index) == measure)], digits),
           col = col,
           cex = size,
           ...)
    } else{
      stop("'measure' must be one of {", paste(colnames(object), collapse = ", "),"}.", call. = FALSE)
    }
  }
}



#' Summary an object index
#'
#' Performs a report of the index between and within objects when `object_index`
#' argument is used in `analyze_objects()`. By using a cut point, the number and
#' proportion of objects with mean value of `index` bellow and above `cut_point`
#' are returned. Additionaly, the number and proportion of pixels bellow and
#' above the cutpoint is shown for each object (id).
#'
#' @param object An object computed with [analyze_objects()].
#' @param index The index desired, e.g., `"B"`. Note that these value must match
#'   the index(es) used in the argument `object_index` of `analyze_objects()`.
#' @param cut_point The cut point.
#' @param select_higher If `FALSE` (default) selects the objects with `index`
#'   smaller than the `cut_point`. Use `select_higher = TRUE` to select the
#'   objects with `index` higher than `cut_point`.
#'
#' @return A list with the following elements:
#' * `ids` The identification of selected objects.
#' * `between_id` A data frame with the following columns
#'    - `n` The number of objects.
#'    - `nsel` The number of selected objects.
#'    - `prop` The proportion of objects selected.
#'    - `mean_index_sel`, and `mean_index_nsel` The mean value of `index` for the
#' selected and non-selected objects, respectively.
#' * `within_id` A data frame with the following columns
#'    - `id` The object identification
#'    - `n_less` The number of pixels with values lesser than or equal to
#'    `cut_point`.
#'    - `n_greater` The number of pixels with values greater than `cut_point`.
#'    - `less_ratio` The proportion of pixels with values lesser than or equal to
#'    `cut_point`.
#'    - `greater_ratio` The proportion of pixels with values greater than
#'    `cut_point`.
#' @importFrom stats setNames
#' @export
#' @name summary_index
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examples
#' library(pliman)
#' soy <- image_pliman("soy_green.jpg")
#' anal <- analyze_objects(soy, object_index = "G")
#' plot_measures(anal, measure = "G")
#'
#' summary_index(anal, index = "G", cut_point = 0.5)
summary_index <- function(object,
                         index,
                         cut_point,
                         select_higher = FALSE){
  if(is.null(object$object_index)){
    stop("'object' was not computed using the `object_index` argument.")
  }
  coords <- object$results[2:3]
  obj_in <- object$object_index

  if(isFALSE(select_higher)){
    ids <- which(obj_in[[index]] <= cut_point)
  } else{
    ids <- which(obj_in[[index]] >= cut_point)
  }

  temp <- object$object_rgb
  indexes <-
    do.call(rbind,
            lapply(  by(temp,
                        INDICES = temp$id,
                        FUN = function(x){
                          x[[index]] <= cut_point
                        }
            ), data.frame)
    )
  res <-
    transform(temp,
              threshold = ifelse(indexes[[1]] == TRUE, "less", "greater"),
              n = 1:nrow(indexes)) %>%
    aggregate(n ~ id + threshold, FUN = length, data = .) %>%
    reshape(direction = "wide",
            timevar = "threshold",
            idvar = "id",
            v.names = "n") %>%
    as.data.frame() %>%
    setNames(c("id", "n_greater", "n_less")) %>%
    transform(less_ratio = round(n_less / (n_less  + n_greater), 3),
              greater_ratio = 1 - round(n_less / (n_less  + n_greater), 3))


  list(ids = ids,
       between_id = data.frame(
         n = nrow(obj_in),
         nsel = length(ids),
         prop = length(ids) / nrow(obj_in),
         mean_index_sel = mean(obj_in[[index]][ids]),
         mean_index_nsel = mean(obj_in[[index]][!obj_in$id %in% ids])
       ),
       within_id = cbind(coords, res)[, c(3, 1, 2, 5, 4, 6, 7)]
  )
}

