#' Utilities for object measures
#'
#'* `get_measures()` computes object measures (area, perimeter, radius) by using
#'either a known resolution (dpi) or an object with known measurements.
#' * `plot_measures()` draws the object measures given in an object to the
#' current plot. The object identification (`"id"`) is drawn by default.
#'
#' @name utils_measures
#' @param object An object computed with [analyze_objects()].
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
#' @param id An object in the image to indicate a known value.
#' @param dpi A known resolution of the image in DPI (dots per inch).
#' @param sep Regular expression to manage file names. The function combines in
#'   the `merge` object the object measures (sum of area and mean of all the
#'   other measures) of all images that share the same filename prefix, defined
#'   as the part of the filename preceding the first hyphen (-) or underscore
#'   (_) (no hyphen or underscore is required). For example, the measures of
#'   images named `L1-1.jpeg`, `L1-2.jpeg`, and `L1-3.jpeg` would be combined
#'   into a single image information (L1). This feature allows the user to treat
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
#' right-hand side of `meae`. If `dpi` is informed, then all the measures
#' will be adjusted to the knosurwn `dpi`.
#'
#' -  If applied to an object of class `anal_obj`, returns a data frame with the
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
#'                  index = "B",
#'                  lower_noise = 0.1)
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
#'
#'
get_measures <- function(object,
                         measure = NULL,
                         id = NULL,
                         dpi = NULL,
                         sep = "\\_|-",
                         verbose = TRUE,
                         digits = 5){
  if(is.data.frame(object)){
    if(any(c("area", "perimeter", "radius_mean") %in% colnames(object) == FALSE)){
      stop("Object informed seems to be not an object computed with pliman.")
    }
    res <- object
  }
  if(any(class(object) %in% c("anal_obj", "anal_obj_ls"))){
    res <- object$results
  }
  if(any(inherits(object, "objects_rgb"))){
    res <- object[["objects"]]
  }
  if(class(object)  %in%  c("plm_disease", "plm_disease_byl")){
    res <- object$shape
  }
  if(!is.null(id) & !is.null(dpi)){
    stop("Only one of 'dpi' or 'id' can be used.", call. = FALSE)
  }
  if(!is.null(id) & is.null(measure) ){
    stop("'measure' must be informed.", call. = FALSE)
  }
  ncols <- ifelse(class(object)  %in%  c("plm_disease", "plm_disease_byl"), 16, 18)

  if(!is.null(id)){
    if(!inherits(measure, "formula")){
      stop("'measure' must be a two-sided formula, e.g., 'area ~ 25'.")
    }
    terms <- as.formula(measure)
    var <- as.character(terms[[2]])
    if(exists(as.character(terms[[3]]), envir = parent.frame())){
      value <- eval(terms[[3]], envir = parent.frame())
    } else{
      value <- as.numeric(terms[[3]])
    }
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
      if(inherits(object, "plm_disease_byl")){
        res[7:18] <- apply(res[7:18], 2, function(x){
          x * px_side
        })
      } else{
        if(inherits(object, "plm_disease")){
          res[5:ncols] <- apply(res[5:ncols], 2, function(x){
            x * px_side
          })
        } else{
          res[6:ncols] <- apply(res[6:ncols], 2, function(x){
            x * px_side
          })
        }
      }

    }
    if(var != "area"){
      id_val <- res[which(res$id == id), var]
      px_side <- value / id_val
      res$area <- res$area * px_side^2
      res$area_ch <- res$area_ch * px_side^2
      if(inherits(object, "plm_disease_byl")){
        res[7:18] <- apply(res[7:18], 2, function(x){
          x * px_side
        })
      } else{
        if(inherits(object, "plm_disease")){
          res[5:ncols] <- apply(res[5:ncols], 2, function(x){
            x * px_side
          })
        } else{
          res[6:ncols] <- apply(res[6:ncols], 2, function(x){
            x * px_side
          })
        }
      }
    }
    res <- res[which(res$id != id),]
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
    if(!inherits(object, c("plm_disease", "plm_disease_byl"))){
      res$area_ch <- res$area_ch * 1/dpc^2
    }
    if(inherits(object, "plm_disease_byl")){
      res[7:18] <- apply(res[7:18], 2, pixels_to_cm, dpi = dpi)
    } else{
      if("img" %in% colnames(res)){
        if(inherits(object, "plm_disease")){
          res[6:(ncols + 1)] <- apply(res[6:(ncols + 1)], 2, pixels_to_cm, dpi = dpi)
        } else{
          res[7:(ncols + 1)] <- apply(res[7:(ncols + 1)], 2, pixels_to_cm, dpi = dpi)
        }
      } else{
        if(inherits(object, "plm_disease")){
          res[5:ncols] <- apply(res[5:ncols], 2, pixels_to_cm, dpi = dpi)
        } else{
          res[6:ncols] <- apply(res[6:ncols], 2, pixels_to_cm, dpi = dpi)
        }
      }
    }
  }

  if("img" %in% names(res)){
    if(!inherits(object, "plm_disease_byl") & !inherits(object, "anal_obj")){

      # bind object_index, if it exists
      if(!is.null(object$object_index)){
        if(ncol(object$object_index) < 4){
          nam_res <- colnames(res)
          nam_ind <- colnames(object$object_index)[3]
          res <- cbind(res, object$object_index[, 3])
          colnames(res) <- c(nam_res, nam_ind)
        } else{
          res <- cbind(res, object$object_index[, -c(1:2)])
        }
      }

      # bind efourier coefficients, if it exists
      if(!is.null(object$efourier)){
        res <- cbind(res, object$efourier_norm[, -c(1:2)])
      }

      # bind apex and base angles if it exists
      if(!is.null(object$angles)){
        res <- cbind(res, object$angles[, -c(1:2)])
      }
      # bind perimeter complexity value if it exists
      if(!is.null(object$pcv)){
        res <- cbind(res, pcv = object[["pcv"]][, 2])
      }
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
      smr$img <- paste0("img", unique(res$img))
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
    }

    if(inherits(object, "plm_disease_byl")){
      smr <-
        do.call(cbind,
                lapply(6:ncol(res), function(i){
                  if(i == 6){
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
      names(smr) <- c("n", "area_sum", "area_mean", "area_sd",  names(res[7:ncol(res)]))
      smr$img <- paste0("img", unique(res$img))
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

    }

    if(inherits(object, "anal_obj")){
      index <- object$object_index
      shapefiles <- object$shapefiles$shapefiles
      coords <-
        do.call(rbind,
                lapply(shapefiles, function(x){
                  data.frame(x = mean(x$x[-1]), y = mean(x$y[-1]))
                }))
      if(!is.null(index)){
        index$img <- as.numeric(gsub(pattern = "shp", x = index$img, replacement = ""))
        aggr <-
          do.call(cbind,
                  lapply(3:ncol(index), function(i){
                    aggregate(index[[i]] ~ img, index, mean, na.rm = TRUE)[2]
                  })
          )
        names(aggr) <- c(names(index[3:ncol(index)]))
        aggr$img <- paste0("obj", unique(index$img))
        aggr <- aggr[,c(ncol(aggr), 1:(ncol(aggr)-1))]
        aggr$x <- coords$x
        aggr$y <- coords$y
        aggr <- aggr[, c(c("img", "x", "y"), setdiff(colnames(aggr), c("img", "x", "y")))]
      } else{
        aggr <- NULL
      }
      res_img <- res$img
      res$img <- as.numeric(gsub(pattern = "shp", x = res$img, replacement = ""))

      # bind object_index, if it exists
      if(!is.null(object$object_index)){
        if(ncol(object$object_index) < 4){
          nam_res <- colnames(res)
          nam_ind <- colnames(object$object_index)[3]
          res <- cbind(res, object$object_index[, 3])
          colnames(res) <- c(nam_res, nam_ind)
        } else{
          res <- cbind(res, object$object_index[, -c(1:2)])
        }
      }

      # bind efourier coefficients, if it exists
      if(!is.null(object$efourier)){
        res <- cbind(res, object$efourier_norm[, -c(1:2)])
      }

      # bind apex and base angles if it exists
      if(!is.null(object$angles)){
        res <- cbind(res, object$angles[, -c(1:2)])
      }
      # bind width_at if it exists
      if(!is.null(object$width_at)){
        res <- cbind(res, object$width_at[, -c(1:2)])
      }
      # bind perimeter complexity value if it exists
      if(!is.null(object$pcv)){
        res <- cbind(res, pcv = object$pcv)
      }

      smr <-
        do.call(cbind,
                lapply(5:ncol(res), function(i){
                  if(i  %in% c(5, 6, 35)){
                    if(i == 5){
                      n <- aggregate(res[[i]] ~ img, res, length)[[2]]
                      a <- aggregate(res[[i]] ~ img, res, sum, na.rm = TRUE)[2]
                      cbind(n, a)
                    } else{
                      aggregate(res[[i]] ~ img, res, sum, na.rm = TRUE)[2]
                    }
                  } else{
                    aggregate(res[[i]] ~ img, res, mean, na.rm = TRUE)[2]
                  }
                })
        )
      names(smr) <- c("n", "area", names(res[6:ncol(res)]))
      smr$img <- paste0("obj", unique(res$img))
      res$img <- res_img
      smr <- smr[,c(ncol(smr), 1:ncol(smr)-1)]
      smr[,3:ncol(smr)] <- apply(smr[,3:ncol(smr)], 2, round, digits)
      res[,3:ncol(res)] <- apply(res[,3:ncol(res)], 2, round, digits)
      smr$x <- coords$x
      smr$y <- coords$y
      smr <- smr[, c(c("img", "x", "y", "n"), setdiff(colnames(smr), c("img", "x", "y", "n")))]
      rownames(res) <- NULL

      class(res) <- c("data.frame", "measures")
      class(smr) <- c("data.frame", "measures")
      out <-
        list(results = res,
             summary = smr,
             index = aggr)
    }
    class(out) <- c("measures_ls")
    invisible(out)
  } else{

    # bind object_index, if it exists
    if(!is.null(object$object_index)){
      if(ncol(object$object_index) < 3){
        nam_res <- colnames(res)
        nam_ind <- colnames(object$object_index)[2]
        res <- cbind(res, object$object_index[, 2])
        colnames(res) <- c(nam_res, nam_ind)
      } else{
        res <- cbind(res, object$object_index[, -1])
      }
    }
    # bind efourier coefficients, if it exists
    if(!is.null(object$efourier)){
      res <- cbind(res, object$efourier_norm[, -1])
    }
    # bind apex and base angles if it exists
    if(!is.null(object$angles)){
      res <- cbind(res, object$angles[, -1])
    }

    # bind width_at if it exists
    if(!is.null(object$width_at)){
      res <- cbind(res, object$width_at[, -1])
    }

    # bind perimeter complexity value if it exists
    if(!is.null(object$pcv)){
      res <- cbind(res, pcv = object[["pcv"]])
    }

    res <- round_cols(res, digits = digits)
    class(res) <- c("data.frame", "measures")
    invisible(res)
  }
}


#' @name utils_measures
#' @export
plot_measures <- function(object,
                          measure = "id",
                          id = NULL,
                          hjust = NULL,
                          vjust = NULL,
                          digits = 2,
                          size = 0.9,
                          col = "white",
                          ...){
  if("shapefiles" %in% names(object)){
    if(inherits(object, "plm_disease_byl")){
      object <- object$severity
    } else{
      meas <- get_measures(object)$summary
      index <- object$object_index
      shapefiles <- object$shapefiles$shapefiles
      coords <-
        do.call(rbind,
                lapply(shapefiles, function(x){
                  data.frame(x = mean(x$x[-1]), y = mean(x$y[-1]))
                }))
      object <- cbind(shp = meas[,1], coords, meas[,2:ncol(meas)])

      if(!is.null(index)){

        index$img <- as.numeric(gsub(pattern = "shp", x = index$img, replacement = ""))
        aggr <-
          do.call(cbind,
                  lapply(3:ncol(index), function(i){
                    aggregate(index[[i]] ~ img, index, mean, na.rm = TRUE)[2]
                  })
          )
        names(aggr) <- c(names(index[3:ncol(index)]))
        aggr$img <- paste0("obj", unique(index$img))
        aggr <- aggr[,c(ncol(aggr), 1:(ncol(aggr)-1))]
        aggr <- cbind(obj = aggr[,1], coords, data.frame(aggr[,2:ncol(aggr)]))
        colnames(aggr) <- c("obj", "x", "y", colnames(index)[3:ncol(index)])
        index <- aggr
      } else{
        index <- NULL
      }
    }


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
        text(x = index[,2],
             y = index[,3],
             labels = round(index[, which(colnames(index) == measure)], digits),
             col = col,
             cex = size,
             ...)
      } else{
        stop("'measure' must be one of {", paste(colnames(object), collapse = ", "),"}.", call. = FALSE)
      }
    }

  } else{
    if("measures"  %in% class(object)){
      object <- object
    } else if(inherits(object, "anal_obj")){
      index <- object$object_index
      if(!is.null(object$pcv)){
        object <- cbind(object$results, pcv = object$pcv)
      } else{
        object <- object$results
      }
    } else if(inherits(object, "objects_rgb")){
      object <- object$objects
    } else if(inherits(object, "plm_disease")){
      object <- object$shape
    } else if(inherits(object, "measures_ls")){
      index <- object$index
      colnames(index)[1] <- "id"
      object <- object$summary
      colnames(object)[1] <- "id"
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
             labels = round(index[which(index$id %in% object$id), which(colnames(index) == measure)], digits),
             col = col,
             cex = size,
             ...)
      } else{
        stop("'measure' must be one of {", paste(colnames(object), collapse = ", "),"}.", call. = FALSE)
      }
    }
  }
}


#' Plot length and width lines on objects
#'
#' This function plots the length and width lines given an `object` computed
#' with [analyze_objects()]. The function does not call `plot.new`, so it must
#' be called after an image is plotted. This can be done either using, e.g.,
#' `plot(img)`, or `analyze_objects(..., plot = TRUE)`.
#'
#' @param object An object computed with [analyze_objects()].
#' @param col_length The color of the length line. Default is `"red"`.
#' @param col_width The color of the width line. Default is `"green"`.
#' @param lwd_length The line width of the length line. Default is 2.
#' @param lwd_width The line width of the width line. Default is 2.
#'
#' @details This function takes an object computed with [analyze_objects()] and
#'   plots the length and width lines of each object onto an image. The length
#'   and width lines are calculated based on the position and orientation of the
#'   object, and are plotted using the specified colors and line widths.
#'
#' @importFrom graphics lines
#' @export
#'
#' @examples
#' img <- image_pliman("flax_leaves.jpg")
#' res <- analyze_objects(img, watershed = FALSE, show_contour = FALSE)
#' plot_lw(res)
plot_lw <- function(object,
                    col_length = "red",
                    col_width = "green",
                    lwd_length = 2,
                    lwd_width = 2){
  if(inherits(object, "anal_obj")){
    rest <- object$results
  } else{
    rest <- object
  }
  if(!all(c("x", "y", "length", "width", "theta") %in% colnames(rest))){
    stop("`object` must be an object computed with `analyze_objects() or a data.frame with the columns `x`, `y`, `length`, `width`, and `theta`", call. = FALSE)
  }
  xc <- rest$x
  yc <- rest$y
  length <- rest$length
  width <- rest$width
  theta <- rest$theta

  theta_degrees <- theta * 180 / pi
  # Calculate the endpoints of the length line
  xls <- xc - (length/2)*cos(theta)
  yls <- yc - (length/2)*sin(theta)
  xle <- xc + (length/2)*cos(theta)
  yle <- yc + (length/2)*sin(theta)
  # Calculate the endpoints of the width line
  xws <- xc - (width/2)*cos(theta + pi/2)
  yws <- yc - (width/2)*sin(theta + pi/2)
  xwe <- xc + (width/2)*cos(theta + pi/2)
  ywe <- yc + (width/2)*sin(theta + pi/2)
  # Plot the lines
  segments(xws, yws, xwe,  ywe, col = col_width, lwd = lwd_width)
  segments(xls, yls, xle, yle, col = col_length, lwd = lwd_length)
}



#' Summary an object index
#'
#' If more than one index is available, the function performs a Principal
#' Component Analysis and produces a plot showing the contribution of the
#' indexes to the PC1 (see [pca()]). If an index is declared in
#' `index` and a cut point in `cut_point`, the number and proportion of objects
#' with mean value of `index` bellow and above `cut_point` are returned.
#' Additionaly, the number and proportion of pixels bellow and above the
#' cutpoint is shown for each object (id).
#'
#' @param object An object computed with [analyze_objects()].
#' @param index The index desired, e.g., `"B"`. Note that these value must match
#'   the index(es) used in the argument `object_index` of `analyze_objects()`.
#' @param cut_point The cut point.
#' @param select_higher If `FALSE` (default) selects the objects with `index`
#'   smaller than the `cut_point`. Use `select_higher = TRUE` to select the
#'   objects with `index` higher than `cut_point`.
#' @param plot Shows the contribution plot when more than one index is
#'   available? Defaults to `TRUE`.
#' @param type The type of plot to produce. Defaults to `"var"`. See more at
#'   [get_biplot()].
#' @param ... Further arguments passed on to [get_biplot()].
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
#' * `pca_res` An object computed with [pca()]
#' @importFrom stats prcomp setNames
#' @importFrom graphics abline barplot
#' @export
#' @name summary_index
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examples
#' library(pliman)
#' soy <- image_pliman("soy_green.jpg")
#' anal <- analyze_objects(soy, object_index = "G", pixel_level_index = TRUE)
#' plot_measures(anal, measure = "G")
#'
#' summary_index(anal, index = "G", cut_point = 0.5)
summary_index <- function(object,
                          index = NULL,
                          cut_point = NULL,
                          select_higher = FALSE,
                          plot = TRUE,
                          type = "var",
                          ...){
  if(is.null(object$object_index)){
    stop("'object' was not computed using the `object_index` argument.")
  }
  coords <- object$results[2:3]
  obj_in <- check_inf(object$object_index)
  obj_in$id <- as.character(obj_in$id)
  if(!is.null(index)){

    if(isFALSE(select_higher)){
      ids <- obj_in$id[which(obj_in[[index]] <= cut_point)]
    } else{
      ids <- obj_in$id[which(obj_in[[index]] >= cut_point)]
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

    between_id <- data.frame(
      n = nrow(obj_in),
      nsel = length(ids),
      prop = length(ids) / nrow(obj_in),
      mean_index_sel = mean(obj_in[[index]][obj_in$id %in% ids], na.rm = TRUE),
      mean_index_nsel = mean(obj_in[[index]][!obj_in$id %in% ids], na.rm = TRUE)
    )
    within_id <- cbind(coords, res)[, c(3, 1, 2, 5, 4, 6, 7)]

  } else{
    ids <- NULL
    between_id <- NULL
    within_id <- NULL
  }

  obj_in_pca <- column_to_rownames(obj_in, "id")
  if (ncol(obj_in_pca) >  1){
    pca_res <- pca(obj_in_pca)
    if (isTRUE(plot)){
      plot(pca_res, type = type, ...)
    }
  } else{
    pca <- NULL
    pca_res <- NULL
  }

  invisible(list(ids = ids,
                 between_id = between_id,
                 within_id = within_id,
                 pca_res = pca_res))
}




names_measures <- function(){
  c("id",
    "x",
    "y",
    "area",
    "area_ch",
    "perimeter",
    "radius_mean",
    "radius_min",
    "radius_max",
    "radius_sd",
    "diam_mean",
    "diam_min",
    "diam_max",
    "major_axis",
    "minor_axis",
    "caliper",
    "length",
    "width",
    "radius_ratio",
    "theta",
    "eccentricity",
    "form_factor",
    "narrow_factor",
    "asp_ratio",
    "rectangularity",
    "pd_ratio",
    "plw_ratio",
    "solidity",
    "convexity",
    "elongation",
    "circularity",
    "circularity_haralick",
    "circularity_norm",
    "coverage")
}

har_names <- function(){
  c( "asm",
     "con",
     "cor",
     "var",
     "idm",
     "sav",
     "sva",
     "sen",
     "ent",
     "dva",
     "den",
     "f12",
     "f13")
}

## helper function to compute the measures based on a mask



features_moment <- function(x){
  mc <- poly_mass(x)
  moms <- t(sapply(x, function(x){
    help_moments(x)
  }))
  res <- cbind(mc, moms)
  colnames(res) <- c("mx", "my", "maj_axis", "min_axis", "eccentricity", "theta")
  invisible(data.frame(res))
}

features_shape <- function(x){
  perimeter <- sapply(x, function(x) {
    sum(help_distpts(x))
  })
  distp <- lapply(x, function(x) {
    help_centdist(x)
  })
  rmean <- mean_list(distp)
  rmin <- min_list(distp)
  rmax <- max_list(distp)
  rsd <- sd_list(distp)
  res <- data.frame(perimeter = perimeter,
                    radius_mean = rmean,
                    radius_min = rmin,
                    radius_max = rmax,
                    radius_sd = rsd)
  invisible(res)
}

## helper function to compute the measures based on a mask
compute_measures <- function(mask,
                             img,
                             haralick =  FALSE,
                             har_nbins = 32,
                             har_scales = 1,
                             har_band = 1){
  ocont <- EBImage::ocontour(mask)
  shape <-
    cbind(features_moment(ocont),
          cbind(area = get_area_mask(mask), features_shape(ocont)))
  valid <- which(shape$mx != "NaN")
  shape <- shape[valid, ]
  coverage <- length(which(mask != 0)) / length(mask)
  ocont <- ocont[valid]
  names(ocont) <- valid
  ch <- conv_hull(ocont)
  area_ch <- help_area(ch)
  caliper = poly_caliper(ocont)
  lw <- help_lw(ocont)
  shape <- transform(shape,
                     id = as.numeric(valid),
                     radius_ratio = radius_max / radius_min,
                     diam_mean = radius_mean * 2,
                     diam_min = radius_min * 2,
                     diam_max = radius_max * 2,
                     length = lw[, 1],
                     width = lw[, 2],
                     coverage = area / length(mask),
                     area_ch =   area_ch,
                     solidity = area / area_ch,
                     caliper = caliper,
                     form_factor = 4 * pi * area / perimeter ^ 2,
                     narrow_factor =  caliper / lw[, 1],
                     asp_ratio = lw[, 1] / lw[, 2],
                     rectangularity = lw[, 1]  * lw[, 2] / area,
                     pd_ratio = perimeter / caliper,
                     plw_ratio = perimeter / (lw[, 1]  + lw[, 2]),
                     convexity = poly_convexity(ocont),
                     elongation = poly_elongation(ocont),
                     circularity = perimeter ^ 2 / area,
                     circularity_haralick = radius_mean / radius_sd,
                     circularity_norm = poly_circularity_norm(ocont))
  shape <- shape[, c("id",
                     "mx",
                     "my",
                     "area",
                     "area_ch",
                     "perimeter",
                     "radius_mean",
                     "radius_min",
                     "radius_max",
                     "radius_sd",
                     "diam_mean",
                     "diam_min",
                     "diam_max",
                     "maj_axis",
                     "min_axis",
                     "caliper",
                     "length",
                     "width",
                     "radius_ratio",
                     "theta",
                     "eccentricity",
                     "form_factor",
                     "narrow_factor",
                     "asp_ratio",
                     "rectangularity",
                     "pd_ratio",
                     "plw_ratio",
                     "solidity",
                     "convexity",
                     "elongation",
                     "circularity",
                     "circularity_haralick",
                     "circularity_norm",
                     "coverage")]
  colnames(shape) <- names_measures()
  if(isTRUE(haralick)){
    hal <- data.frame(
      EBImage::computeFeatures.haralick(mask,
                                        img[,,har_band],
                                        haralick.nbins = har_nbins,
                                        haralick.scales = har_scales)
    )
    shape <- cbind(shape, hal[valid, ])
    colnames(shape) <- c(names_measures(), har_names())
  }
  invisible(list(shape = shape,
              cont = ocont,
              ch = ch))
}

## helper function to compute the measures based on a mask
compute_measures_disease <- function(mask){
  ocont <- EBImage::ocontour(mask)
  shape <-
    cbind(features_moment(ocont),
          cbind(area = get_area_mask(mask), features_shape(ocont)))
  valid <- which(shape$mx != "NaN")
  shape <- shape[valid, ]
  ocont <- ocont[valid]
  names(ocont) <- valid
  lw <- help_lw(ocont)
  shape <- transform(shape,
                     id = as.numeric(valid),
                     radius_ratio = radius_max / radius_min,
                     diam_mean = radius_mean * 2,
                     diam_min = radius_min * 2,
                     diam_max = radius_max * 2,
                     length = lw[, 1],
                     width = lw[, 2],
                     form_factor = 4 * pi * area / perimeter ^ 2)
  shape <- shape[, c("id",
                     "mx",
                     "my",
                     "area",
                     "perimeter",
                     "radius_mean",
                     "radius_min",
                     "radius_max",
                     "radius_sd",
                     "diam_mean",
                     "diam_min",
                     "diam_max",
                     "maj_axis",
                     "min_axis",
                     "length",
                     "width")]
  invisible(list(shape = shape,
              cont = ocont))
}


# Helper functions to apply stats on a list

#' These functions applies common statistics to a list of objects, returning a
#' numeric vector.
#'
#' @param x A data.frame or matrix with numeric values.
#' @param ... Further arguments passed on to the R base function (e.g, mean(),
#'   sd(), etc.)
#'
#' @return A numeric vector.
#' @export
#' @name utils_stats
#' @examples
#' mean_list(list(a = 1:10, b = 2:20))
mean_list <- function(x, ...) {
  if (inherits(x, "list")) {
    sapply(x, mean, ...)
  } else{
    mean(x, ...)
  }
}
#' @export
#' @name utils_stats
sd_list <- function(x, ...) {
  if (inherits(x, "list")) {
    sapply(x, sd, ...)
  } else{
    sd(x, ...)
  }
}
#' @export
#' @name utils_stats
max_list <- function(x, ...) {
  if (inherits(x, "list")) {
    sapply(x, max, ...)
  } else{
    max(x, ...)
  }
}
#' @export
#' @name utils_stats
min_list <- function(x, ...) {
  if (inherits(x, "list")) {
    sapply(x, min, ...)
  } else{
    min(x, ...)
  }
}
