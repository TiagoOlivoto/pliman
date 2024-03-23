#'Analyzes objects in an image
#'
#' A lighter option to [analyze_objects()]
#'
#' @export
#' @name analyze_objects_minimal
#' @inheritParams analyze_objects
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' obj <- analyze_objects(img)
#' obj$statistics
#'
#' }
#'
analyze_objects_minimal <- function(img,
                                    segment_objects = TRUE,
                                    reference = FALSE,
                                    reference_area = NULL,
                                    back_fore_index = "R/(G/B)",
                                    fore_ref_index = "B-R",
                                    reference_larger = FALSE,
                                    reference_smaller = FALSE,
                                    pattern = NULL,
                                    parallel = FALSE,
                                    workers = NULL,
                                    watershed = TRUE,
                                    fill_hull = FALSE,
                                    opening = FALSE,
                                    closing = FALSE,
                                    filter = FALSE,
                                    invert = FALSE,
                                    object_size = "medium",
                                    index = "NB",
                                    r = 1,
                                    g = 2,
                                    b = 3,
                                    re = 4,
                                    nir = 5,
                                    threshold = "Otsu",
                                    tolerance = NULL,
                                    extension = NULL,
                                    lower_noise = 0.10,
                                    lower_size = NULL,
                                    upper_size = NULL,
                                    topn_lower = NULL,
                                    topn_upper = NULL,
                                    lower_eccent = NULL,
                                    upper_eccent = NULL,
                                    lower_circ = NULL,
                                    upper_circ = NULL,
                                    plot = TRUE,
                                    show_original = TRUE,
                                    show_contour = TRUE,
                                    contour_col = "red",
                                    contour_size = 1,
                                    col_foreground = NULL,
                                    col_background = NULL,
                                    marker = FALSE,
                                    marker_col = NULL,
                                    marker_size = NULL,
                                    save_image = FALSE,
                                    prefix = "proc_",
                                    dir_original = NULL,
                                    dir_processed = NULL,
                                    verbose = TRUE){
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)
  if(!object_size %in% c("small", "medium", "large", "elarge")){
    stop("'object_size' must be one of 'small', 'medium', 'large', or 'elarge'")
  }
  if(!missing(img) & !missing(pattern)){
    stop("Only one of `img` or `pattern` arguments can be used.", call. = FALSE)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  help_count <-
    function(img, fill_hull, threshold, opening, closing, filter, tolerance, extension,  plot,
             show_original,  marker, marker_col, marker_size,
             save_image, prefix, dir_original, dir_processed, verbose,
             col_background, col_foreground, lower_noise){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      # when reference is not used
      if(isFALSE(reference)){
        if(isTRUE(segment_objects)){
          img2 <- help_binary(img,
                              index = index,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              invert = invert,
                              fill_hull = fill_hull,
                              threshold = threshold,
                              opening = opening,
                              closing = closing,
                              filter = filter,
                              resize = FALSE)
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img2)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(img2),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(img2)
          }
        } else{
          img2 <- img[,,1]
          img2[img2@.Data == 0 | img2@.Data != 0] <- TRUE
          nmask <- EBImage::bwlabel(img2)
        }

        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
        if(isTRUE(fill_hull)){
          nmask <- EBImage::fillHull(nmask)
        }
        shape <- compute_measures_minimal(mask = nmask)
        object_contour <- shape$cont
        shape <- shape$shape

      } else{
        # when reference is used
        if(is.null(reference_area)){
          stop("A known area must be declared when a template is used.", call. = FALSE)
        }
        if(isFALSE(reference_larger) & isFALSE(reference_smaller)){
          # segment back and fore
          if(!isFALSE(invert)){
            invert1 <- ifelse(length(invert) == 1, invert, invert[1])
          } else{
            invert1 <- FALSE
          }
          img_bf <-
            help_binary(img,
                        threshold = threshold,
                        index = back_fore_index,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        invert = invert1,
                        fill_hull = fill_hull)
          img3 <- img
          img3@.Data[,,1][which(img_bf != 1)] <- 2
          img3@.Data[,,2][which(img_bf != 1)] <- 2
          img3@.Data[,,3][which(img_bf != 1)] <- 2
          ID <-  which(img_bf == 1) # IDs for foreground
          ID2 <- which(img_bf == 0) # IDs for background
          # segment fore and ref
          if(!isFALSE(invert)){
            invert2 <- ifelse(length(invert) == 1, invert, invert[2])
          } else{
            invert2 <- FALSE
          }
          img4 <-
            help_binary(img3,
                        threshold = threshold,
                        index = fore_ref_index,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        invert = invert2)
          mask <- img_bf
          pix_ref <- which(img4 != 1)
          img@.Data[,,1][pix_ref] <- 1
          img@.Data[,,2][pix_ref] <- 0
          img@.Data[,,3][pix_ref] <- 0
          npix_ref <- length(pix_ref)
          mask[pix_ref] <- 0
          if(is.numeric(filter) & filter > 1){
            mask <- EBImage::medianFilter(mask, size = filter)
          }
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(mask),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(mask)
          }

          shape <- compute_measures_minimal(mask = nmask)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape

          # correct measures based on the area of the reference object
          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side^2
          shape[5:8] <- apply(shape[5:8], 2, function(x){
            x * px_side
          })
        } else{
          # correct the measures based on larger or smaller objects
          mask <-
            help_binary(img,
                        threshold = threshold,
                        index = index,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        opening = opening,
                        closing = closing,
                        filter = filter,
                        invert = invert,
                        fill_hull = fill_hull)
          ID <-  which(mask == 1) # IDs for foreground
          ID2 <- which(mask == 0) # IDs for background
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(mask)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(mask),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(mask)
          }

          shape <- compute_measures_minimal(mask = nmask)

          object_contour <- shape$cont
          shape <- shape$shape

          if(isTRUE(reference_larger)){
            id_ref <- which.max(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
          } else{
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
            id_ref <- which.min(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
          }

          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side ^ 2
          shape[5:8] <- apply(shape[5:8], 2, function(x){
            x * px_side
          })
        }
      }


      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        stop("Only one of 'lower_*' or 'topn_*' can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$area > lower_size, ],
             shape <- shape[shape$area > mean(shape$area) * lower_noise, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$area, decreasing = TRUE),][1:topn_upper,]
      }
      if(!is.null(lower_eccent)){
        shape <- shape[shape$eccentricity > lower_eccent, ]
      }
      if(!is.null(upper_eccent)){
        shape <- shape[shape$eccentricity < upper_eccent, ]
      }
      if(!is.null(lower_circ)){
        shape <- shape[shape$circularity > lower_circ, ]
      }
      if(!is.null(upper_circ)){
        shape <- shape[shape$circularity < upper_circ, ]
      }
      object_contour <- object_contour[as.character(shape$id)]



      stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area"),
                          value = c(length(shape$area),
                                    min(shape$area),
                                    mean(shape$area),
                                    max(shape$area)))
      results <- list(results = shape,
                      statistics = stats)
      class(results) <- "anal_obj_minimal"
      if(plot == TRUE | save_image == TRUE){
        backg <- !is.null(col_background)
        # color for background
        if (is.null(col_background)){
          col_background <- col2rgb("white") / 255
        } else{
          ifelse(is.character(col_background),
                 col_background <- col2rgb(col_background) / 255,
                 col_background <- col_background / 255)
        }
        # color for lesions
        if (is.null(col_foreground)){
          col_foreground <- col2rgb("black") / 255
        } else{
          ifelse(is.character(col_foreground),
                 col_foreground <- col2rgb(col_foreground) / 255,
                 col_foreground <- col_foreground / 255)
        }

        if(show_original == TRUE){
          im2 <- img[,,1:3]
          EBImage::colorMode(im2) <- "Color"
          if(backg){
            im3 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- col_background[3]
          }
        }
        show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
        marker <- ifelse(is.null(marker), "id", marker)
        if(!isFALSE(show_mark) & marker != "point" & !marker %in% colnames(shape)){
          warning("Accepted 'marker' are: {", paste(colnames(shape), collapse = ", "),
                  "}. Drawing the object id.", call. = FALSE)
          marker <- "id"
        }
        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
        # correct the contour
        object_contour <- lapply(object_contour, function(x){
          x + 1
        })

        if(plot == TRUE){
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2] + 1,
                   shape[, 3] + 1,
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour)  & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2] + 1,
                     shape[, 3] + 1,
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
        }

        if(save_image == TRUE){
          if(dir.exists(diretorio_processada) == FALSE){
            dir.create(diretorio_processada, recursive = TRUE)
          }
          png(paste0(diretorio_processada, "/",
                     prefix,
                     name_ori, ".",
                     extens_ori),
              width = dim(im2@.Data)[1],
              height = dim(im2@.Data)[2])
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2] + 1,
                   shape[, 3] + 1,
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2] + 1,
                     shape[, 3] + 1,
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          dev.off()
        }
      }
      invisible(results)
    }

  if(missing(pattern)){
    help_count(img, fill_hull, threshold, opening, closing, filter, tolerance, extension,  plot,
               show_original,  marker, marker_col, marker_size,
               save_image, prefix, dir_original, dir_processed, verbose,
               col_background, col_foreground, lower_noise)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("Pattern '", pattern, "' not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      init_time <- Sys.time()
      nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores()*.3), workers)
      future::plan(future::multisession, workers = nworkers)
      on.exit(future::plan(future::sequential))
      `%dofut%` <- doFuture::`%dofuture%`

      if(verbose == TRUE){
        message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
      }

      results <-
        foreach::foreach(i = seq_along(names_plant)) %dofut%{
          help_count(names_plant[i],
                     fill_hull, threshold, opening, closing, filter, tolerance, extension,  plot,
                     show_original,  marker, marker_col, marker_size,
                     save_image, prefix, dir_original, dir_processed, verbose,
                     col_background, col_foreground, lower_noise)
        }

    } else{
      init_time <- Sys.time()
      pb <- progress(max = length(plants), style = 4)
      foo <- function(plants, ...){
        if(verbose == TRUE){
          run_progress(pb, ...)
        }
        help_count(img  = plants,
                   fill_hull, threshold, filter, tolerance, extension,  plot,
                   show_original,  marker, marker_col, marker_size,
                   save_image, prefix, dir_original, dir_processed, verbose,
                   col_background, col_foreground, lower_noise)
      }
      results <-
        lapply(seq_along(names_plant), function(i){
          foo(names_plant[i],
              actual = i,
              text = paste("Processing image", names_plant[i]))
        })
    }

    ## bind the results
    names(results) <- names_plant

    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["statistics"]],
                          id =  names(results[i]))[,c(3, 1, 2)]
              })
      )

    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["results"]],
                          img =  names(results[i]))
              })
      )

    if("img" %in% colnames(results)){
      results <- results[, c(ncol(results), 1:ncol(results) - 1)]
    }
    summ <- stats[stats$stat == "n", c(1, 3)]
    if(verbose == TRUE){
      names(summ) <- c("Image", "Objects")
      cat("--------------------------------------------\n")
      print(summ, row.names = FALSE)
      cat("--------------------------------------------\n")
      message("Done!")
      message("Elapsed time: ", sec_to_hms(as.numeric(difftime(Sys.time(),  init_time, units = "secs"))))

    }

    invisible(
      structure(
        list(statistics = stats,
             count = summ,
             results = results),
        class = "anal_obj_ls_minimal"
      )
    )
  }
}


#' @name analyze_objects_minimal
#' @inheritParams plot.analyze_objects
#' @method plot anal_obj_minimal
#' @export
#'
#' @examples
#' \donttest{
#' library(pliman)
#'
#' img <- image_pliman("soy_green.jpg")
#' # Segment the foreground (grains) using the normalized blue index (NB, default)
#' # Shows the average value of the blue index in each object
#'
#' rgb <- analyze_objects_minimal(img)
#' # density of area
#' plot(rgb)
#'
#' # histogram of area
#' plot(rgb, type = "histogram") # or 'hist'
#' }
plot.anal_obj_minimal <- function(x,
                                  which = "measure",
                                  measure = "area",
                                  type = c("density", "histogram"),
                                  ...){
  if(!which %in% c("measure", "index")){
    stop("'which' must be one of 'measure' or 'index'", call. = FALSE)
  }
  nam <- colnames(x$results)
  if(!measure %in% nam){
    stop("Measure '", measure, "' not available in 'x'. Try one of the '",
         paste0(nam, collapse = ", "), call. = FALSE)
  }
  temp <- x$results[[measure]]
  types <- c("density", "histogram")
  matches <- grepl(type[1], types)
  type <- types[matches]
  if(type == "histogram"){
    hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
  } else{
    density_data <- density(temp)  # Calculate the density for the column
    plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
    points(x = temp, y = rep(0, length(temp)), col = "red")
  }
}
#' @name analyze_objects_minimal
#' @export
plot.anal_obj_ls_minimal <- function(x,
                                     which = "measure",
                                     measure = "area",
                                     type = c("density", "histogram"),
                                     ...){
  if(!which %in% c("measure", "index")){
    stop("'which' must be one of 'measure' or 'index'", call. = FALSE)
  }
  nam <- colnames(x$results)
  if(!measure %in% nam){
    stop("Measure '", measure, "' not available in 'x'. Try one of the '",
         paste0(nam, collapse = ", "), call. = FALSE)
  }
  temp <- x$results[[measure]]
  types <- c("density", "histogram")
  matches <- grepl(type[1], types)
  type <- types[matches]
  if(type == "histogram"){
    hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
  } else{
    density_data <- density(temp)  # Calculate the density for the column
    plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
    points(x = temp, y = rep(0, length(temp)), col = "red")
  }
}


