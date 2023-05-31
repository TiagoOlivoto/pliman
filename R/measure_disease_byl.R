
#' Performs plant disease measurements by leaf
#'
#' Computes the percentage of symptomatic leaf area using color palettes or RGB
#' indexes `by` each `l`eaf of an image. This allows, for example, processing
#' replicates of the same treatment  and obtaining the results for each
#' replication with a single image. To do that, leaf samples are first splitten
#' with [object_split()] and then, [measure_disease()] is applied to the list of
#' leaves.
#'
#' @inheritParams object_split
#' @inheritParams measure_disease
#' @param dir_original The directory containing the original and processed images.
#'   Defaults to `NULL`. In this case, the function will search for the image `img` in the
#'   current working directory.
#' @param index A character value specifying the target mode for
#'   conversion to binary to segment the leaves from background. Defaults to "B"
#'   (blue). See [image_index()] for more details. Personalized indexes can be
#'   informed as, e.g., `index = "R*G/B`.
#' @param lower_size To prevent dust from affecting object segmentation, objects
#'   with lesser than `10%` of the mean of all objects are removed. . One can
#'   set a known area or use `lower_limit = 0` to select all objects (not
#'   advised).
#' @param ... Additional arguments passed on to [measure_disease()].
#' @return
#' * A list with the following objects:
#'   - `severity` A data frame with the percentage of healthy and symptomatic
#'  areas for each leaf in the image(s).
#'   - `shape`,`statistics` If `show_features = TRUE` is used, returns the shape
#'  (area, perimeter, etc.) for each lesion and a summary statistic of the
#'  results.
#' @export
#'
#' @examples
#'library(pliman)
#'img <- image_pliman("mult_leaves.jpg", plot = TRUE)
#'sev <-
#'  measure_disease_byl(img = img,
#'                      index_lb = "B",
#'                      index_dh = "NGRDI",
#'                      workers = 2)
#' sev$severity
#'
#'
measure_disease_byl <- function(img,
                                index = "B",
                                index_lb = "B",
                                index_dh = "NGRDI",
                                lower_size = NULL,
                                watershed = TRUE,
                                invert = FALSE,
                                fill_hull = FALSE,
                                filter = 3,
                                threshold = "Otsu",
                                extension = NULL,
                                tolerance = NULL,
                                object_size = "large",
                                dir_original = NULL,
                                img_healthy = NULL,
                                img_symptoms = NULL,
                                pattern = NULL,
                                parallel = FALSE,
                                workers = NULL,
                                plot = TRUE,
                                show_features = FALSE,
                                save_image = FALSE,
                                verbose = TRUE,
                                ...){

  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if (is.character(img_healthy)){
    all_files <- sapply(list.files(getwd()), file_name)
    imag <- list.files(getwd(), pattern = img_healthy)
    check_names_dir(img_healthy, all_files, getwd())
    name <- file_name(imag)
    extens <- file_extension(imag)
    img_healthy <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
  }
  if (is.character(img_symptoms)){
    all_files <- sapply(list.files(getwd()), file_name)
    imag <- list.files(getwd(), pattern = img_symptoms)
    check_names_dir(img_symptoms, all_files, getwd())
    name <- file_name(imag)
    extens <- file_extension(imag)
    img_symptoms <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
  }
  back <- EBImage::Image(rep(1, 100*300),dim=c(100,300,3), colormode = 'Color')

  help_byl <- function(img,
                       img_healthy,
                       img_symptoms,
                       back,
                       index_dh,
                       index_lb){
    if(is.character(img)){
      all_files <- sapply(list.files(diretorio_original), file_name)
      check_names_dir(img, all_files, diretorio_original)
      imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
      name_ori <- file_name(imag)
      extens_ori <- file_extension(imag)
      img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
    } else{
      name_ori <- match.call()[[2]]
      extens_ori <- "jpg"
    }

    splits <- object_split(img,
                           index = index,
                           watershed = watershed,
                           invert = invert,
                           fill_hull = fill_hull,
                           filter = filter,
                           threshold = threshold,
                           extension = extension,
                           tolerance = tolerance,
                           object_size = object_size,
                           plot = FALSE,
                           verbose = FALSE)
    results <- list()
    tmp_dir <- tempdir()
    if(isTRUE(plot)){
      save_image <- FALSE
      clear_td()
    }
    if(is.null(img_healthy)){
      results <-
        lapply(seq_along(splits),
               function(i){
                 measure_disease(splits[[i]],
                                 index_dh = index_dh,
                                 index_lb = index_lb,
                                 plot = FALSE,
                                 save_image = TRUE,
                                 show_features = show_features,
                                 dir_processed = tmp_dir,
                                 prefix = paste0(name_ori, "_", i),
                                 name = "",
                                 filter = FALSE,
                                 ...)
               })

    } else{
      results <-
        lapply(seq_along(splits),
               function(i){
                 measure_disease(splits[[i]],
                                 img_healthy = img_healthy,
                                 img_symptoms = img_symptoms,
                                 img_background = back,
                                 plot = FALSE,
                                 save_image = TRUE,
                                 show_features = show_features,
                                 dir_processed = tmp_dir,
                                 prefix = paste0(name_ori, "_", i),
                                 name = "",
                                 filter = FALSE,
                                 ...)
               })
    }
    names(results) <- paste0(name_ori, "-", 1:length(results))
    if(isTRUE(plot)){
      imgs <- image_import(pattern = as.character(name_ori), path = tmp_dir)
      image_combine(imgs)
      clear_td()
    }
    severity <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["severity"]],
                          img =  names(results[i]))[, c(3, 1:2)]
              })
      ) |>
      separate_col(img, into = c("img", "leaf"), sep = "-")
    if(isTRUE(show_features)){
      stats <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  if(is.null(results[[i]][["statistics"]])){
                    res <- data.frame(stat = "NA", value = 0)
                  } else{
                    res <- results[[i]][["statistics"]]
                  }
                  transform(res, img =  names(results[i]))[,c(3, 1, 2)]
                })
        ) |>
        separate_col(img, into = c("img", "leaf"), sep = "-")
      shape <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  if(is.null(results[[i]][["shape"]])){
                    names <- names_measures()
                    res <- data.frame(matrix(nrow = 1, ncol = length(names)))
                    res[1, ] <- 0
                    colnames(res) <- names
                  } else{
                    res <- results[[i]][["shape"]]
                  }
                  transform(res, img =  names(results[i]))
                })
        ) |>
        separate_col(img, into = c("img", "leaf"), sep = "-")

    } else{
      shape <- NULL
      stats <- NULL
    }
    return(
      structure(
        list(severity = severity,
             stats = stats,
             shape = shape),
        class = "plm_disease_byl")
    )
  }

  if(missing(pattern)){
    results <- help_byl(img, img_healthy, img_symptoms, back, index_dh, index_lb)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("'", pattern, "' pattern not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_byl", "img_healthy", "img_symptoms", "back"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_byl(img  = x, img_healthy, img_symptoms, back, index_dh, index_lb)
                  })
    } else{
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        if(verbose == TRUE){
          run_progress(pb, actual = i,
                       text = paste("Processing image", names_plant[i]))
        }
        results[[i]] <- help_byl(img  = names_plant[i], img_healthy, img_symptoms, back, index_dh, index_lb)
      }
    }
    names(results) <- names_plant
    if(isTRUE(show_features)){
      stats <- do.call(rbind, lapply(results, function(x){x[["stats"]]}))
      rownames(stats) <- NULL
      shape <- do.call(rbind, lapply(results, function(x){x[["shape"]]}))
      rownames(shape) <- NULL
    } else{
      shape <- NULL
      stats <- NULL
    }
    severity <- do.call(rbind, lapply(results, function(x){x[["severity"]]}))
    rownames(severity) <- NULL

    results <- list(severity = severity,
                    shape = shape,
                    stats = stats,
                    parms = list(
                      pattern = pattern,
                      img_healthy = img_healthy,
                      img_symptoms = img_symptoms,
                      dir_original = diretorio_original,
                      save_image = save_image))
  }
  return(structure(
    results, class = "plm_disease_byl"
  ))
}




