#' Construct a shape file from an image
#'
#' Creates a list of object coordinates given the desired number of nrow and
#' columns. It starts by selecting 4 points at the corners of objects of
#' interest in the plot space. Then, given `nrow` and `ncol`, a grid is drawn
#' and the objects' coordinates are returned.
#'
#' @param img An object of class `Image`
#' @param nrow The number of desired rows in the grid. Defaults to `1`.
#' @param ncol The number of desired columns in the grid. Defaults to `1`.
#' @param buffer_x,buffer_y Buffering factor for the width and height,
#'   respectively, of each individual shape's side. A value between 0 and 0.5 where 0
#'   means no buffering and 0.5 means complete buffering (default: 0). A value of
#'   0.25 will buffer the shape by 25% on each side.
#' @param interactive If `FALSE` (default) the grid is created automatically
#'   based on the image dimension and number of rows/columns. If `interactive =
#'   TRUE`, users must draw points at the diagonal of the desired bounding box
#'   that will contain the grid.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are `"base"` and
#'   `"mapview"`. If set to "base", the base R graphics system is used for
#'   interactive plotting. If set to "mapview", the mapview package is used. To
#'   set this argument globally for all functions in the package, you can use
#'   the [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param col_line,col_text The color of the line/text in the grid. Defaults to
#'   `"red"`.
#' @param size_line,size_text The size of the line/text in the grid. Defaults to
#'   `2.5`.
#' @param plot Plots the grid on the image? Defaults to `TRUE`.
#' @return A list with `row * col` objects containing the plot coordinates.
#' @export
#'
#' @examples
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg")
#' shape <- image_shp(flax, nrow = 3, ncol = 5)
#'
image_shp <- function(img,
                      nrow = 1,
                      ncol = 1,
                      buffer_x = 0,
                      buffer_y = 0,
                      interactive = FALSE,
                      viewer = get_pliman_viewer(),
                      col_line = "red",
                      size_line = 2,
                      col_text = "red",
                      size_text = 1,
                      plot = TRUE){
  viewopt <- c("base", "mapview")
  viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]
  if(isTRUE(interactive)){
    if(viewopt == "base"){
      message("Select 2 points drawing the diagonal that includes the objects of interest.")
      plot(img)
      cord <- locator(type = "p", n = 2, col = "red", pch = 19)
      c1 <- data.frame(do.call(rbind, cord)) |> t()
      c1 <- c(min(c1[,1]), max(c1[,1]), min(c1[,2]), max(c1[,2]))
    } else{
      coords <- mv_rectangle(img)
      c1 <- c(min(coords[,1]), max(coords[,1]), min(coords[,2]), max(coords[,2]))
    }
  } else{
    imgd <- dim(img)
    c1 <- c(0, imgd[1], 0, imgd[2])
  }
  bbox <-
    data.frame(x = c(c1[1], c1[2], c1[2], c1[1], c1[1]),
               y = c(c1[3], c1[3], c1[4], c1[4], c1[3]))
  shps <- help_shp(nrow, ncol, c1, buffer_x = buffer_x, buffer_y = buffer_y)
  shps <- data.frame(plot = paste0(rep(1:(ncol * nrow), each = 5)), shps)
  colnames(shps) <- c("plot", "x", "y")
  coords <- split(shps, shps$plot)
  names(coords) <- paste0("plot_", names(coords))
  coords <- coords[paste0("plot_", 1:length(coords))]
  if(isTRUE(plot)){
    plot(img)
    plot_shp(coords,
             col_line = col_line,
             size_line = size_line,
             col_text = col_text,
             size_text = size_text)
  }
  lst <- list(shapefiles = coords,
              bbox = bbox,
              nrow = nrow,
              ncol = ncol)
  invisible(structure(lst, class = "image_shp"))
}



#' S3 method `plot` for `image_shp` objects
#'
#' Draws the bounding boxes for each object computed with [image_shp()].
#'
#' @inheritParams  image_shp
#' @param x An object computed with [image_shp()].
#' @param img The image that was used to compute the shapefile (optional)
#' @param ... Currently not used.
#' @return A `NULL` object
#' @export
#' @importFrom grDevices dev.list
#'
#' @examples
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg")
#' shape <- image_shp(flax, nrow = 3, ncol = 5)
#'
#' # grid on the existing image
#' plot(flax)
#' plot(shape)
plot.image_shp <- function(x,
                           img = NULL,
                           col_line = "black",
                           size_line = 2,
                           col_text = "black",
                           size_text = 0.75,
                           ...){
  shapefiles <- x$shapefiles
  bbox <- x$bbox
  if(is.null(img)){
    if(is.null(dev.list())){
      plot(1,
           axes = FALSE,
           type = "n",
           asp = 1,
           xlab="",
           ylab="",
           xlim=c(min(bbox[,1]), max(bbox[,1])),
           ylim=c(min(bbox[,2]), max(bbox[,2])))
    }
  } else{
    plot(img)
  }
  con <- 0
  for (i in seq_along(shapefiles)) {
    con <- con + 1
    tmp <- shapefiles[[i]]
    lines(tmp[, -1], col = col_line, type = "l", lwd = size_line)
    text(min(tmp$x[-1]), min(tmp$y[-1]),
         label = con,
         col = col_text,
         cex = size_text,
         adj = c(-0.2, 1.2))
  }
}

#' Splits image objects based on a shape file
#'
#'  Here, [image_shp()] is used to create a shape file based on the desired
#'  number of rows and columns. Then, using the object coordinates, a list of
#'  `Image` objects is created.
#' @inheritParams  image_shp
#' @param only_shp If `TRUE` returns only the shapefiles with the coordinates
#'   for each image. If `FALSE` (default) returns the splitted image according
#'   to `nrow` and `ncol` arguments.
#' @param ... Other arguments passed on to [image_shp()]
#' @return A list of `Image` objects
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' objects <- object_split_shp(flax, nrow = 3, ncol = 5)
#' image_combine(objects$imgs)
#' }
object_split_shp <- function(img,
                             nrow = 1,
                             ncol = 1,
                             buffer_x = 0,
                             buffer_y = 0,
                             interactive = FALSE,
                             viewer = get_pliman_viewer(),
                             only_shp = FALSE,
                             ...){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  shps <- image_shp(img, nrow, ncol,
                    buffer_x = buffer_x,
                    buffer_y = buffer_y,
                    interactive = interactive,
                    plot = FALSE,
                    viewer = vieweropt, ...)
  shapefile <- shps$shapefiles
  if(!isTRUE(only_shp)){
    imgs <- list()
    get_borders <- function(x){
      min_x <- min(x[,1])
      max_x <- max(x[,1])
      min_y <- min(x[,2])
      max_y <- max(x[,2])
      invisible(list(min_x, max_x, min_y, max_y))
    }
    for (i in 1:length(shapefile)) {
      tmp <- shapefile[[i]][-1]
      borders <- get_borders(tmp)
      imgs[[paste0("shp", i)]] <-
        image_crop(img,
                   width = borders[[1]]:borders[[2]],
                   height = borders[[3]]:borders[[4]])
    }
  } else{
    imgs <- img
  }
  invisible(list(imgs = imgs,
              shapefile = shps))
}


#' Export multiple objects from an image to multiple images
#'
#' Givin an image with multiple objects, `object_export_shp()` will split the
#' objects into a list of objects using [object_split_shp()] and then export
#' them to multiple images into the current working directory (or a subfolder).
#' Batch processing is performed by declaring a file name pattern that matches
#' the images within the working directory.
#'
#' @inheritParams object_split_shp
#' @inheritParams utils_image
#' @inheritParams analyze_objects
#'
#' @param pattern A pattern of file name used to identify images to be
#'   processed. For example, if `pattern = "im"` all images in the current
#'   working directory that the name matches the pattern (e.g., img1.-,
#'   image1.-, im2.-) will be imported and processed. Providing any number as
#'   pattern (e.g., `pattern = "1"`) will select images that are named as 1.-,
#'   2.-, and so on. An error will be returned if the pattern matches any file
#'   that is not supported (e.g., img1.pdf).
#'@param dir_original The directory containing the original images. Defaults to
#'  `NULL`. It can be either a full path, e.g., `"C:/Desktop/imgs"`, or a
#'  subfolder within the current working directory, e.g., `"/imgs"`.
#' @param dir_processed Optional character string indicating a subfolder within the
#'   current working directory to save the image(s). If the folder doesn't
#'   exist, it will be created.
#' @param format The format of image to be exported.
#' @param squarize Squarizes the image before the exportation? If `TRUE`,
#'   [image_square()] will be called internally.
#' @return A `NULL` object.
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' object_export_shp(flax)
#'
#' }
object_export_shp <- function(img,
                              pattern = NULL,
                              dir_original = NULL,
                              dir_processed = NULL,
                              format = ".jpg",
                              subfolder = NULL,
                              squarize = FALSE,
                              nrow = 1,
                              ncol = 1,
                              buffer_x = 0,
                              buffer_y = 0,
                              interactive = FALSE,
                              parallel = FALSE,
                              verbose = TRUE,
                              viewer = get_pliman_viewer()){
  if(is.null(pattern)){
    list_objects <- object_split_shp(img,
                                     nrow = nrow,
                                     ncol = ncol,
                                     buffer_x = buffer_x,
                                     buffer_y = buffer_y,
                                     interactive = interactive,
                                     viewer = viewer)[["imgs"]]

    a <- lapply(seq_along(list_objects), function(i){
      tmp <- list_objects[[i]]
      if(isTRUE(squarize)){
        tmp <- image_square(tmp,
                            plot = FALSE,
                            sample_left = 10,
                            sample_top = 10,
                            sample_right = 10,
                            sample_bottom = 10)
      }
      image_export(tmp,
                   name = paste0(names(list_objects[i]), format),
                   subfolder = dir_processed)
    })
  } else{
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


    if(isTRUE(parallel)){

      init_time <- Sys.time()
      nworkers <- trunc(detectCores()*.3)
      cl <- parallel::makePSOCKcluster(nworkers)
      doParallel::registerDoParallel(cl)
      on.exit(stopCluster(cl))

      if(verbose == TRUE){
        message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
      }
      ## declare alias for dopar command
      `%dopar%` <- foreach::`%dopar%`

      results <-
        foreach::foreach(i = seq_along(plants), .packages = c("pliman")) %dopar%{

          tmpimg <- image_import(plants[[i]], path = diretorio_original)

          list_objects <- object_split_shp(tmpimg,
                                           nrow = nrow,
                                           ncol = ncol,
                                           buffer_x = buffer_x,
                                           buffer_y = buffer_y,
                                           interactive = interactive,
                                           viewer = viewer)[["imgs"]]

          a <- lapply(seq_along(list_objects), function(j){
            tmp <- list_objects[[j]]
            if(isTRUE(squarize)){
              tmp <- image_square(tmp,
                                  plot = FALSE,
                                  sample_left = 10,
                                  sample_top = 10,
                                  sample_right = 10,
                                  sample_bottom = 10)
            }
            image_export(tmp,
                         name = paste0(file_name(plants[[i]]), "_", names(list_objects[j]), format),
                         subfolder = diretorio_processada)
          }
          )
        }

      message("Done!")
      message("Elapsed time: ", sec_to_hms(as.numeric(difftime(Sys.time(),  init_time, units = "secs"))))

    } else{

      for(i in seq_along(plants)){
        tmpimg <- image_import(plants[[i]], path = diretorio_original)

        list_objects <- object_split_shp(tmpimg,
                                         nrow = nrow,
                                         ncol = ncol,
                                         buffer_x = buffer_x,
                                         buffer_y = buffer_y,
                                         interactive = interactive,
                                         viewer = viewer)[["imgs"]]

        a <- lapply(seq_along(list_objects), function(j){
          tmp <- list_objects[[j]]
          if(isTRUE(squarize)){
            tmp <- image_square(tmp,
                                plot = FALSE,
                                sample_left = 10,
                                sample_top = 10,
                                sample_right = 10,
                                sample_bottom = 10)
          }
          image_export(tmp,
                       name = paste0(file_name(plants[[i]]), "_", names(list_objects[j]), format),
                       subfolder = diretorio_processada)
        })
      }
    }


  }

}


#' Aligns an `Image` object by hand
#'
#' [image_align()] rotate an image given a line of desired aligment along the y
#' axis that corresponds to the alignment of the objects (e.g., field plots). By
#' default, the aligment will be to the vertical, which means that if the drawed
#' line have an angle < 90º parallel to the x axis, the rotation angle wil be
#' negative (anticlocwise rotation).
#'
#' @details
#' The `image_align` function aligns an image along the vertical or horizontal
#' axis based on user-selected points. The alignment can be performed in either
#' the base plotting system or using the mapview package for interactive
#' visualization. If the viewer option is set to "base", the function prompts
#' the user to select two points on the image to define the alignment line. If
#' the viewer option is set to "mapview", the function opens an interactive map
#' where the user can draw a polyline to define the alignment line. The
#' alignment angle is calculated based on the selected points, and the image is
#' rotated accordingly using the `image_rotate` function. The function returns
#' the aligned image object.
#'
#' @param img An `Image` object
#' @param align The desired alignment. Either `"vertical"` (default) or
#'   `"horizontal"`.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param plot Plots the aligned image? Defaults to `TRUE`.
#'
#' @return The `img` aligned
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' aligned <- image_align(flax)
#' }
image_align <- function(img,
                        align = c("vertical", "horizontal"),
                        viewer = get_pliman_viewer(),
                        plot = TRUE){
  alignopt <- c("vertical", "horizontal")
  alignopt <- alignopt[pmatch(align[1], alignopt)]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(vieweropt == "base"){
    message("Select 2 points drawing a line of desired aligment along the y axis.")
    plot(EBImage::Image(img[,,1:3], colormode = "Color"))
    cord <- locator(type = "p", n = 2, col = "red", pch = 19)
    c1 <- data.frame(do.call(rbind, cord)) |> t()
    lines(c1, col = "red", lty = 2, lwd = 2)
    a <- abs(cord$x[1] - cord$x[2])
    b <- abs(cord$y[1] - cord$y[2])
  } else{
    mv <- mv_two_points(img,
                        title = "Use the 'Draw Polyline' tool to Select 2 points drawing a line of desired aligment")
    a <- abs(mv[[1]] - mv[[3]])
    b <- abs(mv[[2]] - mv[[4]])
    c1 <- data.frame(x = c(mv[[1]], mv[[3]]),
                     y = c(mv[[2]], mv[[4]]))
  }
  angle <- (atan(b / a) * 180) / pi
  if(svd(var(c1))$u[2] >= 0){
    anglev <- angle - 90
    angleh <- angle
    if(alignopt[[1]] == "vertical"){
      img2 <- image_rotate(img, angle = anglev, plot = plot)
    } else{
      img2 <- image_rotate(img, angle = angleh, plot = plot)
    }
    message(paste("Angle to align in the vertical: ", round(anglev, 3)))
    message(paste("Angle to align in the horizontal: ", round(angleh, 3)))
  } else{
    anglev <- 90 - angle
    angleh <- angle * -1
    if(vieweropt == "vertical"){
      img2 <- image_rotate(img, angle = anglev, plot = plot)
    } else{
      img2 <- image_rotate(img, angle = angleh, plot = plot)
    }
    message(paste("Angle to align in the vertical: ", round(anglev, 3)))
    message(paste("Angle to align in the horizontal: ", round(angleh, 3)))
  }
  invisible(img2)
}

#' Analyzes objects using shapefiles
#'
#' @details The `analyze_objects_shp` function performs object analysis on an
#' image and generates shapefiles representing the analyzed objects. The
#' function first prepares the image for analysis using the [image_prepare()]
#' function if the `prepare` argument is set to `TRUE`. If a shapefile object is
#' provided, the number of rows and columns for splitting the image is obtained
#' from the shapefile. Otherwise, the image is split into multiple sub-images
#' based on the specified number of rows and columns using the
#' [object_split_shp()] function. The objects in each sub-image are analyzed
#' using the [analyze_objects()] function, and the results are stored in a list.
#' If parallel processing is enabled, the analysis is performed in parallel
#' using multiple workers.
#'
#' The output object provides access to various components of the analysis
#' results, such as the analyzed object coordinates and properties.
#' Additionally, the shapefiles representing the analyzed objects are included
#' in the output object for further analysis or visualization.
#'
#'
#'
#' @inheritParams analyze_objects
#' @inheritParams  image_shp
#'
#' @param img An `Image` object
#' @param nrow,ncol The number of rows and columns to generate the shapefile
#'   when `shapefile` is not declared. Defaults to `1`.
#' @param prepare Logical value indicating whether to prepare the image for
#'   analysis using [image_prepare()] function. Defaults to `FALSE`. Set to
#'   `TRUE` to interactively align and crop the image before processing.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param shapefile (Optional) An object created with [image_shp()]. If `NULL`
#'   (default), both `nrow` and `ncol` must be declared.
#' @param interactive If `FALSE` (default) the grid is created automatically
#'   based on the image dimension and number of nrow/columns. If `interactive =
#'   TRUE`, users must draw points at the diagonal of the desired bounding box
#'   that will contain the grid.
#' @param plot Plots the processed images? Defaults to `FALSE`.
#' @param object_size Argument to control control the watershed segmentation.
#'   See [analyze_objects()] for more details.
#' @param ... Aditional arguments passed on to [analyze_objects].
#'
#' @return An object of class `anal_obj`. See more details in the `Value`
#'   section of [analyze_objects()].
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#'
#' # Computes the DGCI index for each flax leaf
#' flax <- image_pliman("flax_leaves.jpg", plot =TRUE)
#' res <-
#'    analyze_objects_shp(flax,
#'                        nrow = 3,
#'                        ncol = 5,
#'                        plot = FALSE,
#'                        object_index = "DGCI")
#' plot(flax)
#' plot(res$shapefiles)
#' plot_measures(res, measure = "DGCI")
#' }
analyze_objects_shp <- function(img,
                                nrow = 1,
                                ncol = 1,
                                buffer_x = 0,
                                buffer_y = 0,
                                prepare = FALSE,
                                segment_objects = TRUE,
                                viewer = get_pliman_viewer(),
                                index = "R",
                                r = 1,
                                g = 2,
                                b = 3,
                                re = 4,
                                nir = 5,
                                shapefile = NULL,
                                interactive = FALSE,
                                plot = FALSE,
                                parallel = FALSE,
                                workers = NULL,
                                watershed = TRUE,
                                filter = FALSE,
                                object_size = "medium",
                                efourier = FALSE,
                                object_index = NULL,
                                veins = FALSE,
                                width_at = FALSE,
                                verbose = TRUE,
                                invert = FALSE,
                                ...){
  if(isTRUE(prepare)){
    img <- image_prepare(img, viewer = viewer)
  } else{
    img <- img
  }
  mask <- analyze_objects(img,
                          index = index,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          invert = invert,
                          segment_objects = segment_objects,
                          plot = FALSE,
                          return_mask = TRUE,
                          watershed = watershed,
                          filter = filter,
                          object_size = object_size,
                          object_index = object_index)$mask
  object_index_used <- object_index
  if(is.null(shapefile)){
    tmp <- object_split_shp(img, nrow, ncol,
                            buffer_x = buffer_x,
                            buffer_y = buffer_y,
                            interactive = interactive,
                            only_shp = FALSE)
    imgs <- tmp$imgs
    shapes <- tmp$shapefile$shapefiles
  } else{
    nrow <- shapefile$nrow
    ncol <- shapefile$ncol
    tmp <- object_split_shp(img, nrow, ncol,
                            buffer_x = buffer_x,
                            buffer_y = buffer_y,
                            interactive = FALSE, only_shp = FALSE)
    imgs <- tmp$imgs
    shapes <- tmp$shapefile$shapefiles
  }

  if(isTRUE(plot)){
    op <- par(mfrow = c(nrow, ncol))
    on.exit(par(op))
  }

  if(parallel == TRUE){
    workers <- ifelse(is.null(workers), ceiling(detectCores() * 0.3), workers)
    cl <- parallel::makePSOCKcluster(workers)
    doParallel::registerDoParallel(cl)
    on.exit(stopCluster(cl))

    ## declare alias for dopar command
    `%dopar%` <- foreach::`%dopar%`


    results <-
      foreach::foreach(i = seq_along(imgs), .packages = "pliman") %dopar%{
        analyze_objects(imgs[[i]],
                        index = index,
                        segment_objects = segment_objects,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        plot = plot,
                        object_size = object_size,
                        object_index = object_index,
                        veins = veins,
                        width_at = width_at,
                        efourier = efourier,
                        invert = invert,
                        watershed = watershed,
                        filter = filter,
                        return_mask = FALSE,
                        ...)
      }
  } else{
    results <-
      lapply(seq_along(imgs), function(i){
        analyze_objects(imgs[[i]],
                        index = index,
                        segment_objects = segment_objects,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        plot = plot,
                        object_size = object_size,
                        object_index = object_index,
                        veins = veins,
                        width_at = width_at,
                        efourier = efourier,
                        invert = invert,
                        watershed = watershed,
                        filter = filter,
                        return_mask = FALSE,
                        ...)
      })
  }
  names(results) <- paste0("shp", 1:length(shapes))

  # RESULTS
  res <-
    do.call(rbind,
            lapply(results, function(x){x$results}))
  vect <- rownames(res)

  res$img <-
    sapply(seq_along(vect),
           function(i){
             strsplit(vect[[i]], split = "\\.")[[1]][[1]]
           })
  res <- res[, c(ncol(res), 1:(ncol(res) - 1))]
  rownames(res) <- NULL

  # STATISTICS
  statistics <-
    do.call(rbind,
            lapply(seq_along(results), function(x){
              transform(results[[x]][["statistics"]], img = names(results[x]))[,c(3, 1, 2)]
            }))

  if(!is.null(results[[1]][["object_rgb"]])){
    object_rgb <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["object_rgb"]], img = names(results[i]))
              }))
    object_rgb <- object_rgb[, c(ncol(object_rgb), 1:(ncol(object_rgb) - 1))]
  } else{
    object_rgb <- NULL
  }

  if(!is.null(results[[1]][["object_index"]])){
    object_index <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["object_index"]], img = names(results[i]))
              }))
    object_index <- object_index[, c(ncol(object_index), 1:(ncol(object_index) - 1))]
  } else{
    object_index <- NULL
  }

  # FOURIER COEFFICIENTS
  if(!isFALSE(efourier)){
    efourier <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["efourier"]],
                          img =  names(results[i]))
              })
      )
    efourier <- efourier[, c(ncol(efourier), 1:ncol(efourier)-1)]
    names(efourier)[2] <- "id"

    efourier_norm <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["efourier_norm"]],
                          img =  names(results[i]))
              })
      )
    efourier_norm <- efourier_norm[, c(ncol(efourier_norm), 1:ncol(efourier_norm)-1)]
    names(efourier_norm)[2] <- "id"


    efourier_error <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["efourier_error"]],
                          img =  names(results[i]))
              })
      )
    efourier_error <- efourier_error[, c(ncol(efourier_error), 1:ncol(efourier_error)-1)]
    names(efourier_error)[2] <- "id"

    efourier_power <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["efourier_power"]],
                          img =  names(results[i]))
              })
      )
    efourier_power <- efourier_power[, c(ncol(efourier_power), 1:ncol(efourier_power)-1)]
    names(efourier_power)[2] <- "id"

    efourier_minharm <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["efourier_minharm"]],
                          img =  names(results[i]))
              })
      )
    efourier_minharm <- efourier_minharm[, c(ncol(efourier_minharm), 1:ncol(efourier_minharm)-1)]
    names(efourier_minharm)[2] <- "id"

  } else{
    efourier <- NULL
    efourier_norm <- NULL
    efourier_error <- NULL
    efourier_power <- NULL
    efourier_minharm <- NULL
  }


  # VEINS FEATURES
  if(isTRUE(veins)){
    veins <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["veins"]],
                          img =  names(results[i]))
              })
      )

    veins <- veins[, c(ncol(veins), 1:ncol(veins)-1)]
  } else{
    veins <- NULL
  }

  # WIDTH AT
  if(isTRUE(width_at)){
    width_at <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["width_at"]],
                          img =  names(results[i]))
              })
      )

    width_at <- width_at[, c(ncol(width_at), 1:ncol(width_at)-1)]
  } else{
    width_at <- NULL
  }

  res[, 1:4] <- correct_coords(res[, 1:4],  nrow(img),  ncol(img), nrow, ncol)
  invisible(
    structure(
      list(results = res,
           statistics = statistics,
           object_rgb = object_rgb,
           object_index = object_index,
           efourier = efourier,
           efourier_norm = efourier_norm,
           efourier_error = efourier_error,
           efourier_power = efourier_power,
           efourier_minharm = efourier_minharm,
           veins = veins,
           width_at = width_at,
           shapefiles = tmp$shapefile,
           mask = mask,
           index = index,
           object_index_computed = object_index_used,
           final_image = img),
      class = "anal_obj"
    )
  )
}

#' Map Object Distances
#'
#' Computes distances between objects in an `anal_obj` object and returns a list
#' of distances, coefficient of variation (CV), and means.
#'
#' @param object An `anal_obj` object computed with `analyze_objects_shp()`.
#' @param by_column The column name in the object's results data frame to group
#'   objects by. Default is "img".
#' @param direction The direction of mapping. Should be one of "horizontal" or
#'   "vertical". Default is "horizontal".
#'
#' @return A list with the following components:
#' \item{distances}{A list of distances between objects grouped by unique values
#' in the specified column/row.}
#' \item{cvs}{A vector of coefficient of variation (CV) values for each column/row.}
#' \item{means}{A vector of mean distances for each column/row.}

#' @seealso \code{\link{analyze_objects_shp}}
#'
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot =TRUE)
#' res <-
#'    analyze_objects_shp(flax,
#'                        nrow = 3,
#'                        ncol = 1,
#'                        watershed = FALSE,
#'                        index = "R/(G/B)",
#'                        plot = FALSE)
#' plot(res$final_image_mask)
#' plot(res$shapefiles)
#'
#' # distance from each leave within each row
#' result <- object_map(res)
#' result$distances
#' result$cvs
#' result$means
#' }
object_map <- function(object,
                       by_column = "img",
                       direction = c("horizontal", "vertical")) {
  optdirec <- c("horizontal", "vertical")
  optdirec <- pmatch(direction[[1]], optdirec)
  if(!inherits(object, "anal_obj") | object$results[1,1] != "shp1"){
    stop("Only objects computed with `analyze_objects_shp()` can be used.")
  }
  coordinates <- object$results[, c(1, 3, 4)]
  unique_values <- unique(coordinates[, by_column])
  distances <- vector("list", length(unique_values))
  for (i in 1:length(unique_values)) {
    subset_coords <- coordinates[coordinates[, by_column] == unique_values[i], 2:3]
    n <- nrow(subset_coords)
    nearest <- order(subset_coords[, optdirec])
    subset_distances <- numeric(n - 1)
    for (j in 1:(n - 1)) {
      x1 <- subset_coords[nearest[j], 1]
      y1 <- subset_coords[nearest[j], 2]
      x2 <- subset_coords[nearest[j+1], 1]
      y2 <- subset_coords[nearest[j+1], 2]
      distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
      subset_distances[j] <- distance
    }
    distances[[i]] <- subset_distances
  }
  if(optdirec == 1){
    names(distances) <- paste0("row", 1:length(distances))
  } else{
    names(distances) <- paste0("column", 1:length(distances))
  }
  cvs <- sapply(distances, function(x){
    (sd(x) / mean(x)) * 100
  })
  means <- sapply(distances, mean)
  invisible(list(distances = distances, cvs = cvs, means = means))
}

#' Mark Object Points
#'
#' Marks the coordinates of objects in an `anal_obj` object on a plot.
#'
#' @param object An `anal_obj` object computed with `analyze_objects_shp()` or
#'   `analyze_objects_shp()`.
#' @param col The color of the marked points. Default is "white".
#'
#' @seealso \code{\link{analyze_objects_shp}}
#' @export
#'
#' @examples
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot =TRUE)
#' res <-
#'    analyze_objects(flax,
#'                        watershed = FALSE,
#'                        index = "R/(G/B)",
#'                        plot = FALSE)
#' object_mark(res)
object_mark <- function(object, col = "white"){
  if(!inherits(object, "anal_obj")){
    stop("Only objects computed with `analyze_objects_shp()` or `analyze_objects_shp()` can be used.")
  }
  coordinates <- object$results[, c("x", "y")]
  points(coordinates, col = col, pch = 16)
}


plot_shp <- function(coords,
                     col_line = "red",
                     size_line = 1,
                     col_text =  "red",
                     size_text = 0.7){
  lapply(seq_along(coords), function(i){
    lines(coords[[i]][, -1], col = col_line, type = "l", lwd = size_line)
    text(min(coords[[i]]$x[-1]), min(coords[[i]]$y[-1]),
         label = i,
         col = col_text,
         cex = size_text,
         adj = c(-0.2, 1.2))
  })
}


#' Plot rectangles colored by a quantitative attribute and overlay on an RGB
#' image
#'
#' This function plots rectangles on top of an RGB image, where each rectangle
#' is colored based on a quantitative variable. The quantitative variable is
#' specified in the `attribute` argument and should be present in the
#' `object_index` of the `object` computed using [analyze_objects_shp()]. The
#' rectangles are colored using a color scale.
#' @inheritParams image_shp
#' @inheritParams plot_index
#' @param object An object computed with [analyze_objects_shp()].
#' @param attribute The name of the quantitative variable in the
#'   \code{object_index} to be used for coloring the rectangles.
#' @param color A vector of two colors to be used for the color scale.
#' @param alpha The transparency level of the rectangles' color (between 0 and 1).
#' @param legend.position The position of the color legend, either
#'   \code{"bottom"} or \code{"right"}.
#' @param na.color The color to be used for rectangles with missing values in
#'   the quantitative variable.
#' @param classes The number of classes in the color scale.
#' @param round The number of decimal places to round the legend values.
#' @param horiz Logical, whether the legend should be horizontal (\code{TRUE})
#'   or vertical (\code{FALSE}).
#' @return The function plots rectangles colored by the specified quantitative
#'   variable on top of the RGB image and shows the continuous color legend
#'   outside the plot.
#' @importFrom grDevices colorRamp
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#'
#' # Computes the DGCI index for each flax leaf
#' flax <- image_pliman("flax_leaves.jpg", plot =TRUE)
#' res <-
#'    analyze_objects_shp(flax,
#'                        buffer_x = 0.2,
#'                        buffer_y = 0.2,
#'                        nrow = 3,
#'                        ncol = 5,
#'                        plot = FALSE,
#'                        object_index = "DGCI")
#' plot(res$final_image)
#' plot_index_shp(res)
#' }
#'
plot_index_shp <- function(object,
                           attribute = "coverage",
                           color = c("red","green"),
                           viewer = c("mapview", "base"),
                           max_pixels = 500000,
                           downsample = NULL,
                           downsample_fun = NULL,
                           alpha = 0.5,
                           legend.position = "bottom",
                           na.color = "gray",
                           classes = 6,
                           round = 3,
                           horiz = TRUE) {
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]
  if(vieweropt == "mapview"){
    quant_var <- get_measures(object)$summary
    quant_var <- quant_var[, c(1, 4:ncol(quant_var))]
    quant_var$img <- gsub("obj", "shp", quant_var$img)
    get_numeric_from_img <- function(x) {
      as.numeric(gsub("shp", "", x))
    }
    quant_var <- quant_var[order(get_numeric_from_img(quant_var$img)), ]
    if(!attribute %in% names(quant_var)){
      stop("Attribute not found. Have you included it in the `object_index` argument from `analyze_objects_shp()`?", call. = FALSE)
    }
    num_rows <- object$shapefiles$nrow
    num_cols <- object$shapefiles$ncol

    list_of_polygons <- lapply(object$shapefiles$shapefiles, function(x) {
      sf::st_polygon(list(as.matrix(x[,2:3])))
    })

    if(num_rows > 1 & num_cols > 1){
      total_elements <- num_rows * num_cols
      num_groups <- floor(total_elements / num_cols)
      start_values <- seq(total_elements, by = -num_cols, length.out = num_groups) - (num_cols - 1)
      list_seq <- list()
      for (i in 1:length(start_values)) {
        list_seq[[i]] <-  c(start_values[[i]], seq(start_values[[i]] + 1, start_values[[i]] + (num_cols - 1)))
      }
      quant_var <- quant_var[unlist(list_seq), ]
      # rownames(quant_var) <- unlist(list_seq)
      quant_var$img <- paste0("shp", 1:nrow(quant_var))

      # Convert the list of polygons to an "sf" data frame
      sf_df <- sf::st_sf(
        geometry = list_of_polygons,
        data = data.frame(img = paste0("shp", unlist(list_seq))),
        crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
      )
    } else{
      if(num_cols == 1){
        quant_var <- quant_var[order(get_numeric_from_img(quant_var$img), decreasing = TRUE), ]
        sf_df <- sf::st_sf(
          geometry = list_of_polygons,
          data = data.frame(img = paste0("shp", nrow(quant_var):1)),
          crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
        )
      }
      if(num_rows == 1){
        quant_var <- quant_var[order(get_numeric_from_img(quant_var$img)), ]
        sf_df <- sf::st_sf(
          geometry = list_of_polygons,
          data = data.frame(img = paste0("shp", 1:nrow(quant_var))),
          crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
        )
      }
    }

    sf_df <- cbind(sf_df, quant_var[, 2:ncol(quant_var)])
    sf_df <- sf::st_transform(sf_df, crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m"))
    mp <-
      mapview::mapview(sf_df,
                       map.types = "OpenStreetMap",
                       zcol = attribute,
                       legend = TRUE,
                       alpha.regions = alpha,
                       layer.name = attribute)
    rgb <- stars::st_as_stars(terra::rast(EBImage::transpose(object$final_image)@.Data[,,1:3]))
    dimsto <- dim(rgb[,,,1])
    nr <- dimsto[1]
    nc <- dimsto[2]
    npix <- nc * nr
    if(max_pixels > 500000){
      message("The number of pixels is too high, which might slow the rendering process.")
    }
    if(npix > max_pixels){
      if(is.null(downsample)){
        compute_downsample <- function(nr, nc, n) {
          if (n == 0) {
            invisible(nr * nc)
          } else if (n == 1) {
            invisible(ceiling(nr/2) * ceiling(nc/2))
          } else if (n > 1) {
            invisible(ceiling(nr/(n+1)) * ceiling(nc/(n+1)))
          } else {
            stop("Invalid downsampling factor. n must be a non-negative integer.")
          }
        }
        possible_downsamples <- 0:100
        possible_npix <- sapply(possible_downsamples, function(x){
          compute_downsample(nr, nc, x)
        })
        downsample <- which.min(abs(possible_npix - max_pixels)) - 1
        message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
      }
      rs <- stars::st_downsample(rgb[,,,1], n = downsample)
      gs <- stars::st_downsample(rgb[,,,2], n = downsample)
      bs <- stars::st_downsample(rgb[,,,3], n = downsample)
      rgb <- terra::rast(c(rs, gs, bs, along = 3)) |> stars::st_as_stars()
    }

    sf::st_crs(rgb) <- sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m ")

    leafem::addStarsRGB(map = mp,
                        x = rgb,
                        r = 1,
                        g = 2,
                        b = 3,
                        maxBytes = 16 * 1024 * 1024,
                        na.color = "#00000000")

  } else{
    if(!is.null(object$object_index)){
      quant_var <- aggregate(. ~ img, data = object[["object_index"]], FUN = mean)
      quant_var <- cbind(quant_var, coverage =  aggregate(coverage ~ img, data = object$results, FUN = sum)$coverage)
    } else{
      quant_var <- aggregate(coverage ~ img, data = object$results, FUN = sum)
    }
    get_numeric_from_img <- function(x) {
      as.numeric(gsub("shp", "", x))
    }
    quant_var <- quant_var[order(get_numeric_from_img(quant_var$img)), ]
    if(!attribute %in% names(quant_var)){
      stop("Attribute not found. Have you included it in the `object_index` argument from `analyze_objects_shp()`?", call. = FALSE)
    }
    quant_variable <- quant_var[, attribute]
    coords_list <- object$shapefiles$shapefiles

    # Combine all rectangles into one data frame for plotting
    all_rectangles <- do.call(rbind, coords_list)

    # Define the xmax, xmin, ymax, ymin of the image
    xmax <- max(all_rectangles$x)
    xmin <- min(all_rectangles$x)
    ymax <- max(all_rectangles$y)
    ymin <- min(all_rectangles$y)

    # Normalize the quantitative variable for color scaling
    rr <- range(quant_variable, na.rm = TRUE)
    svals <- (quant_variable - rr[1]) / diff(rr)
    svals[is.na(svals)] <- 0

    # Create the color ramp function
    f <- colorRamp(color)

    # Calculate colors based on the normalized values and alpha
    valcol <- rgb(f(svals)/255, alpha = alpha)
    valcol[is.na(svals)] <- rgb(t(col2rgb(col = na.color, alpha = FALSE))/255, alpha = alpha)

    # Add RGB image (raster) on the plot
    fin_img <- object$final_image[,,1:3]
    EBImage::colorMode(fin_img) <- "Color"
    plot(fin_img)
    for (i in 1:length(coords_list)) {
      rect(min(coords_list[[i]]$x), min(coords_list[[i]]$y), max(coords_list[[i]]$x), max(coords_list[[i]]$y), col = valcol[i], border = NA)
      rect(min(coords_list[[i]]$x), min(coords_list[[i]]$y), max(coords_list[[i]]$x), max(coords_list[[i]]$y), col = NA, border = "black")
    }

    # Generate the legend
    pos <- round(seq(min(quant_variable, na.rm = TRUE),
                     max(quant_variable, na.rm = TRUE),
                     length.out = classes), round)
    if (any(is.na(quant_variable))) {
      pos <- c(pos, "NA")
    }
    col <- rgb(f(seq(0, 1, length.out = classes))/255, alpha = alpha)
    if (any(is.na(quant_variable))) {
      col <- c(col, rgb(t(col2rgb(col = na.color, alpha = FALSE))/255, alpha = alpha))
    }
    legend(legend.position, title = attribute, legend = pos, fill = col, bty = "n", horiz = horiz)
  }
}


#' Measure disease using shapefiles
#'
#' This function calls [measure_disease()] in each image polygon of a shapefile
#' object generated with [image_shp()] and bind the results into read-ready data
#' frames.
#'
#' @inheritParams measure_disease
#' @inheritParams image_shp
#'
#' @param img The image to be analyzed. Either an image of class `Image` or a
#'   character string containing the image name. In the last, the image will be
#'   searched in the root directory. Declare dir_original to inform a subfolder
#'   that contains the images to be processed.
#' @param nrow,ncol The number of rows and columns to generate the shapefile.
#'   Defaults to `1`.
#' @param prepare Logical value indicating whether to prepare the image for
#'   analysis using [image_prepare()] function. This allows to align and crop
#'   the image before processing. Defaults to `FALSE`.
#' @param dir_original The directory containing the original and processed images.
#'   Defaults to `NULL`. In this case, the function will search for the image `img` in the
#'   current working directory.
#' @param interactive If `FALSE` (default) the grid is created automatically
#'   based on the image dimension and number of rows/columns. If `interactive =
#'   TRUE`, users must draw points at the diagonal of the desired bounding box
#'   that will contain the grid.
#' @param ... Aditional arguments passed on to [measure_disease].
#'
#' @return An object of class `plm_disease_byl`. See more details in the `Value`
#'   section of [measure_disease()].
#' @export
#'
#' @examples
#' if(interactive()){
#' # severity for the three leaflets (from left to right)
#' img <- image_pliman("mult_leaves.jpg", plot = TRUE)
#' sev <-
#'  measure_disease_shp(img = img,
#'                      nrow = 1,
#'                      ncol = 3,
#'                      index_lb = "B",
#'                      index_dh = "NGRDI")
#' sev$severity
#' }

measure_disease_shp <- function(img,
                                nrow = 1,
                                ncol = 1,
                                buffer_x = 0,
                                buffer_y = 0,
                                prepare = FALSE,
                                viewer = "mapview",
                                index_lb = "HUE2",
                                index_dh = "NGRDI",
                                pattern = NULL,
                                threshold = NULL,
                                invert = FALSE,
                                dir_original = NULL,
                                show_features = FALSE,
                                interactive = FALSE,
                                plot = TRUE,
                                parallel = FALSE,
                                workers = NULL,
                                verbose = TRUE,
                                ...){
  if(isTRUE(prepare)){
    img <- image_prepare(img, viewer = viewer)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  ## declare alias for dopar command
  `%dopar%` <- foreach::`%dopar%`
  # helper function
  help_meas_shp <- function(img,
                            nrow,
                            ncol,
                            buffer_x,
                            buffer_y,
                            index_lb,
                            index_dh,
                            threshold,
                            invert,
                            show_features,
                            ...){
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
    tmp <- object_split_shp(img, nrow, ncol,
                            buffer_x = buffer_x,
                            buffer_y = buffer_y,
                            interactive = interactive,
                            only_shp = FALSE)
    imgs <- tmp$imgs
    shapes <- tmp$shapefile$shapefiles

    if(isTRUE(plot)){
      op <- par(mfrow = c(nrow, ncol))
      on.exit(par(op))
    }

    results <-
      lapply(seq_along(imgs), function(i){
        measure_disease(imgs[[i]],
                        name = names(imgs[i]),
                        prefix = "",
                        index_lb = index_lb,
                        index_dh = index_dh,
                        invert = invert,
                        threshold = threshold,
                        show_features = show_features,
                        plot = plot,
                        ...)
      })


    names(results) <- paste0("shp", 1:length(shapes))

    # severity
    res <-
      do.call(rbind,
              lapply(results, function(x){x$severity}))
    vect <- rownames(res)

    res$img <-
      sapply(seq_along(vect),
             function(i){
               strsplit(vect[[i]], split = "\\.")[[1]][[1]]
             })
    res <- res[, c(ncol(res), 1:(ncol(res) - 1))]
    rownames(res) <- NULL

    # shape
    if(!is.null(results$shp1$shape[[1]])){
      shape <-
        do.call(rbind,
                lapply(seq_along(results), function(x){
                  transform(results[[x]][["shape"]], img = names(results[x]))[,c(17, 1:16)]
                }))
      statistics <-
        do.call(rbind,
                lapply(seq_along(results), function(x){
                  transform(results[[x]][["statistics"]], img = names(results[x]))[,c(3, 1, 2)]
                }))
    } else{
      shape <- NULL
      statistics <- NULL
    }

    xycoords <- do.call(rbind,
                        lapply(tmp$shapefile$shapefiles, function(x){
                          coords <- x[, 2:3]
                          x <- mean(c(max(coords[, 1]), min(coords[, 1])))
                          y <- mean(c(max(coords[, 2]), min(coords[, 2])))
                          c(x, y)
                        }))
    res <- cbind(res[, 1], xycoords, res[, 2:3])
    colnames(res) <- c("img", "x", "y",  "healthy", "symptomatic")

    invisible(
      structure(
        list(severity = res,
             shape = shape,
             statistics = statistics,
             shapefiles = tmp$shapefile),
        class = "plm_disease_byl"
      )
    )
  }

  ## apply the function to the image list
  if(missing(pattern)){
    results <- help_meas_shp(img,
                             nrow,
                             ncol,
                             buffer_x,
                             buffer_y,
                             index_lb,
                             index_dh,
                             threshold,
                             invert,
                             show_features,
                             ...)
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
      workers2 <- ifelse(is.null(workers), ceiling(detectCores() * 0.2), workers)
      cl2 <- parallel::makePSOCKcluster(workers2)
      doParallel::registerDoParallel(cl2)
      on.exit(stopCluster(cl2))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",workers2, "). Please wait.")
      }
      results <-
        foreach::foreach(i = seq_along(names_plant), .packages = "pliman") %dopar%{
          help_meas_shp(names_plant[[i]],
                        nrow,
                        ncol,
                        buffer_x,
                        buffer_y,
                        index_lb,
                        index_dh,
                        threshold,
                        invert,
                        show_features,
                        ...)
        }
    } else{
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        if(verbose == TRUE){
          run_progress(pb, actual = i,
                       text = paste("Processing image", names_plant[i]))
        }
        results[[i]] <- help_meas_shp(img  = names_plant[i],
                                      nrow,
                                      ncol,
                                      buffer_x,
                                      buffer_y,
                                      index_lb,
                                      index_dh,
                                      threshold,
                                      invert,
                                      show_features,
                                      ...)
      }
    }
    names(results) <- names_plant
    if(isTRUE(show_features)){
      stats <-
        do.call(rbind,
                lapply(seq_along(results), function(x){
                  transform(results[[x]][["statistics"]],
                            shp = img,
                            img = names(results[x]))
                }))[, c(1, 4, 2, 3)]
      shape <-
        do.call(rbind,
                lapply(seq_along(results), function(x){
                  transform(results[[x]][["shape"]],
                            shp = img,
                            img = names(results[x]))
                }))[, c(1, 18, 2:17)]
    } else{
      shape <- NULL
      stats <- NULL
    }
    severity <-
      do.call(rbind,
              lapply(seq_along(results), function(x){
                transform(results[[x]][["severity"]],
                          shp = img,
                          img = names(results[x]))
              }))[, c(1, 4, 2, 3)]

    results <- list(severity = severity,
                    shape = shape,
                    statistics = stats)
  }
  invisible(structure(
    results, class = "plm_disease_byl"
  ))

}
