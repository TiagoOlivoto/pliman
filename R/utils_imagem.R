#'Combines images to a grid
#'
#'Combines several images to a grid
#' @param ... a comma-separated name of image objects or a list containing image
#'   objects.
#' @param labels A character vector with the same length of the number of
#'   objects in `...` to indicate the plot labels.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param col The color for the plot labels. Defaults to `col = "black"`.
#' @param verbose Shows the name of objects declared in `...` or a numeric
#'   sequence if a list with no names is provided. Set to `FALSE` to supress the
#'   text.
#' @importFrom stats reshape
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A grid with the images in `...`
#' @examples
#' library(pliman)
#'img1 <- image_pliman("sev_leaf.jpg")
#'img2 <- image_pliman("sev_leaf_nb.jpg")
#'image_combine(img1, img2)
image_combine <- function(...,
                          labels = NULL,
                          nrow = NULL,
                          ncol = NULL,
                          col = "black",
                          verbose = TRUE){
  if(is.list(c(...))){
    plots <- as.list(...)
    if(class(plots) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      plots <- lapply(plots, function(x){x[[1]]})
    }
    if(!is.null(labels)){
      names(plots) <- labels
    }
  }else{
    plots <- list(...)
    if(is.null(labels)){
      names(plots) <- unlist(strsplit(gsub("c\\(|\\)",  "", substitute(c(...))), "\\s*(\\s|,)\\s*"))[-1]
    } else{
      names(plots) <- labels
    }
  }
  num_plots <- length(plots)
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
  ifelse(is.null(names(plots)), index <- 1:length(plots), index <- names(plots))
  for(i in 1:length(plots)){
    plot(plots[[i]])
    if(verbose == TRUE){
      dim <- image_dimension(plots[[i]], verbose = FALSE)
      text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = col)
    }
  }
}

#'Import and export images
#'
#'Import images from files and URLs and write images to files, possibly with
#'batch processing.
#' @name utils_image
#' @param image
#' * For `image_import()`, a character vector of file names or URLs.
#' * For `image_export()`, an Image object, an array or a list of images.
#' * For `image_pliman()`, a charactere value specifying the image example. See
#' `?pliman_images` for more details.
#' @param name An string specifying the name of the image. It can be either a
#'   character with the image name (e.g., "img1") or name and extension (e.g.,
#'   "img1.jpg"). If none file extension is provided, the image will be saved as
#'   a *.jpg file.
#' @param prefix A prefix to include in the image name when exporting a list of
#'   images. Defaults to `""`, i.e., no prefix.
#' @param extension When `image` is a list, `extension` can be used to define
#'   the extension of exported files. This will overwrite the file extensions
#'   given in `image`.
#' @param pattern A pattern of file name used to identify images to be imported.
#'   For example, if `pattern = "im"` all images in the current working
#'   directory that the name matches the pattern (e.g., img1.-, image1.-, im2.-)
#'   will be imported as a list. Providing any number as pattern (e.g., `pattern
#'   = "1"`) will select images that are named as 1.-, 2.-, and so on. An error
#'   will be returned if the pattern matches any file that is not supported
#'   (e.g., img1.pdf).
#' @param subfolder Optional character string indicating a subfolder within the
#'   current working directory to save the image(s). If the folder doesn't
#'   exist, it will be created.
#' @param path A character vector of full path names; the default corresponds to
#'   the working directory, [getwd()]. It will overwrite (if given) the path
#'   informed in `image` argument.
#' @param plot Plots the image after importing? Defaults to `FALSE`.
#' @param nrow,ncol Passed on to [image_combine()]. The number of rows and
#'   columns to use in the composite image when `plot = TRUE`.
#' @param ... Alternative arguments passed to the corresponding functions from
#'   the `jpeg`, `png`, and `tiff` packages.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_import()` returns a new `Image` object.
#' * `image_export()` returns an invisible vector of file names.
#' * `image_pliman()` returns a new `Image` object with the example image
#' required. If an empty call is used, the path to the `tmp_images` directory
#' installed with the package is returned.
#' @examples
#' library(pliman)
#' folder <- image_pliman()
#' full_path <- paste0(folder, "/sev_leaf.jpg")
#' (path <- file_dir(full_path))
#' (file <- basename(full_path))
#' image_import(image = full_path)
#' image_import(image = file, path = path)
image_import <- function(image,
                         ...,
                         pattern = NULL,
                         path = NULL,
                         plot = FALSE,
                         nrow = NULL,
                         ncol = NULL){
  check_ebi()
  valid_extens <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF")
  if(!is.null(pattern)){
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    path <- ifelse(is.null(path), getwd(), path)
    imgs <- list.files(pattern = pattern, path)
    extensions <- as.character(sapply(imgs, file_extension))
    all_valid <- extensions %in% valid_extens
    if(any(all_valid == FALSE)){
      warning("'", paste(imgs[which(all_valid == FALSE)], collapse = ", "),
              "' of invalid format ignored.", call. = FALSE)
    }
    imgs <- paste0(path, "/", imgs[all_valid])
    if(length(grep(pattern, imgs)) == 0){
      stop(paste("'", pattern, "' pattern not found in '",
                 paste0(dir)),
           call. = FALSE)
    }
    list_img <-
      lapply(imgs, function(x){
        EBImage::readImage(x)
      })
    names(list_img) <- basename(imgs)
    if(isTRUE(plot)){
      image_combine(list_img, nrow = nrow, ncol = ncol)
    }
    return(list_img)
  } else{
    img_dir <- ifelse(is.null(path), file_dir(image), path)
    all_files <- sapply(list.files(img_dir), file_name)
    img_name <- file_name(image)
    test <- img_name %in% all_files
    if(!any(grepl("http", img_dir, fixed = TRUE)) & !all(test)){
      stop(" '",img_name[which(test == FALSE)],"' not found in ", img_dir[which(test == FALSE)],  call. = FALSE)
    }
    img_name <- paste0(img_dir, "/",img_name , ".", file_extension(image))
    if(length(image) > 1){
      ls <-
        lapply(seq_along(img_name),
               function(x){
                 EBImage::readImage(img_name[x], ...)
               })
      names(ls) <- basename(img_name)
      if(isTRUE(plot)){
        image_combine(ls, nrow = nrow, ncol = ncol)
      }
      return(ls)
    } else{
      img <- EBImage::readImage(img_name, ...)
      if(isTRUE(plot)){
        plot(img)
      }
      return(img)
    }
  }
}

#' @export
#' @name utils_image
image_export <- function(image,
                         name,
                         prefix = "",
                         extension = NULL,
                         subfolder = NULL,
                         ...){
  check_ebi()
  if(class(image) %in% c("binary_list", "index_list",
                         "img_mat_list", "palette_list")){
    image <- lapply(image, function(x){x[[1]]})
  }
  if(class(image) == "segment_list"){
    image <- lapply(image, function(x){x[[1]][[1]]})
  }
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    name <- file_name(names(image))
    extens <- file_extension(names(image))
    if(any(sapply(extens, length)) ==  0 & is.null(extension)){
      extens <- rep("jpg", length(image))
      message("Image(s) exported as *.jpg file(s).")
    }
    if(!is.null(extension)){
      extens <- rep(extension, length(image))
    }
    names(image) <- paste0(name, ".", extens)
    if(!missing(subfolder)){
      dir_out <- paste0(getwd(), "/", subfolder)
      if(dir.exists(dir_out) == FALSE){
        dir.create(dir_out, recursive = TRUE)
      }
      names(image) <- paste0(dir_out, "/", prefix, name, ".", extens)
      a <-
        lapply(seq_along(image), function(i){
          EBImage::writeImage(x = image[[i]], files = names(image[i]), ...)
        })
    } else{
      a <-
        lapply(seq_along(image), function(i){
          EBImage::writeImage(x = image[[i]], files = paste0(prefix, names(image[i])), ...)
        })
    }

  } else{
    filname <- file_name(name)
    extens <- unlist(file_extension(name))
    dir_out <- file_dir(name)
    if(length(extens) ==  1){
      extens <- extens
    } else if(length(extens) ==  0 & is.null(extension)){
      extens <- "jpg"
      message("Image(s) exported as *.jpg file(s).")
    } else if(!is.null(extension)){
      extens <- extension
    }
    if(!missing(subfolder) & nchar(dir_out) == 2){
      dir_out <- paste0("./", subfolder)
    }
    if(dir.exists(dir_out) == FALSE){
      dir.create(dir_out, recursive = TRUE)
    }
    name <- paste0(dir_out, "/", filname, ".", extens)
    EBImage::writeImage(image, name)
  }
}
#' @export
#' @name utils_image
image_pliman <- function(image, plot = FALSE){
  path <- system.file("tmp_images", package = "pliman")
  files <- list.files(path)
  if(!missing(image)){
    if(!image %in% files){
      stop("Image not available in pliman.\nAvaliable images: ", paste(files, collapse = ", "), call. = FALSE)
    }
    im <- image_import(system.file(paste0("tmp_images/", image), package = "pliman"))
    if(isTRUE(plot)){
      plot(im)
    }
    return(im)
  } else{
    path
  }
}

##### Spatial transformations
#'Spatial transformations
#'
#' Performs image rotation and reflection
#' * `image autocrop()` Crops automatically  an image to the area of objects.
#' * `image_crop()` Crops an image to the desired area.
#' * `image_trim()` Remove pixels from the edges of an image (20 by default).
#' * `image_dimension()` Gives the dimension (width and height) of an image.
#' * `image_rotate()` Rotates the image clockwise by the given angle.
#' * `image_horizontal()` Converts (if needed) an image to a horizontal image.
#' * `image_vertical()` Converts (if needed) an image to a vertical image.
#' * `image_hreflect()` Performs horizontal reflection of the `image`.
#' * `image_vreflect()` Performs vertical reflection of the `image`.
#' * `image_resize()` Resize the `image`. See more at [EBImage::resize()].
#' * `image_contrast()` Improve contrast locally by performing adaptive
#' histogram equalization. See more at [EBImage::clahe()].
#' * `image_dilate()` Performs image dilatation. See more at [EBImage::dilate()].
#' * `image_erode()` Performs image erosion. See more at [EBImage::erode()].
#' * `image_opening()` Performs an erosion followed by a dilation. See more at
#' [EBImage::opening()].
#' * `image_closing()` Performs a dilation followed by an erosion. See more at
#' [EBImage::closing()].
#' * `image_filter()` Performs median filtering in constant time. See more at
#' [EBImage::medianFilter()].
#' * `image_blur()` Performs blurring filter of images. See more at
#' [EBImage::gblur()].
#' * `image_skeleton()` Performs image skeletonization.
#' @name utils_transform
#' @param image An image or a list of images of class `Image`.
#' @param index The index to segment the image. See [image_index()] for more
#'   details. Defaults to `"NB"` (normalized blue).
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param edge
#' * for [image_autocrop()] the number of pixels in the edge of the cropped
#' image. If `edge = 0` the image will be cropped to create a bounding rectangle
#' (x and y coordinates) around the image objects.
#' * for [image_trim()], the number of pixels removed from the edges. By
#' default, 20 pixels are removed from all the edges.
#' @param top,bottom,left,right The number of pixels removed from `top`,
#'   `bottom`, `left`, and `right` when using [image_trim()].
#' @param angle The rotation angle in degrees.
#' @param bg_col Color used to fill the background pixels, defaults to `"white"`.
#' @param rel_size The relative size of the resized image. Defaults to 100. For
#'   example, setting `rel_size = 50` to an image of width `1280 x 720`, the new
#'   image will have a size of `640 x 360`.
#' @param width,height
#'  * For `image_resize()` the Width and height of the resized image. These arguments
#'   can be missing. In this case, the image is resized according to the
#'   relative size informed in `rel_size`.
#'  * For `image_crop()` a numeric vector indicating the pixel range (x and y,
#' respectively) that will be maintained in the cropped image, e.g., width =
#' 100:200
#' @param kern An `Image` object or an array, containing the structuring
#'   element. Defaults to a brushe generated with [EBImage::makeBrush()].
#' @param shape A character vector indicating the shape of the brush. Can be
#'   `box`, `disc`, `diamond`, `Gaussian` or `line`. Default is `disc`.
#' @param size
#' * For `image_filter()` is the median filter radius (integer). Defaults to `3`.
#' * For `image_dilate()` and `image_erode()` is an odd number containing the
#' size of the brush in pixels. Even numbers are rounded to the next odd one.
#' The default depends on the image resolution and is computed as the image
#' resolution (megapixels) times 20.
#' @param sigma A numeric denoting the standard deviation of the Gaussian filter
#'   used for blurring. Defaults to `3`.
#' @param cache The the L2 cache size of the system CPU in kB (integer).
#'   Defaults to `512`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param plot If `TRUE` plots the modified image. Defaults to `FALSE`.
#' @param ... Additional arguments passed on to [image_binary()].
#' @md
#' @importFrom parallel detectCores clusterExport makeCluster parLapply
#'   stopCluster
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_skeleton()` returns a binary `Image` object.
#' * All other functions returns a  modified version of `image` depending on the
#' `image_*()` function used.
#' * If `image` is a list, a list of the same length will be returned.
#' @examples
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'plot(img)
#'img <- image_resize(img, 50)
#'img1 <- image_rotate(img, 45)
#'img2 <- image_hreflect(img)
#'img3 <- image_vreflect(img)
#'img4 <- image_vertical(img)
#'image_combine(img1, img2, img3, img4)
image_autocrop <- function(image,
                           index = "NB",
                           edge = 5,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_autocrop, index, edge)
    } else{
      res <- lapply(image, image_autocrop, index, edge)
    }
    return(structure(res, class = "autocrop_list"))
  } else{
    conv_hull <- object_coord(image,
                              index = index,
                              id = NULL,
                              edge = edge,
                              show_image = FALSE)
    segmented <- image[conv_hull$row_min:conv_hull$row_max,
                       conv_hull$col_min:conv_hull$col_max,
                       1:3]
    if(isTRUE(plot)){
      plot(segmented)
    }
    return(segmented)
  }
}
#' @name utils_transform
#' @export
image_crop <- function(image,
                       width = NULL,
                       height = NULL,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_crop, width, height)
    } else{
      res <- lapply(image, image_crop, width, height)
    }
    return(res)
  } else{
    if (!is.null(width) | !is.null(height)) {
      dim <- dim(image)[1:2]
      if (!is.null(width)  & is.null(height)) {
        height <- 1:dim[2]
      }
      if (is.null(width) & !is.null(height)) {
        width <- 1:dim[1]
      }
      if(!is.null(height) & !is.null(width)){
        width <- width
        height <- height
      }
      if (!is.numeric(width) | !is.numeric(height)) {
        stop("Vectors must be numeric.")
      }
      image@.Data <- image@.Data[width, height, ]
    }
    if (is.null(width) & is.null(height)) {
      message("Use the left mouse buttom to crop the image.")
      plot(image)
      cord <- locator(type = "p", n = 2, col = "red", pch = 22)
      w <- round(cord$x[[1]], 0):round(cord$x[[2]], 0)
      h <- round(cord$y[[1]], 0):round(cord$y[[2]], 0)
      cord <- apply(data.frame(do.call(rbind, cord)), 2, round, digits = 0)
      rownames(cord) <- c("width", "height")
      colnames(cord) <- c("min", "max")
      image@.Data <- image@.Data[w, h, ]
      if(isTRUE(verbose)){
        print(cord)
      }
    }
    if (isTRUE(plot)) {
      plot(image)
    }
    return(image)
  }
}
#' @name utils_transform
#' @export
image_dimension <- function(image,
                            parallel = FALSE,
                            workers = NULL,
                            verbose = TRUE){
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <-
        as.data.frame(
          do.call(rbind,
                  parLapply(clust, image,  image_dimension, verbose =  FALSE))
        )
      res <- transform(res, image = rownames(res))[,c(3, 1, 2)]
    } else{
      res <-
        do.call(rbind,
                lapply(image, function(x){
                  dim <- image_dimension(x, verbose = FALSE)
                  data.frame(width = dim[[1]],
                             height = dim[[2]])
                }))
      res <- transform(res, image = rownames(res))[,c(3, 1, 2)]
      rownames(res) <- NULL
    }
    if(verbose == TRUE){
      cat("\n----------------------\n")
      cat("Image dimension\n")
      cat("----------------------\n")
      print(res, row.names = FALSE)
      cat("\n")
    }
    invisible(res)
  } else{
    width <- dim(image)[[1]]
    height <- dim(image)[[2]]
    if(verbose == TRUE){
      cat("\n----------------------\n")
      cat("Image dimension\n")
      cat("----------------------\n")
      cat("Width : ", width, "\n")
      cat("Height: ", height, "\n")
      cat("\n")
    }
    invisible(list(width = width, height = height))
  }
}
#' @name utils_transform
#' @export
image_rotate <- function(image,
                         angle,
                         bg_col = "white",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = TRUE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_rotate, angle, bg_col)
    } else{
      lapply(image, image_rotate, angle, bg_col)
    }
  } else{
    img <- EBImage::rotate(image, angle, bg.col = bg_col)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_horizontal <- function(image,
                             parallel = FALSE,
                             workers = NULL,
                             verbose = TRUE,
                             plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_horizontal)
    } else{
      lapply(image, image_horizontal)
    }
  } else{
    width <- dim(image)[[1]]
    height <- dim(image)[[2]]
    if(width < height){
      img <- EBImage::rotate(image, 90)
    } else{
      img <- image
    }
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_vertical <- function(image,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_vertical)
    } else{
      lapply(image, image_vertical)
    }
  } else{
    width <- dim(image)[[1]]
    height <- dim(image)[[2]]
    if(width > height){
      img <- EBImage::rotate(image, 90)
    } else{
      img <- image
    }
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_hreflect <- function(image,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_hreflect)
    } else{
      lapply(image, image_hreflect)
    }
  } else{
    img <- EBImage::flop(image)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_vreflect <- function(image,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_vreflect)
    } else{
      lapply(image, image_vreflect)
    }
  } else{
    img <- EBImage::flip(image)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' @name utils_transform
#' @export
image_resize <- function(image,
                         rel_size = 100,
                         width,
                         height,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_resize, rel_size)
    } else{
      lapply(image, image_resize, rel_size, width, height)
    }
  } else{
    nrow <- dim(image)[[1]]
    new_row <- nrow * rel_size / 100
    width <- ifelse(missing(width), new_row, width)
    img <- EBImage::resize(image, width, height)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' @name utils_transform
#' @export
image_trim <- function(image,
                       edge = NULL,
                       top  = NULL,
                       bottom  = NULL,
                       left = NULL,
                       right = NULL,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  check_ebi()
  if(is.null(edge) & all(sapply(list(top, bottom, left, right), is.null))){
    edge <- 20
  }
  if(is.null(edge) & !all(sapply(list(top, bottom, left, right), is.null))){
    edge <- 0
  }
  top <- ifelse(is.null(top), edge, top)
  bottom <- ifelse(is.null(bottom), edge, bottom)
  left <- ifelse(is.null(left), edge, left)
  right <- ifelse(is.null(right), edge, right)
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_trim, edge, top, bottom, left, right)
    } else{
      lapply(image, image_trim, edge, top, bottom, left, right)
    }
  } else{
    image <- image[, -c(1:top) ,]
    image <- image[, -c((dim(image)[2] - bottom + 1):dim(image)[2]) ,]
    image <- image[-c((dim(image)[1] - right + 1):dim(image)[1]) ,  ,]
    image <- image[-c(1:left), ,]
    if (isTRUE(plot)) {
      plot(image)
    }
    return(image)
  }
}
#' @name utils_transform
#' @export
image_dilate <- function(image,
                         kern = NULL,
                         size = NULL,
                         shape = "disc",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_dilate, kern, size, shape)
    } else{
      lapply(image, image_dilate, kern, size, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(image)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::dilate(image, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_erode <- function(image,
                        kern = NULL,
                        size = NULL,
                        shape = "disc",
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_erode, kern, size, shape)
    } else{
      lapply(image, image_erode, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(image)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::erode(image, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_opening <- function(image,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_opening, kern, size, shape)
    } else{
      lapply(image, image_opening, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(image)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::opening(image, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_closing <- function(image,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_closing, kern, size, shape)
    } else{
      lapply(image, image_closing, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(image)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::closing(image, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_skeleton <- function(image,
                           kern = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE,
                           ...){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_erode)
    } else{
      lapply(image, image_erode)
    }
  } else{
    if(EBImage::colorMode(image) != 0){
      image <- image_binary(image, ..., resize = FALSE, show_image = FALSE)[[1]]
    }
    s <- matrix(1, nrow(image), ncol(image))
    skel <- matrix(0, nrow(image), ncol(image))
    if(is.null(kern)){
      kern <- suppressWarnings(EBImage::makeBrush(2, shape = "diamond"))
    } else{
      kern <- kern
    }
    while (max(s) == 1) {
      k <- EBImage::opening(image, kern)
      s <- image - k
      skel <- skel | s
      image <- EBImage::erode(image, kern)
    }
    img <- EBImage::Image(skel)
    if (plot == TRUE) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_filter <- function(image,
                         size = 2,
                         cache = 512,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_filter, size, cache)
    } else{
      lapply(image, image_filter, size, cache)
    }
  } else{
    img <- EBImage::medianFilter(image, size, cache)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_blur <- function(image,
                       sigma = 3,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_blur, sigma)
    } else{
      lapply(image, image_blur, sigma)
    }
  } else{
    img <- EBImage::gblur(image, sigma)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_contrast <- function(image,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(image)){
    if(class(image) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      image <- lapply(image, function(x){x[[1]]})
    }
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_contrast)
    } else{
      lapply(image, image_contrast)
    }
  } else{
    get_factors <- function(x) {
      factors <- vector()
      for(i in 1:x) {
        if((x %% i) == 0) {
          factors[i] <- i
        }
      }
      return(factors[!is.na(factors)])
    }
    img_width <- dim(image)[1]
    img_height <- dim(image)[2]
    fx <- get_factors(img_width)
    nx <- suppressWarnings(fx[max(which(fx > 1 & fx < 100))])
    fy <- get_factors(img_height)
    ny <- suppressWarnings(fy[max(which(fy > 1 & fy < 100))])
    testx <- any(fx > 1 & fx < 100) == FALSE
    if(testx == TRUE){
      while(testx == TRUE){
        img_width <- img_width + 1
        fx <- get_factors(img_width)
        testx <- !any(fx > 1 & fx < 100)
        if(any(fx) > 100){
          break
        }
      }
      image <- EBImage::resize(image, w = img_width, h = img_height)
      nx <- suppressWarnings(fx[max(which(fx > 1 & fx < 100))])
    }
    testy <- any(fy > 1 & fy < 100) == FALSE
    if(testy == TRUE){
      while(testy == TRUE){
        img_height <- img_height + 1
        fy <- get_factors(img_height)
        testy <- !any(fy > 1 & fy < 100)
        if(any(fy) > 100){
          break
        }
      }
      image <- EBImage::resize(image, w = img_width, h = img_height)
      ny <- suppressWarnings(fy[max(which(fy > 1 & fy < 100))])
    }
    img <- EBImage::clahe(image, nx = nx, ny = ny, bins = 256)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}


#' Creates a binary image
#'
#' Reduce a color, color near-infrared, or grayscale images to a binary image
#' using a given color channel (red, green blue) or even color indexes. The
#' Otsu's thresholding method (Otsu, 1979) is used to automatically perform
#' clustering-based image thresholding.
#'
#' @param image An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. See the available indexes with
#'   [pliman_indexes()] and [image_index()] for more details.
#' @param my_index User can calculate a different index using the band names,
#'   e.g. `my_index = "R+B/G"`.
#' @param threshold By default (`threshold = "Otsu"`), a threshold value based
#'   on Otsu's method is used to reduce the grayscale image to a binary image.
#'   If a numeric value is informed, this value will be used as a threshold.
#'   Inform any non-numeric value different than "Otsu" to iteratively chosen
#'   the threshold based on a raster plot showing pixel intensity of the index.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value as the percentage of desired resizing. For example, if
#'   `resize = 30`, the resized image will have 30% of the size of original
#'   image.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param invert Inverts the binary image, if desired.
#' @param show_image Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing binary images. The length will depend on the number
#'   of indexes used.
#' @importFrom utils read.csv
#' @examples
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg")
#'image_binary(img, index = c("R, G"))
image_binary <- function(image,
                         index = NULL,
                         my_index = NULL,
                         threshold = "Otsu",
                         resize = 30,
                         fill_hull = FALSE,
                         re = NULL,
                         nir = NULL,
                         invert = FALSE,
                         show_image = TRUE,
                         nrow = NULL,
                         ncol = NULL,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust,
                       image,
                       image_binary,
                       index,
                       my_index,
                       threshold,
                       resize,
                       fill_hull,
                       re,
                       nir,
                       invert,
                       show_image,
                       nrow,
                       ncol)
    } else{
      res <- lapply(image,
                    image_binary,
                    index,
                    my_index,
                    threshold,
                    resize,
                    fill_hull,
                    re,
                    nir,
                    invert,
                    show_image,
                    nrow,
                    ncol)
    }
    return(structure(res, class = "binary_list"))
  } else{
    bin_img <- function(imgs,
                        invert,
                        fill_hull,
                        threshold){
      no_inf <- imgs[!is.infinite(imgs)]
      if(threshold == "Otsu"){
        threshold <- EBImage::otsu(imgs, range = c(min(no_inf, na.rm = TRUE),
                                                   max(no_inf, na.rm = TRUE)))
      } else{
        if(is.numeric(threshold)){
          threshold <- threshold
        } else{
          pixels <- data.frame(imgs@.Data)
          colnames(pixels) <- 1:ncol(pixels)
          pixels$id <- 1:nrow(pixels)
          pixels <-
            reshape(pixels,
                    direction = "long",
                    varying = list(names(pixels)[1:ncol(pixels)-1]),
                    v.names = "value",
                    idvar = "id",
                    timevar = "y",
                    times = names(pixels)[1:ncol(pixels)-1])
          pixels$y <- as.numeric(pixels$y)
          p <-
            levelplot(value ~ id * y,
                      data = pixels,
                      xlab = NULL,
                      ylab = NULL,
                      useRaster = TRUE,
                      col.regions = terrain.colors(300),
                      colorkey = list(interpolate = TRUE,
                                      raster = TRUE))
          plot(p)
          threshold <- readline("Selected threshold: ")
        }
      }
      if(invert == FALSE){
        imgs <- EBImage::Image(imgs < threshold)
      } else{
        imgs <- EBImage::Image(imgs > threshold)
      }
      if(isTRUE(fill_hull)){
        imgs <- EBImage::fillHull(imgs)
      }
      return(imgs)
    }
    imgs <- lapply(image_index(image, index, my_index, resize, re, nir, show_image = FALSE, nrow, ncol),
                   bin_img,
                   invert,
                   fill_hull,
                   threshold)
    if(show_image == TRUE){
      num_plots <- length(imgs)
      if (is.null(nrow) && is.null(ncol)){
        ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
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
      index <- names(imgs)
      for(i in 1:length(imgs)){
        plot(imgs[[i]])
        if(verbose == TRUE){
          dim <- image_dimension(imgs[[i]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
        }
      }
    }
    invisible(imgs)
  }
}

#' Image indexes
#'
#' `image_index()` Builds image indexes using Red, Green, Blue, Red-Edge, and NIR bands.
#'
#' @details
#' The following indexes are available in pliman.
#'
#' * `R` red
#' * `G` green
#' * `B` blue
#' * `NR` normalized red `R/(R+G+B)`.
#' * `NG` normalized green `G/(R+G+B)`
#' * `NB` normalized blue `B/(R+G+B)`
#' * `GB` green blue ratio `G/B`
#' * `RB` red blue ratio `R/B`
#' * `GR` green red ratio `G/R`
#' * `BI` brightness Index `sqrt((R^2+G^2+B^2)/3)`
#' * `BIM` brightness Index 2 `sqrt((R*2+G*2+B*2)/3)`
#' * `SCI` Soil Colour Index `(R-G)/(R+G)`
#' * `GLI` Green leaf index Vis Louhaichi et al. (2001) `(2*G-R-B)/(2*G+R+B)`
#' * `HI` Primary colours Hue Index    (2*R-G-B)/(G-B)
#' * `NDGRI` Normalized green red difference index (Tucker, 1979) `(G-R)/(G+R)`
#' * `NDGBI` Normalized green blue difference index `(G-B)/(G+B)`
#' * `NDRBI` Normalized red blue difference index `(R-B)/(R+B)`
#' * `I`     R+G+B
#' * `S`     ((R+G+B)-3*B)/(R+G+B)
#' * `L`     R+G+B/3
#' * `VARI` A Visible Atmospherically Resistant Index `(G-R)/(G+R-B)`
#' * `HUE` Overall Hue Index `atan(2*(B-G-R)/30.5*(G-R))`
#' * `HUE2`  atan(2*(R-G-R)/30.5*(G-B))
#' * `BGI`   B/G
#' * `GRAY`	`0.299*R + 0.587*G + 0.114*B`
#' * `GLAI` `(25*(G-R)/(G+R-B)+1.25)`
#' * `CI` Coloration Index `(R-B)/R`
#' * `SAT` Overhall Saturation Index `(max(R,G,B) - min(R,G,B)) / max(R,G,B)`
#' * `SHP` Shape Index `2*(R-G-B)/(G-B)`
#' * `RI` Redness Index `R**2/(B*G**3)`
#'
#' @name image_index
#' @param image An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. Use [pliman_indexes()] or the
#'   `details` section to see the available indexes.  Defaults to `NULL`
#'   ((normalized) Red, Green and Blue).  One can also use "RGB" for RGB only,
#'   "NRGB" for normalized RGB, or "all" for all indexes.
#' @param my_index User can calculate a different index using the bands names,
#'   e.g. `my_index = "R+B/G"`.
#' @param resize Resize the image before processing? Defaults to `30`, which
#'   resizes the image to 30% of the original size to speed up image processing.
#'   Set `resize = FALSE` to keep the original size of the image.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param show_image Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing Grayscale images. The length will depend on the
#'   number of indexes used.
#' @examples
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg")
#'image_index(img, index = c("R, NR"))
image_index <- function(image,
                        index = NULL,
                        my_index = NULL,
                        resize = FALSE,
                        re = NULL,
                        nir = NULL,
                        show_image = TRUE,
                        nrow = NULL,
                        ncol = NULL,
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE){
  check_ebi()
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_index, index, my_index, resize, re, nir, show_image, nrow, ncol)
    } else{
      res <- lapply(image, image_index, index, my_index, resize, re, nir, show_image, nrow, ncol)
    }
    return(structure(res, class = "index_list"))
  } else{
    if(resize != FALSE){
      image <- image_resize(image, resize)
    }
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    if(is.null(my_index)){
      ifelse(is.null(index),
             index <- c("R", "G", "B", "NR", "NG", "NB"),
             if(index %in% c("RGB", "NRGB", "all")){
               index <-  switch (index,
                                 RGB = c("R", "G", "B"),
                                 NRGB = c("NR", "NG", "NB"),
                                 all = ind$Index
               )} else{
                 index <- strsplit(index, "\\s*(\\s|,)\\s*")[[1]]
               })
    } else{
      index <- my_index
    }
    nir_ind <- as.character(ind$Index[ind$Band %in% c("RedEdge","NIR")])
    imgs <- list()
    for(i in 1:length(index)){
      indx <- index[[i]]
      if(indx %in% c("R", "G", "B")){
        indx <-
          switch (indx,
                  R = "red",
                  G = "green",
                  B = "blue",
                  GR = "gray"
          )
        imgs[[i]] <- EBImage::channel(image, indx)
      } else{
        if(is.null(my_index) & !indx %in% ind$Index){
          stop(paste("Index '",indx,"' is not available in pliman",sep = ""), call. = FALSE)
        }
        R <- try(image@.Data[,,1], TRUE)
        G <- try(image@.Data[,,2], TRUE)
        B <- try(image@.Data[,,3], TRUE)
        test_band <- any(sapply(list(R, G, B), class) == "try-error")
        if(isTRUE(test_band)){
          stop("At least 3 bands (RGB) are necessary to calculate indices available in pliman.", call. = FALSE)
        }
        if(!is.null(re)|!is.null(nir)){
          if(indx %in% nir_ind & is.null(nir)){
            stop(paste("Index ", indx, " need NIR/RedEdge band to be calculated."), call. = FALSE)
          }
          if(!is.null(re)){
            RE <-  try(image@.Data[,,re], TRUE)
          }
          if(!is.null(nir)){
            NIR <- try(image@.Data[,,nir], TRUE)
          }
          test_nir_ne <- any(lapply(list(RE, NIR), class)  == "try-error" )
          if(isTRUE(test_nir_ne)){
            stop("RE and/or NIR is/are not available in your image.", call. = FALSE)
          }
        }
        if(is.null(my_index)){
          imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==indx]))))
        } else{
          imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(my_index))))
        }
      }
    }
    names(imgs) <- index
    if(show_image == TRUE){
      num_plots <- length(imgs)
      if (is.null(nrow) && is.null(ncol)){
        ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
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
      for(i in 1:length(imgs)){
        plot(imgs[[i]])
        if(verbose == TRUE){
          dim <- image_dimension(imgs[[i]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
        }
      }
    }
    invisible(structure(imgs, class = "image_index"))
  }
}


#' Plots an `image_index` object
#'
#' `plot.image_index()` produces a raster (`type = "raster"`, default) or a
#' density (`type = "density"`) plot of the index values computed with
#' `image_index()`.
#'
#' @name image_index
#' @param x An object of class `image_index`.
#' @param type The type of plot. Use `type = "raster"` (default) to produce a
#'   raster plot showing the intensity of the pixels for each image index or
#'   `type = "density"` to produce a density plot with the pixels' intensity.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param ... Currently not used
#' @method plot image_index
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A `trellis` object containing the distribution of the pixels for each
#'   index.
#' @examples
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#'
#' # resize the image to 30% of the original size
#' ind <- image_index(img, resize = 30, show_image = FALSE)
#' plot(ind)
plot.image_index <- function(x,
                             type = "raster",
                             nrow = NULL,
                             ncol = NULL,
                             ...){
  if(!type %in% c("raster", "density")){
    stop("`type` must be one of the 'raster' or 'density'. ")
  }
  if(type == "density"){
    mat <-
      as.data.frame(
        do.call(cbind,
                lapply(x, function(i){
                  as.vector(i)}
                ))
      )
    colnames(mat) <- names(x)
    mat$id <- rownames(mat)
    if(length(x) == 1){
      mat$Spectrum <- colnames(mat)[1]
      colnames(mat)[1] <- "value"
      a <- mat
    } else{
      a <-
        reshape(data.frame(mat),
                idvar = "id",
                varying = list(1:length(x)),
                times = names(x),
                timevar = "Spectrum",
                v.names = "value",
                direction = "long")
    }
    a$Spectrum <- factor(a$Spectrum, levels = names(x))
    num_plots <- nlevels(a$Spectrum)
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
    p <-
      densityplot(~value | factor(Spectrum),
                  data = a,
                  groups = Spectrum,
                  scales=list(relation="free"),
                  xlab = "Pixel value",
                  layout = c(ncol, nrow),
                  plot.points = FALSE)
    return(p)
  } else{
    get_pixels <- function(x, spectrum){
      pixels <- data.frame(x@.Data)
      colnames(pixels) <- 1:ncol(pixels)
      pixels$id <- 1:nrow(pixels)
      pixels <-
        reshape(pixels,
                direction = "long",
                varying = list(names(pixels)[1:ncol(pixels)-1]),
                v.names = "value",
                idvar = "id",
                timevar = "y",
                times = names(pixels)[1:ncol(pixels)-1])
      pixels$y <- as.numeric(pixels$y)
      pixels$spectrum <- spectrum
      return(pixels)
    }
    pixels <-
      do.call(rbind,
              lapply(seq_along(x), function(i){
                get_pixels(x[i], names(x[i]))
              })
      )
    num_plots <-length(unique(pixels$spectrum))
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
    p <-
      levelplot(value ~ id * y | spectrum,
                layout = c(ncol, nrow),
                data = pixels,
                xlab = NULL,
                ylab = NULL,
                useRaster = TRUE,
                col.regions = terrain.colors(300),
                colorkey = list(interpolate = TRUE,
                                raster = TRUE))
    return(p)
  }
}


#' Image segmentation
#' @description
#' * `image_segment()` reduces a color, color near-infrared, or grayscale images
#' to a segmented image using a given color channel (red, green blue) or even
#' color indexes (See [image_index()] for more details). The Otsu's thresholding
#' method (Otsu, 1979) is used to automatically perform clustering-based image
#' thresholding.
#'
#' * `image_segment_iter()` Provides an iterative image segmentation, returning
#' the proportions of segmented pixels.
#' @param image An image object or a list of image objects.
#' @param index
#'  * For `image_segment()`, a character value (or a vector of characters)
#'  specifying the target mode for conversion to binary image. See the available
#'  indexes with [pliman_indexes()].  See [image_index()] for more details.
#' * For `image_segment_iter()` a character or a vector of characters with the
#' same length of `nseg`. It can be either an available index (described above)
#' or any operation involving the RGB values (e.g., `"B/R+G"`).
#' @param my_index User can calculate a different index using the bands names,
#'   e.g. `my_index = "R+B/G"`.
#' @param threshold By default (`threshold = "Otsu"`), a threshold value based
#'   on Otsu's method is used to reduce the grayscale image to a binary image.
#'   If a numeric value is informed, this value will be used as a threshold.
#'   Inform any non-numeric value different than `"Otsu"` to iteratively chosen
#'   the threshold based on a raster plot showing pixel intensity of the index.
#'   For `image_segmentation_iter()`, use a vector (allows a mixed (numeric and
#'   character) type) with the same length of `nseg`.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param invert Inverts the binary image, if desired. For
#'   `image_segmentation_iter()` use a vector with the same length of `nseg`.
#' @param show_image Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param nseg The number of iterative segmentation steps to be performed.
#' @param ... Additional arguments passed on to `image_segment()`.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @export
#' @name image_segment
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_segment()` returns list containing `n` objects where `n` is the
#' number of indexes used. Each objects contains:
#'    * `image` an image with the RGB bands (layers) for the segmented object.
#'    * `mask` A mask with logical values of 0 and 1 for the segmented image.
#'
#' * `image_segment_iter()` returns a list with (1) a data frame with the
#' proportion of pixels in the segmented images and (2) the segmented images.
#'

#' @examples
#' library(pliman)
#'img <- image_pliman("soybean_touch.jpg", plot = TRUE)
#'image_segment(img, index = c("R, G, B"))
#'
image_segment <- function(image,
                          index = NULL,
                          my_index = NULL,
                          threshold = "Otsu",
                          fill_hull = FALSE,
                          re = NULL,
                          nir = NULL,
                          invert = FALSE,
                          show_image = TRUE,
                          nrow = NULL,
                          ncol = NULL,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  check_ebi()
  if(class(image) == "img_segment"){
    image <- image[[1]][["image"]]
  }
  if(is.list(image)){
    if(!all(sapply(image, class)  %in% c("Image", "img_segment"))){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_segment, index, my_index, threshold, fill_hull, re, nir, invert, show_image, nrow, ncol)
    } else{
      res <- lapply(image, image_segment, index, my_index, threshold, fill_hull, re, nir, invert, show_image, nrow, ncol)
    }
    return(structure(res, class = "segment_list"))
  } else{
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    if(is.null(my_index)){
      ifelse(is.null(index),
             index <- c("R", "G", "B", "NR", "NG", "NB"),
             if(index %in% c("RGB", "NRGB", "all")){
               index <-  switch (index,
                                 RGB = c("R", "G", "B"),
                                 NRGB = c("NR", "NG", "NB"),
                                 all = ind$Index
               )} else{
                 index <- strsplit(index, "\\s*(\\s|,)\\s*")[[1]]
               })
    } else{
      index <- my_index
    }
    imgs <- list()
    for(i in 1:length(index)){
      indx <- index[[i]]
      img2 <- image_binary(image,
                           index = indx,
                           my_index = my_index,
                           threshold = threshold,
                           resize = FALSE,
                           fill_hull = fill_hull,
                           re = re,
                           nir = nir,
                           show_image = FALSE,
                           invert = invert)[[1]]
      ID <- which(img2@.Data == FALSE)
      img <- image
      img@.Data[,,1][ID] <- 1
      img@.Data[,,2][ID] <- 1
      img@.Data[,,3][ID] <- 1
      mask <- image@.Data[,,1]
      mask[ID] <- 1
      mask[!ID] <- 0
      mask <- EBImage::as.Image(mask)
      imgs[[i]] <- list(image = img, mask = mask)
    }
    names(imgs) <- index
    num_plots <- length(imgs)
    if (is.null(nrow) && is.null(ncol)){
      ncol <- ifelse(num_plots == 3, 3, ceiling(sqrt(num_plots)))
      nrow <- ceiling(num_plots/ncol)
    }
    if (is.null(ncol)){
      ncol <- ceiling(num_plots/nrow)
    }
    if (is.null(nrow)){
      nrow <- ceiling(num_plots/ncol)
    }
    if(show_image == TRUE){
      op <- par(mfrow = c(nrow, ncol))
      on.exit(par(op))
      for(i in 1:length(imgs)){
        plot(imgs[[i]][[1]])
        if(verbose == TRUE){
          dim <- image_dimension(imgs[[i]][[1]], verbose = FALSE)
          text(0, dim[[2]]*0.075, index[[i]], pos = 4, col = "red")
        }
      }
    }
    invisible(structure(imgs, class = "img_segment"))
  }
}



#' @export
#' @name image_segment
image_segment_iter <- function(image,
                               nseg = 1,
                               index = NULL,
                               invert = NULL,
                               threshold = NULL,
                               show_image = TRUE,
                               verbose = TRUE,
                               nrow = NULL,
                               ncol = NULL,
                               parallel = FALSE,
                               workers = NULL,
                               ...){
  check_ebi()
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, c("image", "image_segment", "image_combine"))
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      a <- parLapply(clust, image, image_segment_iter, nseg, index, invert, threshold, show_image, verbose, nrow, ncol,  ...)
    } else{
      a <- lapply(image, image_segment_iter, nseg, index, invert, threshold, show_image, verbose, nrow, ncol, ...)
    }
    results <-
      do.call(rbind, lapply(a, function(x){
        x$results
      }))
    images <-
      lapply(a, function(x){
        x$images
      })
    return(list(results = results,
                images = images))
  } else{
    avali_index <- pliman_indexes()
    if(nseg == 1){
      if(is.null(invert)){
        invert <- FALSE
      } else{
        invert <- invert
      }
      if(is.null(threshold)){
        threshold <- "Otsu"
      } else{
        threshold <- threshold
      }
      if(is.null(index)){
        image_segment(image,
                      invert = invert[1],
                      index = "all",
                      ...)
        index <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI")
        my_index <- NULL
      } else{
        index <- index[1]
        if(!index %in% avali_index){
          my_index <- index
          index <- NULL
        } else{
          my_index <- NULL
          index <- index
        }
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      segmented <-
        image_segment(image,
                      index = index,
                      my_index = my_index,
                      threshold = my_thresh,
                      invert = invert[1],
                      show_image = FALSE,
                      ...)
      total <- length(image)
      segm <- length(which(segmented[[1]][["image"]] != 1))
      prop <- segm / total * 100
      results <- data.frame(total = total,
                            segmented = segm,
                            prop = prop)
      imgs <- list(image, segmented[[1]][["image"]])
      if(verbose){
        print(results)
      }
      if(show_image == TRUE){
        image_combine(imgs, ...)
      }
      invisible(list(results = results,
                     images = imgs))
    } else{
      if(is.null(index)){
        image_segment(image,
                      index = "all",
                      ...)
        indx <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI")
        my_index <- NULL
      } else{
        if(length(index) != nseg){
          stop("Length of 'index' must be equal 'nseg'.", call. = FALSE)
        }
        indx <- index[1]
        if(!indx %in% avali_index){
          my_index <- indx
          indx <- NULL
        } else{
          my_index <- NULL
          indx <- indx
        }
      }
      if(is.null(invert)){
        invert <- rep(FALSE, nseg)
      } else{
        invert <- invert
      }
      segmented <- list()
      total <- length(image)
      if(is.null(threshold)){
        threshold <- rep("Otsu", nseg)
      } else{
        threshold <- threshold
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      first <-
        image_segment(image,
                      index = indx,
                      my_index = my_index,
                      invert = invert[1],
                      threshold = my_thresh[1],
                      show_image = FALSE,
                      ...)
      segmented[[1]] <- first
      for (i in 2:(nseg)) {
        if(is.null(index)){
          image_segment(first,
                        index = "all",
                        show_image = TRUE,
                        ncol = ncol,
                        nrow = nrow,
                        ...)
          indx <-
            switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                   "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                   "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                   "GRAY", "GLAI", "SAT", "CI", "SHP", "RI")
          if(is.null(indx)){
            break
          }
          my_index <- NULL
        } else{
          indx <- index[i]
          if(!indx %in% avali_index){
            my_index <- indx
            indx <- NULL
          } else{
            my_index <- NULL
            indx <- indx
          }
        }
        my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[i]))),
                            as.character(threshold[i]),
                            as.numeric(threshold[i]))
        second <-
          image_segment(first,
                        index = indx,
                        my_index = my_index,
                        threshold = my_thresh,
                        invert = invert[i],
                        show_image = FALSE,
                        ...)
        segmented[[i]] <- second
        first <- second
      }
      pixels <-
        rbind(total,
              do.call(rbind,
                      lapply(segmented, function(x){
                        length(which(x[[1]][["image"]] != 1))
                      })
              )
        )
      rownames(pixels) <- NULL
      colnames(pixels) <- "pixels"
      prop <- NULL
      for(i in 2:nrow(pixels)){
        prop[1] <- 100
        prop[i] <- pixels[i] / pixels[i - 1] * 100
      }
      pixels <- data.frame(pixels)
      pixels$percent <- prop
      imgs <- lapply(segmented, function(x){
        x[[1]][["image"]]
      })
      imgs <- c(list(image), imgs)
      names <- paste("seg", 1:length(segmented), sep = "")
      names(imgs) <- c("original", names)
      pixels <- transform(pixels, image = c("original",names))
      pixels <- pixels[,c(3, 1, 2)]
      if(verbose){
        print(pixels)
      }
      if(show_image == TRUE){
        image_combine(imgs, ncol = ncol, nrow = nrow, ...)
      }
      invisible(list(results = pixels,
                     images = imgs))
    }
  }
}



#' Convert an image to numerical matrices
#'
#' Given an object image, converts it into three matrices (RGB) and a data frame
#' where each column corresponds to the RGB values.
#' @param image An image object.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing three matrices (R, G, and B), and a data frame
#'   containing four columns: the name of the image in `image` and the R, G, B
#'   values.
#' @examples
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'dim(img)
#'mat <- image_to_mat(img)
#'dim(mat[[1]])
image_to_mat <- function(image,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_to_mat)
    } else{
      res <- lapply(image, image_to_mat)
    }
    return(structure(res, class = "img_mat_list"))
  } else{
    d <- match.call()
    ncols <- ncol(image@.Data[,,1])
    im <- cbind(c(image@.Data[,,1]), c(image@.Data[,,2]), c(image@.Data[,,3]))
    df_in <- transform(data.frame(im), code = paste(d[["image"]]))[c(4, 1, 2, 3)]
    colnames(df_in) <-  c("CODE", "R", "G", "B")
    rbg <- list(R = matrix(im[, 1], ncol = ncols),
                G = matrix(im[, 2], ncol = ncols),
                B = matrix(im[, 3], ncol = ncols),
                df_in = df_in)
    return(rbg)
  }
}


#' Create image palettes
#'
#' `image_palette()`  creates image palettes by applying the k-means algorithm
#' to the RGB values.
#' @param image An image object.
#' @param npal The number of color palettes.
#' @param filter Performs median filtering. This can be useful to reduce the
#'   noise in produced palettes. Defaults to `TRUE`. See more at
#'   [image_filter()].
#' @param blur Performs blurring filter of palettes?  Defaults to `FALSE`. See
#'   more at [image_blur()].
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return
#' * `image_palette()` returns a list with `npal` color palettes of class `Image`.
#' *
#' @name palettes
#' @export
#' @examples
#' \donttest{
#' library(pliman)
#'img <- image_pliman("sev_leaf_nb.jpg")
#'pal <- image_palette(img, npal = 4)
#'image_combine(pal)
#'
#'
#'# runs only in an iterative section
#' if(FALSE){
#' image_palette_pick(img)
#' }
#'}
image_palette <- function(image,
                          npal,
                          filter = TRUE,
                          blur = FALSE,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  check_ebi()
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, c("image", "image_to_mat"))
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, image, image_palette, npal)
    } else{
      res <- lapply(image, image_palette, npal)
    }
    return(structure(res, class = "palette_list"))
  } else{
    imgRGB <-
      data.frame(R = as.vector(image[,,1]),
                 G = as.vector(image[,,2]),
                 B = as.vector(image[,,3]))
    kms <- kmeans(imgRGB, centers = npal)
    rgbs <- list()
    for (i in 1:npal) {
      ID <- which(kms$cluster == i)
      dim_mat <- trunc(sqrt(length(ID))*0.9)
      ID <- sample(ID[1:dim_mat^2])
      R <- image@.Data[,,1][ID]
      G <- image@.Data[,,2][ID]
      B <- image@.Data[,,3][ID]
      pal <- EBImage::Image(c(R, G, B), dim = c(dim_mat, dim_mat, 3), colormode = "Color")
      if(filter == TRUE){
        pal <- image_filter(pal, size = 2)
      }
      if(blur == TRUE){
        pal <- image_blur(pal, sigma = 1)
      }
      rgbs[[i]] <-  pal
    }
    return(rgbs)
  }
}


#' Utilities for image resolution
#'
#' Provides useful conversions between size (cm), number of pixels (px) and
#' dots per inch (dpi).
#' * [dpi_to_cm()] converts a known dpi value to centimeters.
#' * [cm_to_dpi()] converts a known centimeter values to dpi.
#' * [pixels_to_cm()] converts the number of pixels to centimeters, given a
#' known resolution (dpi).
#' * [cm_to_pixels()] converts a distance (cm) to number of pixels, given a
#' known resolution (dpi).
#' * [distance()] Computes the distance between two points in an image based on
#' the Pythagorean theorem.
#' * [dpi()] An interactive function to compute the image resolution given a
#' known distance informed by the user. See more information in the **Details**
#' section.
#' * [npixels()] returns the number of pixels of an image.
#' @details [dpi()] only run in an interactive section. To compute the image
#'   resolution (dpi) the user must use the left button mouse to create a line
#'   of known distance. This can be done, for example, using a template with
#'   known distance in the image (e.g., `la_leaves.jpg`).
#'
#' @name utils_dpi
#' @param image An image object.
#' @param plot Call a new plot to `image`? Defaults to `TRUE`.
#' @param dpi The image resolution in dots per inch.
#' @param px The number of pixels.
#' @param cm The size in centimeters.
#' @return
#' * [dpi_to_cm()], [cm_to_dpi()], [pixels_to_cm()], and [cm_to_pixels()] return
#' a numeric value or a vector of numeric values if the input data is a vector.
#' * [dpi()] returns the computed dpi (dots per inch) given the known distance
#' informed in the plot.
#' @export
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics locator
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' library(pliman)
#' # Convert  dots per inch to centimeter
#' dpi_to_cm(c(1, 2, 3))
#'
#' # Convert centimeters to dots per inch
#' cm_to_dpi(c(1, 2, 3))
#'
#' # Convert centimeters to number of pixels with resolution of 96 dpi.
#' cm_to_pixels(c(1, 2, 3), 96)
#'
#' # Convert number of pixels to cm with resolution of 96 dpi.
#' pixels_to_cm(c(1, 2, 3), 96)
#'
#' if(isTRUE(interactive())){
#' #### compute the dpi (dots per inch) resolution ####
#' # only works in an interactive section
#' # objects_300dpi.jpg has a known resolution of 300 dpi
#' img <- image_pliman("objects_300dpi.jpg")
#' # Higher square: 10 x 10 cm
#' # 1) Run the function dpi()
#' # 2) Use the left mouse button to create a line in the higher square
#' # 3) Declare a known distance (10 cm)
#' # 4) See the computed dpi
#' dpi(img)
#'
#'
#' img2 <- image_pliman("la_leaves.jpg")
#' # square leaf sample (2 x 2 cm)
#' dpi(img2)
#' }
dpi_to_cm <- function(dpi){
  2.54 / dpi
}
#' @name utils_dpi
#' @export
cm_to_dpi <- function(cm){
  cm / 2.54
}
#' @name utils_dpi
#' @export
pixels_to_cm <- function(px, dpi){
  px * (2.54 / dpi)
}
#' @name utils_dpi
#' @export
cm_to_pixels <- function(cm, dpi){
  cm / (2.54 / dpi)
}
#' @name utils_dpi
#' @export
npixels <- function(image){
  if(class(image) != "Image"){
    stop("Image must be of class 'Image'.")
  }
  dim <- dim(image)
  dim[[1]] * dim[[2]]
}
#' @name utils_dpi
#' @export
dpi <- function(image, plot = TRUE){
  if(isTRUE(interactive())){
    pix <- distance(image, plot = plot)
    known <- as.numeric(readline("known distance (cm): "))
    pix / (known / 2.54)
  }
}
#' @name utils_dpi
#' @export
distance <- function(image, plot = TRUE){
  if(isTRUE(interactive())){
    if(isTRUE(plot)){
      plot(image)
    }
    message("Use the first mouse button to create a line in the plot.")
    coords <- locator(type = "l",
                      n = 2,
                      lwd = 2,
                      col = "red")
    pix <- sqrt((coords$x[1] - coords$x[2])^2 + (coords$y[1] - coords$y[2])^2)
    return(pix)
  }
}



#' Color spaces
#'
#' Convert RGB to LAB color space.
#' @param image An image object.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing the image in the new color space.
#' @examples
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'img2 <- rgb_to_hsv(img)
#'image_combine(img, img2)
rgb_to_hsv <- function(image){
  hsv <- rgb2hsv(r = c(image[,,1]),
                 g = c(image[,,2]),
                 b = c(image[,,3]),
                 maxColorValue = 1)
  img <- EBImage::Image(array(c(hsv[1,], hsv[2,], hsv[3,]),
                              c(dim(image)[1], dim(image)[2], 3)),
                        colormode = "Color")
  return(img)
}
