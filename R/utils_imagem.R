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
#' @param img
#' * For `image_import()`, a character vector of file names or URLs.
#' * For `image_export()`, an Image object, an array or a list of images.
#' * For `image_pliman()`, a charactere value specifying the image example. See
#' `?pliman_images` for more details.
#' @param which logical scalar or integer vector to indicate which image are
#'   imported if a TIFF files is informed. Defaults to `1` (the first image is
#'   returned).
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
#' @param resize Resize the image after importation? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param plot Plots the image after importing? Defaults to `FALSE`.
#' @param nrow,ncol Passed on to [image_combine()]. The number of rows and
#'   columns to use in the composite image when `plot = TRUE`.
#' @param ... Alternative arguments passed to the corresponding functions from
#'   the `jpeg`, `png`, and `tiff` packages.
#' @md
#' @importFrom raster stack as.array
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
#' image_import(img = full_path)
#' image_import(img = file, path = path)
image_import <- function(img,
                         ...,
                         which = 1,
                         pattern = NULL,
                         path = NULL,
                         resize = FALSE,
                         plot = FALSE,
                         nrow = NULL,
                         ncol = NULL){
  # check_ebi()
  valid_extens <- c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF", "TIF", "tif", "gri", "grd")
  if(!is.null(pattern)){
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    path <- ifelse(is.null(path), getwd(), path)
    imgs <- list.files(pattern = pattern, path)
    if(length(grep(pattern, imgs)) == 0){
      stop(paste("'", pattern, "' pattern not found in '",
                 paste0(dir)),
           call. = FALSE)
    }
    extensions <- as.character(sapply(imgs, file_extension))
    all_valid <- extensions %in% valid_extens
    if(any(all_valid == FALSE)){
      warning("'", paste(imgs[which(all_valid == FALSE)], collapse = ", "),
              "' of invalid format ignored.", call. = FALSE)
    }
    imgs <- paste0(path, "/", imgs[all_valid])

    list_img <-
      lapply(imgs, function(x){
        EBImage::readImage(x)
      })
    names(list_img) <- basename(imgs)
    if(isTRUE(plot)){
      image_combine(list_img, nrow = nrow, ncol = ncol)
    }
    if(resize != FALSE){
      if(!is.numeric(resize)){
        stop("Argument `resize` must be numeric.", call. = FALSE)
      }
      list_img <- image_resize(list_img, resize)
    }
    return(list_img)
  } else{
    img_dir <- ifelse(is.null(path), file_dir(img), path)
    all_files <- sapply(list.files(img_dir), file_name)
    img_name <- file_name(img)
    test <- img_name %in% file_name(list.files(img_dir))
    if(!any(grepl("http", img_dir, fixed = TRUE)) & !all(test)){
      stop(" '",img_name[which(test == FALSE)],"' not found in ", img_dir[which(test == FALSE)],  call. = FALSE)
    }
    img_name <- paste0(img_dir, "/",img_name , ".", file_extension(img))
    if(length(img) > 1){
      ls <-
        lapply(seq_along(img_name),
               function(x){
                 if(file_extension(img_name[[1]]) %in% c("tif", "TIF", "tiff", "TIFF", "gri", "grd")){
                   tmp <- raster::stack(img_name[x])
                   nlayers <- length(tmp@layers)
                   dims <- dim(tmp@layers[[1]])
                   list_layers <- raster::as.array(tmp)
                   if(max(list_layers[,,1]) > 1){
                     list_layers <- apply(list_layers, c(1, dim(list_layers)[3]), function(x){x / 255})
                   }
                   img <- EBImage::Image(list_layers, colormode='Color')
                 } else{
                   EBImage::readImage(img_name[x], ...)
                 }

               })
      names(ls) <- basename(img_name)
      if(isTRUE(plot)){
        image_combine(ls, nrow = nrow, ncol = ncol)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          stop("Argument `resize` must be numeric.", call. = FALSE)
        }
        ls <- image_resize(ls, resize)
      }
      return(ls)
    } else{
      if(file_extension(img_name) %in% c("tif", "TIF", "tiff", "TIFF", "gri", "grd")){
        tmp <- raster::stack(img_name)
        tmp_dir <- tempdir()
        png(paste0(tmp_dir, "tmpimg.png"), height = nrow(tmp), width = ncol(tmp))
        raster::plotRGB(tmp, maxpixels = raster::ncell(tmp))
        dev.off()
        img <- EBImage::readImage(paste0(tmp_dir, "tmpimg.png"))
        aaa <- file.remove(paste0(tmp_dir, "tmpimg.png"))
      } else{
        img <- EBImage::readImage(img_name, ...)
      }
      if(isTRUE(plot)){
        plot(img)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          stop("Argument `resize` must be numeric.", call. = FALSE)
        }
        img <- image_resize(img, resize)
      }
      return(img)
    }
  }
}

#' @export
#' @name utils_image
image_export <- function(img,
                         name,
                         prefix = "",
                         extension = NULL,
                         subfolder = NULL,
                         ...){
  check_ebi()
  if(class(img) %in% c("binary_list", "index_list",
                         "img_mat_list", "palette_list")){
    img <- lapply(img, function(x){x[[1]]})
  }
  if(inherits(img, "segment_list")){
    img <- lapply(img, function(x){x[[1]][[1]]})
  }
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    name <- file_name(names(img))
    extens <- file_extension(names(img))
    if(any(sapply(extens, length)) ==  0 & is.null(extension)){
      extens <- rep("jpg", length(img))
      message("Image(s) exported as *.jpg file(s).")
    }
    if(!is.null(extension)){
      extens <- rep(extension, length(img))
    }
    names(img) <- paste0(name, ".", extens)
    if(!missing(subfolder)){
      dir_out <- paste0(getwd(), "/", subfolder)
      if(dir.exists(dir_out) == FALSE){
        dir.create(dir_out, recursive = TRUE)
      }
      names(img) <- paste0(dir_out, "/", prefix, name, ".", extens)
      a <-
        lapply(seq_along(img), function(i){
          EBImage::writeImage(x = img[[i]], files = names(img[i]), ...)
        })
    } else{
      a <-
        lapply(seq_along(img), function(i){
          EBImage::writeImage(x = img[[i]], files = paste0(prefix, names(img[i])), ...)
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
    EBImage::writeImage(img, name)
  }
}
#' @export
#' @name utils_image
image_pliman <- function(img, plot = FALSE){
  path <- system.file("tmp_images", package = "pliman")
  files <- list.files(path)
  if(!missing(img)){
    if(!img %in% files){
      stop("Image not available in pliman.\nAvaliable images: ", paste(files, collapse = ", "), call. = FALSE)
    }
    im <- image_import(system.file(paste0("tmp_images/", img), package = "pliman"))
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
#'
#'
#' @name utils_transform
#' @param img An image or a list of images of class `Image`.
#' @param index The index to segment the image. See [image_index()] for more
#'   details. Defaults to `"NB"` (normalized blue).
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param show How to plot in mapview viewer, either `"rgb"` or `"index"`.
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
#' @param filter Performs median filtering in the binary image. This is useful
#'   to remove noise (like dust) and improve the image autocropping method. See
#'   more at [image_filter()]. Set to `FALSE` to remove the median filtering.
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
#' @param niter The number of iterations to perform in the thinning procedure.
#'   Defaults to 3. Set to `NULL` to iterate until the binary image is no longer
#'   changing.
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
image_autocrop <- function(img,
                           index = "NB",
                           edge = 5,
                           filter = 3,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, img, image_autocrop, index, edge)
    } else{
      res <- lapply(img, image_autocrop, index, edge)
    }
    return(structure(res, class = "autocrop_list"))
  } else{
    conv_hull <- object_coord(img,
                              index = index,
                              id = NULL,
                              edge = edge,
                              plot = FALSE,
                              filter = filter)
    segmented <- img[conv_hull[1]:conv_hull[2],
                       conv_hull[3]:conv_hull[4],
                       1:3]
    if(isTRUE(plot)){
      plot(segmented)
    }
    return(segmented)
  }
}
#' @name utils_transform
#' @export

image_crop <- function(img,
                       width = NULL,
                       height = NULL,
                       viewer = get_pliman_viewer(),
                       show = "rgb",
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, img, image_crop, width, height, viewer)
    } else{
      res <- lapply(img, image_crop, width, height, viewer)
    }
    return(res)
  } else{
    if (!is.null(width) | !is.null(height)) {
      dim <- dim(img)[1:2]
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
      img@.Data <- img@.Data[width, height, ]
    }
    if (is.null(width) & is.null(height)) {
      if(vieweropt == "base"){
        message("Use the left mouse buttom to crop the image.")
        plot(img)
        cord <- locator(type = "p", n = 2, col = "red", pch = 19)
        minw <- min(cord$x[[1]], cord$x[[2]])
        maxw <- max(cord$x[[1]], cord$x[[2]])
        minh <- min(cord$y[[1]], cord$y[[2]])
        maxh <- max(cord$y[[1]], cord$y[[2]])
        w <- round(minw, 0):round(maxw, 0)
        h <- round(minh, 0):round(maxh, 0)
      } else{
        nc <- ncol(img)
        mv <- mv_rectangle(img, show = show)
        w <- round(min(mv[,1]):max(mv[,1]))
        h <- round((min(mv[,2]))):round(max(mv[,2]))
      }
      img@.Data <- img@.Data[w, h, ]
      if(isTRUE(verbose)){
        cat(paste0("width = ", w[1], ":", w[length(w)]), "\n")
        cat(paste0("height = ", h[1], ":", h[length(h)]), "\n")
      }
    }
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}


#' @name utils_transform
#' @export
image_dimension <- function(img,
                            parallel = FALSE,
                            workers = NULL,
                            verbose = TRUE){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <-
        as.data.frame(
          do.call(rbind,
                  parLapply(clust, img,  image_dimension, verbose =  FALSE))
        )
      res <- transform(res, image = rownames(res))[,c(3, 1, 2)]
    } else{
      res <-
        do.call(rbind,
                lapply(img, function(x){
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
    width <- dim(img)[[1]]
    height <- dim(img)[[2]]
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
image_rotate <- function(img,
                         angle,
                         bg_col = "white",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = TRUE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_rotate, angle, bg_col)
    } else{
      lapply(img, image_rotate, angle, bg_col)
    }
  } else{
    img <- EBImage::rotate(img, angle, bg.col = bg_col)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_horizontal <- function(img,
                             parallel = FALSE,
                             workers = NULL,
                             verbose = TRUE,
                             plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_horizontal)
    } else{
      lapply(img, image_horizontal)
    }
  } else{
    width <- dim(img)[[1]]
    height <- dim(img)[[2]]
    if(width < height){
      img <- EBImage::rotate(img, 90)
    } else{
      img <- img
    }
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_vertical <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_vertical)
    } else{
      lapply(img, image_vertical)
    }
  } else{
    width <- dim(img)[[1]]
    height <- dim(img)[[2]]
    if(width > height){
      img <- EBImage::rotate(img, 90)
    } else{
      img <- img
    }
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_hreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_hreflect)
    } else{
      lapply(img, image_hreflect)
    }
  } else{
    img <- EBImage::flop(img)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_vreflect <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_vreflect)
    } else{
      lapply(img, image_vreflect)
    }
  } else{
    img <- EBImage::flip(img)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' @name utils_transform
#' @export
image_resize <- function(img,
                         rel_size = 100,
                         width,
                         height,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_resize, rel_size)
    } else{
      lapply(img, image_resize, rel_size, width, height)
    }
  } else{
    nrow <- dim(img)[[1]]
    new_row <- nrow * rel_size / 100
    width <- ifelse(missing(width), new_row, width)
    img <- EBImage::resize(img, width, height)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' @name utils_transform
#' @export
image_trim <- function(img,
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
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_trim, edge, top, bottom, left, right)
    } else{
      lapply(img, image_trim, edge, top, bottom, left, right)
    }
  } else{
    img <- img[, -c(1:top) ,]
    img <- img[, -c((dim(img)[2] - bottom + 1):dim(img)[2]) ,]
    img <- img[-c((dim(img)[1] - right + 1):dim(img)[1]) ,  ,]
    img <- img[-c(1:left), ,]
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_dilate <- function(img,
                         kern = NULL,
                         size = NULL,
                         shape = "disc",
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_dilate, kern, size, shape)
    } else{
      lapply(img, image_dilate, kern, size, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(img)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::dilate(img, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_erode <- function(img,
                        kern = NULL,
                        size = NULL,
                        shape = "disc",
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        plot = FALSE){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_erode, kern, size, shape)
    } else{
      lapply(img, image_erode, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(img)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::erode(img, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_opening <- function(img,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_opening, kern, size, shape)
    } else{
      lapply(img, image_opening, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(img)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::opening(img, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_closing <- function(img,
                          kern = NULL,
                          size = NULL,
                          shape = "disc",
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE,
                          plot = FALSE){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_closing, kern, size, shape)
    } else{
      lapply(img, image_closing, size, kern, shape)
    }
  } else{
    if(is.null(kern)){
      dim <- dim(img)
      size <- ifelse(is.null(size), round(dim[[1]]*dim[[2]] / 1e06 * 5, 0), size)
      size <- ifelse(size == 0, 2, size)
      kern <- suppressWarnings(EBImage::makeBrush(size, shape = shape))
    } else{
      kern <- kern
    }
    img <- EBImage::closing(img, kern)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' @name utils_transform
#' @export
image_skeleton <- function(img,
                           kern = NULL,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE,
                           ...){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_skeleton)
    } else{
      lapply(img, image_skeleton)
    }
  } else{
    if(EBImage::colorMode(img) != 0){
      img <- help_binary(img, ..., resize = FALSE)
    }
    s <- matrix(1, nrow(img), ncol(img))
    skel <- matrix(0, nrow(img), ncol(img))
    if(is.null(kern)){
      kern <- suppressWarnings(EBImage::makeBrush(2, shape = "diamond"))
    } else{
      kern <- kern
    }
    while (max(s) == 1) {
      k <- EBImage::opening(img, kern)
      s <- img - k
      skel <- skel | s
      img <- EBImage::erode(img, kern)
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
image_thinning <- function(img,
                           niter = 3,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE,
                           ...){
  check_ebi()
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_thinning, niter)
    } else{
      lapply(img, image_thinning, niter)
    }
  } else{
    if(EBImage::colorMode(img) != 0){
      img <- help_binary(img, ..., resize = FALSE)
    }

    if(is.null(niter)){
      li <- sum(img)
      lf <- 1
      while((li - lf) != 0){
        li <- sum(img)
        tin <- help_edge_thinning(img)
        img <- tin
        lf <- sum(img)
      }
    } else{
      for(i in 1:niter){
        tin <- help_edge_thinning(img)
        img <- tin
      }
    }
    img <- EBImage::Image(img)
    if (plot == TRUE) {
      plot(img)
    }
    return(img)
  }
}





#' @name utils_transform
#' @export
image_filter <- function(img,
                         size = 2,
                         cache = 512,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE,
                         plot = FALSE){
  if(size < 2){
    stop("Using `size` < 2 will crash you R section. Please, consider using 2 or more.")
  }
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_filter, size, cache)
    } else{
      lapply(img, image_filter, size, cache)
    }
  } else{
    img <- EBImage::medianFilter(img, size, cache)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_blur <- function(img,
                       sigma = 3,
                       parallel = FALSE,
                       workers = NULL,
                       verbose = TRUE,
                       plot = FALSE){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_blur, sigma)
    } else{
      lapply(img, image_blur, sigma)
    }
  } else{
    img <- EBImage::gblur(img, sigma)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}
#' @name utils_transform
#' @export
image_contrast <- function(img,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE,
                           plot = FALSE){
  if(is.list(img)){
    if(class(img) %in% c("binary_list", "segment_list", "index_list",
                           "img_mat_list", "palette_list")){
      img <- lapply(img, function(x){x[[1]]})
    }
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, img, image_contrast)
    } else{
      lapply(img, image_contrast)
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
    img_width <- dim(img)[1]
    img_height <- dim(img)[2]
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
      img <- EBImage::resize(img, w = img_width, h = img_height)
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
      img <- EBImage::resize(img, w = img_width, h = img_height)
      ny <- suppressWarnings(fy[max(which(fy > 1 & fy < 100))])
    }
    img <- EBImage::clahe(img, nx = nx, ny = ny, bins = 256)
    if (isTRUE(plot)) {
      plot(img)
    }
    return(img)
  }
}

#' Create an `Image` object of a given color
#'
#' image_create() can be used to create an `Image` object with a desired color and size.
#'
#' @param color either a color name (as listed by [grDevices::colors()]), or a hexadecimal
#'   string of the form `"#rrggbb"`.
#' @param width,heigth The width and heigth of the image in pixel units.
#' @param plot Plots the image after creating it? Defaults to `FALSE`.
#'
#' @return An object of class `Image`.
#' @export
#'
#' @examples
#' image_create("red")
#' image_create("#009E73", width = 300, heigth = 100)

image_create <- function(color,
                         width = 200,
                         heigth = 200,
                         plot = FALSE){
  width <- as.integer(width)
  heigth <- as.integer(heigth)
  rgb <- col2rgb(color) / 255
  r <- rep(rgb[1], width*heigth)
  g <- rep(rgb[2], width*heigth)
  b <- rep(rgb[3], width*heigth)
  img <- EBImage::Image(c(r, g, b),
                        dim = c(width, heigth, 3),
                        colormode = "color")
  if(isTRUE(plot)){
    plot(img)
  }
  return(img)
}

#' Creates a binary image
#'
#' Reduce a color, color near-infrared, or grayscale images to a binary image
#' using a given color channel (red, green blue) or even color indexes. The
#' Otsu's thresholding method (Otsu, 1979) is used to automatically perform
#' clustering-based image thresholding.
#'
#' @param img An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. See the available indexes with
#'   [pliman_indexes()] and [image_index()] for more details.
#' @param threshold The theshold method to be used.
#'  * By default (`threshold = "Otsu"`), a threshold value based
#'  on Otsu's method is used to reduce the grayscale image to a binary image. If
#'  a numeric value is informed, this value will be used as a threshold.
#'
#'  * If `threshold = "adaptive"`, adaptive thresholding (Shafait et al. 2008)
#'  is used, and will depend on the `k` and `windowsize` arguments.
#'
#'  * If any non-numeric value different than `"Otsu"` and `"adaptive"` is used,
#'  an iterative section will allow you to choose the threshold based on a
#'  raster plot showing pixel intensity of the index.
#' @param k a numeric in the range 0-1. when `k` is high, local threshold
#'   values tend to be lower. when `k` is low, local threshold value tend to be
#'   higher.
#' @param windowsize windowsize controls the number of local neighborhood in
#'   adaptive thresholding. By default it is set to `1/3 * minxy`, where
#'   `minxy` is the minimum dimension of the image (in pixels).
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If `TRUE`, pixels that have R, G, and B values equals to 1 will be
#'   considered as `NA`. This may be useful to compute an image index for
#'   objects that have, for example, a white background. In such cases, the
#'   background will not be considered for the threshold computation.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value as the percentage of desired resizing. For example, if
#'   `resize = 30`, the resized image will have 30% of the size of original
#'   image.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param filter Performs median filtering in the binary image? (Defaults to
#'   `FALSE`). Provide a positive integer > 1 to indicate the size of the median
#'   filtering. Higher values are more efficient to remove noise in the
#'   background but can dramatically impact the perimeter of objects, mainly for
#'   irregular perimeters such as leaves with serrated edges.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param invert Inverts the binary image, if desired.
#' @param plot Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references
#' Otsu, N. 1979. Threshold selection method from gray-level histograms. IEEE
#' Trans Syst Man Cybern SMC-9(1): 62–66. \doi{10.1109/tsmc.1979.4310076}
#'
#' Shafait, F., D. Keysers, and T.M. Breuel. 2008. Efficient implementation of
#' local adaptive thresholding techniques using integral images. Document
#' Recognition and Retrieval XV. SPIE. p. 317–322 \doi{10.1117/12.767755}
#'
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
#'
image_binary <- function(img,
                         index = NULL,
                         threshold = c("Otsu", "adaptive"),
                         k = 0.1,
                         windowsize = NULL,
                         has_white_bg = FALSE,
                         resize = FALSE,
                         fill_hull = FALSE,
                         filter = FALSE,
                         re = NULL,
                         nir = NULL,
                         invert = FALSE,
                         plot = TRUE,
                         nrow = NULL,
                         ncol = NULL,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE){
  threshold <- threshold[[1]]
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust,
                       img,
                       image_binary,
                       index,
                       threshold,
                       k,
                       windowsize,
                       has_white_bg,
                       resize,
                       fill_hull,
                       filter,
                       re,
                       nir,
                       invert,
                       plot,
                       nrow,
                       ncol)
    } else{
      res <- lapply(img,
                    image_binary,
                    index,
                    threshold,
                    k,
                    windowsize,
                    has_white_bg,
                    resize,
                    fill_hull,
                    filter,
                    re,
                    nir,
                    invert,
                    plot,
                    nrow,
                    ncol)
    }
    return(structure(res, class = "binary_list"))
  } else{
    bin_img <- function(imgs,
                        invert,
                        fill_hull,
                        threshold,
                        filter){
      # adapted from imagerExtra  https://bit.ly/3Wp4pwv
      if(threshold == "adaptive"){
        if(is.null(windowsize)){
          windowsize <- min(dim(imgs)) / 3
          if(windowsize %% 2 == 0){
            windowsize <- as.integer(windowsize + 1)
          }
        }
        if (windowsize <= 2){
          stop("windowsize must be greater than or equal to 3", call. = FALSE)
        }
        if (windowsize %% 2 == 0){
          warning(sprintf("windowsize is even (%d). windowsize will be treated as %d", windowsize, windowsize + 1), call. = FALSE)
          windowsize <- as.integer(windowsize + 1)
        }
        if (windowsize >= dim(imgs)[[1]] || windowsize >= dim(imgs)[[2]]){
          warning("windowsize is too large. Setting to `min(dim(img)) / 3`", call. = FALSE)
          windowsize <- min(dim(imgs)) / 3
        }
        if (k > 1){
          stop("k is out of range. k must be in [0, 1].", call. = FALSE)
        }
        imgs <- EBImage::Image(threshold_adaptive(as.matrix(imgs), k, windowsize, 0.5))
      }

      if(threshold != "adaptive"){
        if(threshold == "Otsu"){
          if(any(is.infinite(imgs)) | any(is.na(imgs))){
            threshold <- help_otsu(imgs@.Data[!is.infinite(imgs@.Data) & !is.na(imgs@.Data)])
          } else{
            threshold <- help_otsu(imgs@.Data)
          }
        } else{
          if(is.numeric(threshold)){
            threshold <- threshold
          } else{
            pixels <- raster::raster(t(imgs@.Data))
            raster::plot(pixels, col = custom_palette(),  axes = FALSE, asp = NA)
            threshold <- readline("Selected threshold: ")
          }
        }
        imgs <- EBImage::Image(imgs < threshold)
      }

      if(invert == TRUE){
        imgs <- 1 - imgs
      }

      imgs[which(is.na(imgs))] <- FALSE
      if(isTRUE(fill_hull)){
        imgs <- EBImage::fillHull(imgs)
      }
      if(is.numeric(filter) & filter > 1){
        imgs <- EBImage::medianFilter(imgs, filter)
      }
      return(imgs)
    }

    imgs <- lapply(image_index(img, index, resize, re, nir, has_white_bg, plot = FALSE, nrow, ncol, verbose = verbose),
                   bin_img,
                   invert,
                   fill_hull,
                   threshold,
                   filter)
    if(plot == TRUE){
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
#' * RGB color space
#' - `R` red
#' - `G` green
#' - `B` blue
#' - `NR` normalized red `R/(R+G+B)`.
#' - `NG` normalized green `G/(R+G+B)`
#' - `NB` normalized blue `B/(R+G+B)`
#' - `GB` green blue ratio `G/B`
#' - `RB` red blue ratio `R/B`
#' - `GR` green red ratio `G/R`
#' - `BI` brightness Index `sqrt((R^2+G^2+B^2)/3)`
#' - `BIM` brightness Index 2 `sqrt((R*2+G*2+B*2)/3)`
#' - `SCI` Soil Colour Index `(R-G)/(R+G)`
#' - `GLI` Green leaf index Vis Louhaichi et al. (2001) `(2*G-R-B)/(2*G+R+B)`
#' - `HI` Primary colours Hue Index    (2*R-G-B)/(G-B)
#' - `NDGRI` Normalized green red difference index (Tucker, 1979) `(G-R)/(G+R)`
#' - `NDGBI` Normalized green blue difference index `(G-B)/(G+B)`
#' - `NDRBI` Normalized red blue difference index `(R-B)/(R+B)`
#' - `I`     R+G+B
#' - `S`     `((R+G+B)-3*B)/(R+G+B)`
#' - `L`     R+G+B/3
#' - `VARI` A Visible Atmospherically Resistant Index `(G-R)/(G+R-B)`
#' - `HUE` Overall Hue Index `atan(2*(B-G-R)/30.5*(G-R))`
#' - `HUE2`  atan(2*(R-G-R)/30.5*(G-B))
#' - `BGI`   B/G
#' - `GRAY`	`0.299*R + 0.587*G + 0.114*B`
#' - `GRAY2` `((R^2.2+(1.5*G)^2.2+(0.6*B)^2.2)/(1+1.5^2.2+0.6^2.2))^1/2.2`
#' - `GLAI` `(25*(G-R)/(G+R-B)+1.25)`
#' - `CI` Coloration Index `(R-B)/R`
#' - `SAT` Overhall Saturation Index `(max(R,G,B) - min(R,G,B)) / max(R,G,B)`
#' - `SHP` Shape Index `2*(R-G-B)/(G-B)`
#' - `RI` Redness Index `R**2/(B*G**3)`
#'
#' * HSB color space
#' * `DGCI` Dark Green Color Index, based on HSB color space `60\*((G - B) / (max(R, G, B) - min(R, G, B)))`
#'
#' * CIE-Lab color space
#' - `L*`: relative luminance `(0.2126 * R + 0.7152 * G + 0.0722 * B)`
#' - `a*`: `0.55*( (R - (0.2126 * R + 0.7152 * G + 0.0722 * B)) / (1.0 - 0.2126))`
#'
#' @name image_index
#' @param img An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. Use [pliman_indexes()] or the
#'   `details` section to see the available indexes.  Defaults to `NULL`
#'   ((normalized) Red, Green and Blue).  One can also use "RGB" for RGB only,
#'   "NRGB" for normalized RGB, or "all" for all indexes. User can also
#'   calculate your own index using the bands names, e.g. `index = "R+B/G"`.
#' @param resize Resize the image before processing? Defaults to `resize =
#'   FALSE`. Use `resize = 50`, which resizes the image to 50% of the original
#'   size to speed up image processing.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If TRUE, pixels that have R, G, and B values equals to 1 will be considered
#'   as NA. This may be useful to compute an image index for objects that have,
#'   for example, a white background. In such cases, the background will not be
#'   considered for the threshold computation.
#' @param plot Show image after processing?
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 70% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param ... Additional arguments passed on to [plot index()].
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
image_index <- function(img,
                        index = NULL,
                        resize = FALSE,
                        re = NULL,
                        nir = NULL,
                        has_white_bg = FALSE,
                        plot = TRUE,
                        nrow = NULL,
                        ncol = NULL,
                        parallel = FALSE,
                        workers = NULL,
                        verbose = TRUE,
                        ...){
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, img, image_index, index, resize, re, nir, has_white_bg, plot, nrow, ncol)
    } else{
      res <- lapply(img, image_index, index, resize, re, nir, has_white_bg, plot, nrow, ncol)
    }
    return(structure(res, class = "index_list"))
  } else{
    if(resize != FALSE){
      img <- image_resize(img, resize)
    }
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")

    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    }else{
      if(index[[1]] %in% c("RGB", "NRGB", "all")){
        index <-  switch (index,
                          RGB = c("R", "G", "B"),
                          NRGB = c("NR", "NG", "NB"),
                          all = ind$Index
        )} else{
          if(length(index) > 1){
            index <- index
          } else{
            index <- strsplit(index, "\\s*(,)\\s*")[[1]]
          }
        }
    }
    nir_ind <- as.character(ind$Index[ind$Band %in% c("RedEdge","NIR")])
    hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])

    R <- try(img@.Data[,,1], TRUE)
    G <- try(img@.Data[,,2], TRUE)
    B <- try(img@.Data[,,3], TRUE)
    test_band <- any(sapply(list(R, G, B), class) == "try-error")

    if(any(index %in% hsb_ind)){
      hsb <- rgb_to_hsb(data.frame(R = c(R), G = c(G), B = c(B)))
      h <- matrix(hsb$h, nrow = nrow(img), ncol = ncol(img))
      s <- matrix(hsb$s, nrow = nrow(img), ncol = ncol(img))
      b <- matrix(hsb$b, nrow = nrow(img), ncol = ncol(img))
    }
    if(isTRUE(test_band)){
      stop("At least 3 bands (RGB) are necessary to calculate indices available in pliman.", call. = FALSE)
    }

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
        imgs[[i]] <- EBImage::channel(img, indx)
      } else{
        if(!indx %in% ind$Index){
          if(isTRUE(verbose)){
            message(paste("Index '",indx,"' is not available. Trying to compute your own index.",sep = ""))
          }
        }
        if(!is.null(re)|!is.null(nir)){
          if(indx %in% nir_ind & is.null(nir)){
            stop(paste("Index ", indx, " need NIR/RedEdge band to be calculated."), call. = FALSE)
          }
          if(!is.null(re)){
            RE <-  try(img@.Data[,,re], TRUE)
          }
          if(!is.null(nir)){
            NIR <- try(img@.Data[,,nir], TRUE)
          }
          test_nir_ne <- any(lapply(list(RE, NIR), class)  == "try-error" )
          if(isTRUE(test_nir_ne)){
            stop("RE and/or NIR is/are not available in your image.", call. = FALSE)
          }
        }
        if(isTRUE(has_white_bg)){
          R[which(R == 1 & G == 1 & B == 1)] <- NA
          G[which(R == 1 & G == 1 & B == 1)] <- NA
          B[which(R == 1 & G == 1 & B == 1)] <- NA
        }

        if(indx %in% ind$Index){
          imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==indx]))))
        } else{
          imgs[[i]] <- EBImage::Image(eval(parse(text = as.character(indx))))
        }
      }
    }
    names(imgs) <- index
    class(imgs) <- "image_index"
    if(plot == TRUE){
      plot_index(imgs, nrow = nrow, ncol = ncol, ...)
    }
    invisible(imgs)
  }
}




#' Plots an `image_index` object
#'
#' Generates a raster or density plot of the index values computed with
#' `image_index()`.

#'
#' @details When `type = "raster"` (default), the function calls [plot_index()]
#' to create a raster plot for each index present in `x`. If `type = "density"`,
#' a for loop is used to create a density plot for each index. Both types of
#' plots can be arranged in a grid controlled by the `ncol` and `nrow`
#' arguments.
#'
#'
#' @name image_index
#' @param x An object of class `image_index`.
#' @param type The type of plot. Use `type = "raster"` (default) to produce a
#'   raster plot showing the intensity of the pixels for each image index or
#'   `type = "density"` to produce a density plot with the pixels' intensity.
#' @param ... Additional arguments passed to [plot_index()] for customization.
#' @method plot image_index
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A `NULL` object
#' @examples
#'
#' # Example for S3 method plot()
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' # compute the index
#' ind <- image_index(img, index = c("R, G, B, NGRDI"), plot = FALSE)
#' plot(ind)
#'
#' # density plot
#' plot(ind, type = "density")
plot.image_index <- function(x,
                             type = c("raster", "density"),
                             nrow = NULL,
                             ncol = NULL,
                             ...){
  typeop <- c("raster", "density")
  typeop <- typeop[pmatch(type[1], typeop)]

  if(!typeop %in% c("raster", "density")){
    stop("`type` must be one of the 'raster' or 'density'. ")
  }
  if(typeop == "density"){
    mat <-
      as.data.frame(
        do.call(cbind,
                lapply(x, function(i){
                  as.vector(i)}
                ))
      )
    mat <- data.frame(mat[sample(1:nrow(mat), 70000, replace = TRUE),])
    colnames(mat) <- names(x)
    num_plots <- ncol(mat)

    if (is.null(nrow) && is.null(ncol)){
      ncols <- ceiling(sqrt(num_plots))
      nrows <- ceiling(num_plots/ncols)
    }
    if (is.null(ncol)){
      ncols <- ceiling(num_plots/nrows)
    }
    if (is.null(nrow)){
      nrows <- ceiling(num_plots/ncols)
    }
    op <- par(mfrow = c(nrows, ncols),
              mar = c(3, 2.5, 3, 3))
    on.exit(par(op))

    for (col in names(mat)) {
      density_data <- density(mat[[col]])  # Calculate the density for the column
      plot(density_data, main = col, col = "red", lwd = 2, xlab = NA, ylab = "Density")  # Create the density plot
    }

  } else{
    plot_index(x, ncol = ncol, nrow = nrow, ...)
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
#'
#' @inheritParams image_binary
#' @param img An image object or a list of image objects.
#' @param index
#'  * For `image_segment()`, a character value (or a vector of characters)
#'  specifying the target mode for conversion to binary image. See the available
#'  indexes with [pliman_indexes()].  See [image_index()] for more details.
#' * For `image_segment_iter()` a character or a vector of characters with the
#' same length of `nseg`. It can be either an available index (described above)
#' or any operation involving the RGB values (e.g., `"B/R+G"`).
#' @param col_background The color of the segmented background. Defaults to
#'   `NULL` (white background).
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If `TRUE`, pixels that have R, G, and B values equals to 1 will be
#'   considered as `NA`. This may be useful to compute an image index for
#'   objects that have, for example, a white background. In such cases, the
#'   background will not be considered for the threshold computation.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param filter Performs median filtering in the binary image? See more at
#'   [image_filter()]. Defaults to `FALSE`. Use a positive integer to define the
#'   size of the median filtering. Larger values are effective at removing
#'   noise, but adversely affect edges.
#' @param re Respective position of the red-edge band at the original image
#'   file.
#' @param nir Respective position of the near-infrared band at the original
#'   image file.
#' @param invert Inverts the binary image, if desired. For
#'   `image_segmentation_iter()` use a vector with the same length of `nseg`.
#' @param plot Show image after processing?
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
#'# adaptive thresholding
#'
image_segment <- function(img,
                          index = NULL,
                          threshold = c("Otsu", "adaptive"),
                          k = 0.1,
                          windowsize = NULL,
                          col_background = NULL,
                          has_white_bg = FALSE,
                          fill_hull = FALSE,
                          filter = FALSE,
                          re = NULL,
                          nir = NULL,
                          invert = FALSE,
                          plot = TRUE,
                          nrow = NULL,
                          ncol = NULL,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  threshold <- threshold[[1]]
  if(inherits(img, "img_segment")){
    img <- img[[1]]
  }
  if(is.list(img)){
    if(!all(sapply(img, class)  %in% c("Image", "img_segment"))){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, img, image_segment, index, threshold, k, windowsize, col_background, has_white_bg, fill_hull, filter, re, nir, invert, plot = plot, nrow, ncol)
    } else{
      res <- lapply(img, image_segment, index, threshold, k, windowsize, col_background, has_white_bg, fill_hull, filter, re, nir, invert, plot = plot, nrow, ncol)
    }
    return(structure(res, class = "segment_list"))
  } else{
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    if(is.null(index)){
      index <- c("R", "G", "B", "NR", "NG", "NB")
    }else{
      if(index %in% c("RGB", "NRGB", "all")){
        index <-  switch (index,
                          RGB = c("R", "G", "B"),
                          NRGB = c("NR", "NG", "NB"),
                          all = ind$Index
        )} else{
          index <- strsplit(index, "\\s*(,)\\s*")[[1]]
        }
    }
    imgs <- list()
    # color for background
    if (is.null(col_background)){
      col_background <- col2rgb("white") / 255
    } else{
      ifelse(is.character(col_background),
             col_background <- col2rgb(col_background) / 255,
             col_background <- col_background / 255)
    }
    for(i in 1:length(index)){
      imgmask <- img
      indx <- index[[i]]
      img2 <- help_binary(img,
                          index = indx,
                          threshold = threshold,
                          k = k,
                          windowsize = windowsize,
                          has_white_bg = has_white_bg,
                          resize = FALSE,
                          fill_hull = fill_hull,
                          filter = filter,
                          re = re,
                          nir = nir,
                          invert = invert)
      ID <- which(img2@.Data == FALSE)
      imgmask@.Data[,,1][ID] <- col_background[1]
      imgmask@.Data[,,2][ID] <- col_background[2]
      imgmask@.Data[,,3][ID] <- col_background[3]
      imgs[[i]] <- imgmask
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
    if(plot == TRUE){
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
    if(length(imgs) == 1){
      invisible(imgs[[1]])
    } else{
      invisible(structure(imgs, class = "img_segment"))
    }
  }
}




#' @export
#' @name image_segment
image_segment_iter <- function(img,
                               nseg = 2,
                               index = NULL,
                               invert = NULL,
                               threshold = NULL,
                               k = 0.1,
                               windowsize = NULL,
                               has_white_bg = FALSE,
                               plot = TRUE,
                               verbose = TRUE,
                               nrow = NULL,
                               ncol = NULL,
                               parallel = FALSE,
                               workers = NULL,
                               ...){
  check_ebi()
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
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
      a <- parLapply(clust, img, image_segment_iter, nseg, index, invert, threshold, has_white_bg, plot, verbose, nrow, ncol,  ...)
    } else{
      a <- lapply(img, image_segment_iter, nseg, index, invert, threshold, has_white_bg, plot, verbose, nrow, ncol, ...)
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
        image_segment(img,
                      invert = invert[1],
                      index = "all",
                      has_white_bg = has_white_bg,
                      ...)
        index <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
      } else{
        index <- index[1]
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      segmented <-
        image_segment(img,
                      index = index,
                      threshold = my_thresh,
                      invert = invert[1],
                      plot = FALSE,
                      has_white_bg = has_white_bg,
                      ...)
      total <- length(img)
      segm <- length(which(segmented != 1))
      prop <- segm / total * 100
      results <- data.frame(total = total,
                            segmented = segm,
                            prop = prop)
      imgs <- list(img, segmented)
      if(verbose){
        print(results)
      }
      if(plot == TRUE){
        image_combine(imgs, ...)
      }
      invisible(list(results = results,
                     images = imgs))
    } else{
      if(is.null(index)){
        image_segment(img,
                      index = "all",
                      ...)
        indx <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                 "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                 "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
      } else{
        if(length(index) != nseg){
          stop("Length of 'index' must be equal 'nseg'.", call. = FALSE)
        }
        indx <- index[1]
      }
      if(is.null(invert)){
        invert <- rep(FALSE, nseg)
      } else{
        invert <- invert
      }
      segmented <- list()
      total <- length(img)
      if(is.null(threshold)){
        threshold <- rep("Otsu", nseg)
      } else{
        threshold <- threshold
      }
      my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[1]))),
                          as.character(threshold[1]),
                          as.numeric(threshold[1]))
      first <-
        image_segment(img,
                      index = indx,
                      invert = invert[1],
                      threshold = my_thresh[1],
                      plot = FALSE,
                      has_white_bg = has_white_bg,
                      ...)
      segmented[[1]] <- first
      for (i in 2:(nseg)) {
        if(is.null(index)){
          image_segment(first,
                        index = "all",
                        plot = TRUE,
                        has_white_bg = has_white_bg,
                        ncol = ncol,
                        nrow = nrow,
                        ...)
          indx <-
            switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                   "R", "G", "B", "NR", "NG", "NB", "GB", "RB", "GR", "BI", "BIM", "SCI", "GLI",
                   "HI", "NGRDI", "NDGBI", "NDRBI", "I", "S", "VARI", "HUE", "HUE2", "BGI", "L",
                   "GRAY", "GLAI", "SAT", "CI", "SHP", "RI", "G-B", "G-R", "R-G", "R-B", "B-R", "B-G", "DGCI", "GRAY2")
          if(is.null(indx)){
            break
          }
        } else{
          indx <- index[i]
        }
        my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold[i]))),
                            as.character(threshold[i]),
                            as.numeric(threshold[i]))
        second <-
          image_segment(first,
                        index = indx,
                        threshold = my_thresh,
                        invert = invert[i],
                        plot = FALSE,
                        ...)
        segmented[[i]] <- second
        first <- second
      }
      pixels <-
        rbind(total,
              do.call(rbind,
                      lapply(segmented, function(x){
                        length(which(x != 1))
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
        x[[1]]
      })
      imgs <- c(list(img), segmented)
      names <- paste("seg", 1:length(segmented), sep = "")
      names(imgs) <- c("original", names)
      pixels <- transform(pixels, image = c("original",names))
      pixels <- pixels[,c(3, 1, 2)]
      if(verbose){
        print(pixels)
      }
      if(plot == TRUE){
        image_combine(imgs, ncol = ncol, nrow = nrow, ...)
      }
      invisible(list(results = pixels,
                     images = imgs))
    }
  }
}



#' Image segmentation using k-means clustering
#'
#' Segments image objects using clustering by the k-means clustering algorithm
#'
#' @param img An `Image` object.
#' @param bands A numeric integer/vector indicating the RGB band used in the
#'   segmentation. Defaults to `1:3`, i.e., all the RGB bands are used.
#' @param nclasses The number of desired classes after image segmentation.
#' @param invert Invert the segmentation? Defaults to `FALSE`. If `TRUE` the
#'   binary matrix is inverted.
#' @param filter Applies a median filtering in the binary matrix? Defaults to
#'   `FALSE`. Use a numeric integer to indicate the size of the median filter.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param plot Plot the segmented image?
#' @return A list with the following values:
#' * `image` The segmented image considering only two classes (foreground and
#' background)
#' * `clusters` The class of each pixel. For example, if `ncluster = 3`,
#' `clusters` will be a two-way matrix with values ranging from 1 to 3.
#' `masks` A list with the binary matrices showing the segmentation.
#' @export
#' @references Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A
#'   K-means clustering algorithm. Applied Statistics, 28, 100–108.
#'   \doi{10.2307/2346830}
#'
#' @examples
#' img <- image_pliman("la_leaves.jpg", plot = TRUE)
#' seg <- image_segment_kmeans(img)
#' seg <- image_segment_kmeans(img, fill_hull = TRUE, invert = TRUE, filter = 10)

image_segment_kmeans <-   function (img,
                                    bands = 1:3,
                                    nclasses = 2,
                                    invert = FALSE,
                                    filter = FALSE,
                                    fill_hull = FALSE,
                                    plot = TRUE){
  imm <- img@.Data[, , bands]
  if(length(dim(imm)) < 3){
    imb <- data.frame(B1 = image_to_mat(imm)[,3])
  } else{
    imb <- image_to_mat(imm)[, -c(1, 2)]
  }
  x <- suppressWarnings(stats::kmeans(imb, nclasses))
  x2 <- x3 <- x$cluster
  nm <- names(sort(table(x2)))
  for (i in 1:length(nm)) {
    x3[x2 == nm[i]] <- i
  }
  m <- matrix(x3, nrow = dim(img)[1])
  LIST <- list()
  for (i in 1:length(nm)) {
    list <- list(m == i)
    LIST <- c(LIST, list)
  }
  if(isTRUE(fill_hull)){
    LIST <- lapply(LIST, EBImage::fillHull)
  }
  if(is.numeric(filter) & filter > 1){
    LIST <- lapply(LIST, EBImage::medianFilter, size = filter)
  }
  mask <- LIST[[1]]
  if(isFALSE(invert)){
    id <- which(mask == 1)
  } else{
    id <- which(mask != 1)
  }
  im2 <- img
  im2@.Data[, , 1][id] <- 1
  im2@.Data[, , 2][id] <- 1
  im2@.Data[, , 3][id] <- 1
  if(isTRUE(plot)){
    if(nclasses == 2){
      plot(im2)
    } else{
      suppressWarnings(image(m, useRaster = TRUE))
    }
  }
  return(list(img = im2,
              clusters = m,
              masks = LIST))
}


#' Image segmentation by hand
#'
#' This R code is a function that allows the user to manually segment an image based on the parameters provided. This only works in an interactive section.
#'
#' @details If the shape is "free", it allows the user to draw a perimeter to
#'   select/remove objects. If the shape is "circle", it allows the user to
#'   click on the center and edge of the circle to define the desired area. If
#'   the shape is "rectangle", it allows the user to select two points to define
#'   the area.
#'
#' @param img An `Image` object.
#' @param shape The type of shape to use. Defaults to "free". Other possible
#'   values are "circle" and "rectangle". Partial matching is allowed.
#' @param type The type of segmentation. By default (`type = "select"`) objects
#'   are selected. Use `type = "remove"` to remove the selected area from the
#'   image.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param resize By default, the segmented object is resized to fill the
#'   original image size. Use `resize = FALSE` to keep the segmented object in
#'   the original scale.
#' @param edge Number of pixels to add in the edge of the segmented object when
#'   `resize = TRUE`. Defaults to 5.
#' @param plot Plot the segmented object? Defaults to `TRUE`.
#'
#' @return A list with the segmented image and the mask used for segmentation.
#' @export
#'
#' @examples
#' if (interactive()) {
#' img <- image_pliman("la_leaves.jpg")
#' seg <- image_segment_manual(img)
#' plot(seg$mask)
#'
#' }
image_segment_manual <-  function(img,
                                  shape = c("free", "circle", "rectangle"),
                                  type = c("select", "remove"),
                                  viewer = get_pliman_viewer(),
                                  resize = TRUE,
                                  edge = 5,
                                  plot = TRUE){
  vals <- c("free", "circle", "rectangle")
  shape <- vals[[pmatch(shape[1], vals)]]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if (isTRUE(interactive())) {
    if(shape == "free"){
      if(vieweropt == "base"){
        plot(img)
        message("Please, draw a perimeter to select/remove objects. Click 'Esc' to finish.")
        stop <- FALSE
        n <- 1e+06
        coor <- NULL
        a <- 0
        while (isFALSE(stop)) {
          if (a > 1) {
            if (nrow(coor) > 1) {
              lines(coor[(nrow(coor) - 1):nrow(coor), 1], coor[(nrow(coor) -
                                                                  1):nrow(coor), 2], col = "red")
            }
          }
          x = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
          if (is.null(x)){
            stop <- TRUE
          }
          coor <- rbind(coor, x)
          a <- a + 1
          if (a >= n) {
            stop = TRUE
          }
        }
        coor <- rbind(coor, coor[1, ])
      } else{
        coor <- mv_polygon(img)
        plot(img)

      }
    }

    if(shape == "circle"){
      if(vieweropt == "base"){
        plot(img)
        message("Click on the center of the circle")
        cent = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
        message("Click on the edge of the circle")
        ext = unlist(locator(type = "p", n = 1, col = "red", pch = 19))
        radius = sqrt(sum((cent - ext)^2))
        x1 = seq(-1, 1, l = 2000)
        x2 = x1
        y1 = sqrt(1 - x1^2)
        y2 = (-1) * y1
        x = c(x1, x2) * radius + cent[1]
        y = c(y1, y2) * radius + cent[2]
      } else{
        mv <- mv_two_points(img)
        radius = sqrt(sum((c(mv$x1, mv$y1) - c(mv$x2, mv$y2))^2))
        x1 = seq(-1, 1, l = 2000)
        x2 = x1
        y1 = sqrt(1 - x1^2)
        y2 = (-1) * y1
        x = c(x1, x2) * radius + mv$x1
        y = c(y1, y2) * radius + mv$y1
      }
      coor = cbind(x, y)
      plot(img)
    }

    if(shape == "rectangle"){
      if(vieweropt == "base"){
        plot(img)
        message("Select 2 points drawing the diagonal that includes the area of interest.")
        cord <- unlist(locator(type = "p", n = 2, col = "red", pch = 19))
        coor <-
          rbind(c(cord[1], cord[3]),
                c(cord[2], cord[3]),
                c(cord[2], cord[4]),
                c(cord[1], cord[4]))
      } else{
        coor <- mv_rectangle(img)
        plot(img)
      }
    }
    mat <- NULL
    for (i in 1:(nrow(coor) - 1)) {
      c1<-  coor[i, ]
      c2 <- coor[i + 1, ]
      a <- c1[2]
      b <- (c2[2] - c1[2])/(c2[1] - c1[1])
      Xs <- round(c1[1], 0):round(c2[1], 0) - round(c1[1], 0)
      Ys <- round(a + b * Xs, 0)
      mat <- rbind(mat, cbind(Xs + round(c1[1], 0), Ys))
      lines(Xs + round(c1[1], 0), Ys, col = "red")
    }
    n = dim(img)
    imF = matrix(0, n[1], n[2])
    id = unique(mat[, 1])
    for (i in id) {
      coorr <- mat[mat[, 1] == i, ]
      imF[i, min(coorr[, 2], na.rm = T):max(coorr[, 2], na.rm = T)] = 1
    }
    mask <- EBImage::fillHull(EBImage::bwlabel(imF))
    # return(mask)
    if(type[1] == "select"){
      id <- mask != 1
    } else{
      id <- mask == 1
    }
    img@.Data[, , 1][id] = 1
    img@.Data[, , 2][id] = 1
    img@.Data[, , 3][id] = 1

    if(isTRUE(resize)){
      nrows <- nrow(mask)
      ncols <- ncol(mask)
      a <- apply(mask, 2, function(x) {
        any(x != 0)
      })
      col_min <- min(which(a == TRUE))
      col_min <- ifelse(col_min < 1, 1, col_min) - edge
      col_max <- max(which(a == TRUE))
      col_max <- ifelse(col_max > ncols, ncols, col_max) + edge
      b <- apply(mask, 1, function(x) {
        any(x != 0)
      })
      row_min <- min(which(b == TRUE))
      row_min <- ifelse(row_min < 1, 1, row_min) - edge
      row_max <- max(which(b == TRUE))
      row_max <- ifelse(row_max > nrows, nrows, row_max) + edge
      img <- img[row_min:row_max, col_min:col_max, 1:3]
    }
    if(isTRUE(plot)){
      plot(img)
    }
    return(list(img = image, mask = EBImage::Image(mask)))
  }
}



#' Convert an image to a data.frame
#'
#' Given an object image, converts it into a data frame where each row corresponds to the intensity values of each pixel in the image.
#' @param img An image object.
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
image_to_mat <- function(img,
                         parallel = FALSE,
                         workers = NULL,
                         verbose = TRUE){
  if(is.list(img)){
    if(!all(sapply(img, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.7), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "img")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      res <- parLapply(clust, img, image_to_mat)
    } else{
      res <- lapply(img, image_to_mat)
    }
    return(structure(res, class = "img_mat_list"))
  } else{
    mat <- cbind(expand.grid(Row = 1:dim(img)[1], Col = 1:dim(img)[2]))
    if(length(dim(img)) == 3){
      for (i in 1:dim(img)[3]) {
        mat <- cbind(mat, c(img[, , i]))
      }
      colnames(mat) = c("row", "col", paste0("B", 1:dim(img)[3]))
    } else{
      mat <- cbind(mat, c(img))
      colnames(mat) = c("row", "col", "B1")
    }
    return(mat)
  }
}


#' Create image palettes
#'
#' `image_palette()`  creates image palettes by applying the k-means algorithm
#' to the RGB values.
#' @param img An image object.
#' @param npal The number of color palettes.
#' @param proportional Creates a joint palette with proportional size equal to
#'   the number of pixels in the image? Defaults to `TRUE`.
#' @param plot Plot the generated palette? Defaults to `TRUE`.
#' @return `image_palette()` returns a list with two elements:
#' * `palette_list` A list with `npal` color palettes of class `Image`.
#' * `joint` An object of class `Image` with the color palettes
#' * `proportions` The proportion of the entire image corresponding to each color in the palette
#' * `rgbs` The average RGB value for each palette
#' @name palettes
#' @export
#' @importFrom stats na.omit
#' @examples
#' \donttest{
#' library(pliman)
#'img <- image_pliman("sev_leaf.jpg")
#'pal <- image_palette(img, npal = 4)
#'
#'image_combine(pal$palette_list)
#'
#'}
image_palette <- function (img,
                           npal = 5,
                           proportional = TRUE,
                           plot = TRUE) {
  id <- matrix(TRUE,
               nrow = nrow(img@.Data[, , 1]),
               ncol = ncol(img@.Data[, , 1]))
  ck <- image_segment_kmeans(img, nclasses = npal, plot = FALSE)[["masks"]]
  layers = length(ck)
  ck2 <- 1 * ck[[1]]
  for (i in 2:layers) {
    ck2 <- ck2 + i * ck[[i]]
  }
  ck <- ck2
  MAT <- NULL
  for (i in unique(na.omit(c(ck)))) {
    r = mean(img@.Data[, , 1][ck == i], na.rm = T)
    g = mean(img@.Data[, , 2][ck == i], na.rm = T)
    b = mean(img@.Data[, , 3][ck == i], na.rm = T)
    MAT = cbind(MAT, c(r = r, g = g, b = b))
  }
  pal_list <- list()
  pal_rgb <- list()
  for(i in 1:ncol(MAT)){
    R <- matrix(rep(MAT[[1, i]], 10000), 100, 100)
    G <- matrix(rep(MAT[[2, i]], 10000), 100, 100)
    B <- matrix(rep(MAT[[3, i]], 10000), 100, 100)
    pal_list[[paste0("pal_", i)]] <- EBImage::rgbImage(R, G, B)
    pal_rgb[[paste0("pal_", i)]] <- c(R = R[1], G = G[1], B = B[1])
  }
  MATn <- NULL
  for (i in unique(na.omit(c(ck)))) {
    r <- length(na.omit(img@.Data[, , 1][ck == i]))
    MATn <- cbind(MATn, r = r)
  }
  props <- data.frame(class = paste0("c", 1:length(MATn)),
                      pixels = t(MATn),
                      prop = t(MATn/sum(MATn)))
  rownames(props) <- NULL
  if (proportional == FALSE) {
    n <- ncol(MAT)
    ARR <- array(NA, dim = c(150, 66 * n, 3))
    c = 1
    f = 66
    for (i in 1:n) {
      ARR[1:150, c:f, 1] <- MAT[1, i]
      ARR[1:150, c:f, 2] <- MAT[2, i]
      ARR[1:150, c:f, 3] <- MAT[3, i]
      c = f + 1
      f = f + 66
    }
  }
  if (proportional == TRUE) {
    n <- ncol(MAT)
    ARR <- array(NA, dim = c(150, 66 * n, 3))
    nn <- round(66 * n * (MATn/sum(MATn)), 0)
    a <- 1
    b <- nn[1]
    nn <- c(nn, 0)
    for (i in 1:n) {
      ARR[1:150, a:b, 1] <- MAT[1, i]
      ARR[1:150, a:b, 2] <- MAT[2, i]
      ARR[1:150, a:b, 3] <- MAT[3, i]
      a <- b + 1
      b <- b + nn[i + 1]
      if (b > (66 * n)) {
        b <- 66 * n
      }
    }
  }
  im2 <- EBImage::as.Image(ARR)
  EBImage::colorMode(im2) <- 2
  if (plot == TRUE) {
    plot(im2)
  }
  return(list(palette_list = pal_list,
              joint = im2,
              proportions = props,
              rgbs = pal_rgb))
}








#' Expands an image
#'
#' Expands an image towards the left, top, right, or bottom by sampling pixels
#' from the image edge. Users can choose how many pixels (rows or columns) are
#' sampled and how many pixels the expansion will have.
#'
#' @param img An `Image` object.
#' @param left,top,right,bottom The number of pixels to expand in the left, top,
#'   right, and bottom directions, respectively.
#' @param edge The number of pixels to expand in all directions. This can be
#'   used to avoid calling all the above arguments
#' @param sample_left,sample_top,sample_right,sample_bottom The number of pixels
#'   to sample from each side. Defaults to 20.
#' @param random Randomly sampling of the edge's pixels? Defaults to `FALSE`.
#' @param filter Apply a median filter in the sampled pixels? Defaults to
#'   `FALSE`.
#' @param plot Plots the extended image? defaults to `FALSE`.
#'
#' @return An `Image` object
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' image_expand(img, left = 200)
#' image_expand(img, right = 150, bottom = 250, filter = 5)
#'
image_expand <- function(img,
                         left = NULL,
                         top = NULL,
                         right = NULL,
                         bottom = NULL,
                         edge = NULL,
                         sample_left = 10,
                         sample_top = 10,
                         sample_right = 10,
                         sample_bottom = 10,
                         random = FALSE,
                         filter = NULL,
                         plot = TRUE){
  if(!is.null(edge)){
    left <- edge
    top <- edge
    right <- edge
    bottom <- edge
  }
  if(sample_left < 2){
    warning("`sample_left` must be > 1. Setting to 2", call. = FALSE)
    sample_left <- 2
  }
  if(sample_top < 2){
    warning("`sample_top` must be > 1. Setting to 2", call. = FALSE)
    sample_top <- 2
  }
  if(sample_right < 2){
    warning("`sample_right` must be > 1. Setting to 2", call. = FALSE)
    sample_right <- 2
  }
  if(sample_bottom < 2){
    warning("`sample_bottom` must be > 1. Setting to 2", call. = FALSE)
    sample_bottom <- 2
  }
  if(!is.null(left)){
    left_img <- img@.Data[1:sample_left,,] |> EBImage::Image(colormode = "Color")
    left_img <- EBImage::resize(left_img, w = left, h = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(left_img)
      for (i in 1:nc[1]) {
        left_img@.Data[i,,] <- left_img@.Data[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      left_img <- EBImage::medianFilter(left_img, size = filter)
    }
    img <- EBImage::abind(left_img, img, along = 1)
  }
  if(!is.null(top)){
    top_img <- img@.Data[,1:sample_top,] |> EBImage::Image(colormode = "Color")
    top_img <- EBImage::resize(top_img, w = dim(img)[1], h = top)
    if(isTRUE(random)){
      nc <- dim(top_img)
      for (i in 1:nc[2]) {
        top_img@.Data[,i,] <- top_img@.Data[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      top_img <- EBImage::medianFilter(top_img, size = filter)
    }
    img <- EBImage::abind(top_img, img, along = 2)
  }
  if(!is.null(right)){
    dimx <- dim(img)[1]
    right_img <- img@.Data[(dimx-sample_right):dimx,,] |> EBImage::Image(colormode = "Color")
    right_img <- EBImage::resize(right_img, w = right, h = dim(img)[2])
    if(isTRUE(random)){
      nc <- dim(right_img)
      for (i in 1:nc[1]) {
        right_img@.Data[i,,] <- right_img@.Data[i,sample(1:nc[2], nc[2]),]
      }
    }
    if(!is.null(filter)){
      right_img <- EBImage::medianFilter(right_img, size = filter)
    }
    img <- EBImage::abind(img, right_img, along = 1)
  }
  if(!is.null(bottom)){
    dimy <- dim(img)[2]
    bot_img <- img@.Data[,(dimy-sample_bottom):dimy,] |> EBImage::Image(colormode = "Color")
    bot_img <- EBImage::resize(bot_img, w = dim(img)[1], h = bottom)
    if(isTRUE(random)){
      nc <- dim(bot_img)
      for (i in 1:nc[2]) {
        bot_img@.Data[,i,] <- bot_img@.Data[sample(1:nc[1], nc[1]),i,]
      }
    }
    if(!is.null(filter)){
      bot_img <- EBImage::medianFilter(bot_img, size = filter)
    }
    img <- EBImage::abind(img, bot_img, along = 2)
  }
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}


#' Squares an image
#'
#' Converts a rectangular image into a square image by expanding the
#' rows/columns using [image_expand()].
#'
#' @inheritParams image_expand
#'
#' @return The modified `Image` object.
#' @param ... Further arguments passed on to [image_expand()].
#' @export
#'
#' @examples
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' dim(img)
#' square <- image_square(img)
#' dim(square)
image_square <- function(img, plot = TRUE, ...){
  len <- dim(img)
  n <- max(len[1], len[2])
  if (len[1] > len[2]) {
    ni1 <- ceiling((n - len[2])/2)
    if((ni1*2 + len[2]) != n){
      ni2 <- ni1 - 1
    } else{
      ni2 <- ni1
    }
    img <- image_expand(img, bottom = ni1, top = ni2, plot = FALSE, ...)
  }
  if (len[2] > len[1]) {
    ni1 <- ceiling((n - len[1])/2)
    if((ni1*2 + len[1]) != n){
      ni2 <- ni1 - 1
    } else{
      ni2 <- ni1
    }
    img <- image_expand(img, left = ni1, right = ni2, plot = FALSE, ...)
  }
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
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
#' @param img An image object.
#' @param dpi The image resolution in dots per inch.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param px The number of pixels.
#' @param cm The size in centimeters.
#' @return
#' * [dpi_to_cm()], [cm_to_dpi()], [pixels_to_cm()], and [cm_to_pixels()] return
#' a numeric value or a vector of numeric values if the input data is a vector.
#' * [dpi()] returns the computed dpi (dots per inch) given the known distance
#' informed in the plot.
#' @export
#' @importFrom grDevices rgb2hsv convertColor
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
npixels <- function(img){
  if(!inherits(img, "Image")){
    stop("Image must be of class 'Image'.")
  }
  dim <- dim(img)
  dim[[1]] * dim[[2]]
}
#' @name utils_dpi
#' @export
dpi <- function(img,
                viewer = get_pliman_viewer()){
  if(isTRUE(interactive())){
    pix <- distance(img, viewer = viewer)
    known <- as.numeric(readline("known distance (cm): "))
    pix / (known / 2.54)
  }
}

#' @name utils_dpi
#' @export
distance <- function(img,
                     viewer = get_pliman_viewer()){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(isTRUE(interactive())){
    if(vieweropt == "base"){
      plot(img)
      message("Use the first mouse button to create a line in the plot.")
      coords <- locator(type = "l",
                        n = 2,
                        lwd = 2,
                        col = "red")
      pix <- sqrt((coords$x[1] - coords$x[2])^2 + (coords$y[1] - coords$y[2])^2)
    } else{
      coords2 <- mv_two_points(img)
      pix <- sqrt((coords2$x1 - coords2$x2)^2 + (coords2$y1 - coords2$y2)^2)
    }
    return(pix)
  }
}



#' Convert between colour spaces
#' @description
#'  `rgb_to_srgb()` Transforms colors from RGB space (red/green/blue) to
#'  Standard Red Green Blue (sRGB), using a gamma correction of 2.2.
#' * `rgb_to_hsb()` Transforms colors from RGB space (red/green/blue) to HSB
#' space (hue/saturation/brightness).
#' * `rgb_to_lab()` Transforms colors from RGB space (red/green/blue) to
#' CIE-LAB space
#'
#' It is assumed that
#' @param object An `Image` object, an object computed with
#'   `analyze_objects()` with a valid `object_index` argument, or a
#'   `data.frame/matrix`. For the last, a three-column data (R, G, and B, respectively)
#'   is required.
#' @export
#' @name utils_colorspace
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A data frame with the columns of the converted color space
#' @examples
#' if(interactive()){
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' rgb_to_lab(img)
#'
#' # analyze the object and convert the pixels
#' anal <- analyze_objects(img, object_index = "B", pixel_level_index = TRUE)
#' rgb_to_lab(anal)
#' }
rgb_to_hsb <- function(object){
  if (any(class(object) %in%  c("data.frame", "matrix"))){
    hsb <-
      rgb_to_hsb_help(r = object[,1],
                      g = object[,2],
                      b = object[,3])
    colnames(hsb) <- c("h", "s", "b")
  }
  if (any(class(object)  %in% c("anal_obj", "anal_obj_ls"))){
    if(!is.null(object$object_rgb)){
      tmp <- object$object_rgb
      if ("img" %in% colnames(tmp)){
        hsb <-
          rgb_to_hsb_help(r = c(tmp[,3]),
                          g = c(tmp[,4]),
                          b = c(tmp[,5]))
        hsb <- data.frame(cbind(tmp[,1:2], hsb))
        colnames(hsb)[1:2] <- c("img", "id")
        colnames(hsb)[3:5] <- c("h", "s", "b")
      }
      hsb <-
        rgb_to_hsb_help(r = c(tmp[,2]),
                        g = c(tmp[,3]),
                        b = c(tmp[,4]))
      hsb <- data.frame(cbind(tmp[,1], hsb))
      colnames(hsb)[1] <- "id"
      colnames(hsb)[2:4] <- c("h", "s", "b")
    } else{
      stop("Cannot obtain the RGB for each object since `object_index` argument was not used. \nHave you accidentally missed the argument `pixel_level_index = TRUE`?")
    }
  }
  if (any(class(object) == "Image")){
    hsb <-
      rgb_to_hsb_help(r = c(object[,,1]),
                      g = c(object[,,2]),
                      b = c(object[,,3]))
    colnames(hsb) <- c("h", "s", "b")
  }
  return(data.frame(hsb))
}

#' @export
#' @name utils_colorspace
rgb_to_srgb <- function(object){
  if (any(class(object) %in%  c("data.frame", "matrix"))){
    srgb <- rgb_to_srgb_help(object[, 1:3])
    colnames(srgb) <- c("sR", "sG", "sB")
  }
  if (any(class(object)  %in% c("anal_obj", "anal_obj_ls"))){
    if(!is.null(object$object_rgb)){
      tmp <- object$object_rgb
      if ("img" %in% colnames(tmp)){
        srgb <- rgb_to_srgb_help(as.matrix(tmp[, 3:5]))
        srgb <- data.frame(cbind(tmp[,1:2], srgb))
        colnames(srgb)[1:2] <- c("img", "id")
        colnames(srgb)[3:5] <- c("sR", "sG", "sB")
      } else{
        srgb <- rgb_to_srgb_help(as.matrix(tmp[,2:4]))
        srgb <- data.frame(cbind(tmp[,1], srgb))
        colnames(srgb)[1] <- "id"
        colnames(srgb)[2:4] <- c("sR", "sG", "sB")
      }
    } else{
      stop("Cannot obtain the RGB for each object since `object_index` argument was not used. \nHave you accidentally missed the argument `pixel_level_index = TRUE`?")
    }
  }
  if (any(class(object) == "Image")){
    srgb <- rgb_to_srgb_help(cbind(c(object[,,1]), c(object[,,2]), c(object[,,3])))
    colnames(srgb) <- c("sR", "sG", "sB")
  }
  return(data.frame(srgb))
}


#' @export
#' @name utils_colorspace
rgb_to_lab <- function(object){
  object <- rgb_to_srgb(object)
  srgb <- data.frame(r = object[, 1],
                     g = object[, 2],
                     b = object[, 3])
  lab <- convertColor(srgb, from = "sRGB", to = "Lab")
  return(lab)
}




# Faster alternatives (makes only the needed)
help_segment <- function(img,
                         index = NULL,
                         threshold = c("Otsu", "adaptive"),
                         k = 0.1,
                         windowsize = NULL,
                         col_background = NULL,
                         has_white_bg = FALSE,
                         fill_hull = FALSE,
                         filter = FALSE,
                         re = NULL,
                         nir = NULL,
                         invert = FALSE){
  img2 <- help_binary(img,
                      index = index,
                      threshold = threshold,
                      k = k,
                      windowsize = windowsize,
                      has_white_bg = has_white_bg,
                      resize = FALSE,
                      fill_hull = fill_hull,
                      filter = filter,
                      re = re,
                      nir = nir,
                      invert = invert)
  ID <- which(img2@.Data == FALSE)
  img@.Data[,,1][ID] <- 1
  img@.Data[,,2][ID] <- 1
  img@.Data[,,3][ID] <- 1
  invisible(img)
}



help_binary <- function(img,
                        index = NULL,
                        threshold = c("Otsu", "adaptive"),
                        k = 0.1,
                        windowsize = NULL,
                        has_white_bg = FALSE,
                        resize = FALSE,
                        fill_hull = FALSE,
                        filter = FALSE,
                        re = NULL,
                        nir = NULL,
                        invert = FALSE){
  threshold <- threshold[[1]]

  bin_img <- function(imgs,
                      invert,
                      fill_hull,
                      threshold,
                      filter){
    # adapted from imagerExtra  https://bit.ly/3Wp4pwv
    if(threshold == "adaptive"){
      if(is.null(windowsize)){
        windowsize <- min(dim(imgs)) / 3
        if(windowsize %% 2 == 0){
          windowsize <- as.integer(windowsize + 1)
        }
      }
      if (windowsize <= 2){
        stop("windowsize must be greater than or equal to 3", call. = FALSE)
      }
      if (windowsize %% 2 == 0){
        warning(sprintf("windowsize is even (%d). windowsize will be treated as %d", windowsize, windowsize + 1), call. = FALSE)
        windowsize <- as.integer(windowsize + 1)
      }
      if (windowsize >= dim(imgs)[[1]] || windowsize >= dim(imgs)[[2]]){
        warning("windowsize is too large. Setting to `min(dim(img)) / 3`", call. = FALSE)
        windowsize <- min(dim(imgs)) / 3
      }
      if (k > 1){
        stop("k is out of range. k must be in [0, 1].", call. = FALSE)
      }
      imgs <- EBImage::Image(threshold_adaptive(as.matrix(imgs), k, windowsize, 0.5))
    }
    if(threshold != "adaptive"){
      if(threshold == "Otsu"){
        if(any(is.infinite(imgs)) | any(is.na(imgs))){
          threshold <- help_otsu(imgs@.Data[!is.infinite(imgs@.Data) & !is.na(imgs@.Data)])
        } else{
          threshold <- help_otsu(imgs@.Data)
        }
      } else{
        if(is.numeric(threshold)){
          threshold <- threshold
        } else{
          pixels <- raster::raster(t(imgs@.Data))
          raster::plot(pixels, col = custom_palette(),  axes = FALSE, asp = NA)
          threshold <- readline("Selected threshold: ")
        }
      }
      imgs <- EBImage::Image(imgs < threshold)
    }

    if(invert == TRUE){
      imgs <- 1 - imgs
    }

    imgs[which(is.na(imgs))] <- FALSE
    if(isTRUE(fill_hull)){
      imgs <- EBImage::fillHull(imgs)
    }
    if(is.numeric(filter) & filter > 1){
      imgs <- EBImage::medianFilter(imgs, filter)
    }
    return(imgs)
  }

  gray_img <- help_imageindex(img, index, resize, re, nir, has_white_bg)
  bin_img <- bin_img(gray_img,
                     invert,
                     fill_hull,
                     threshold,
                     filter)
  invisible(bin_img)
}



help_imageindex <- function(img,
                            index = NULL,
                            resize = FALSE,
                            re = NULL,
                            nir = NULL,
                            has_white_bg = FALSE){
  if(resize != FALSE){
    img <- image_resize(img, resize)
  }
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")

  nir_ind <- as.character(ind$Index[ind$Band %in% c("RedEdge","NIR")])
  hsb_ind <- as.character(ind$Index[ind$Band == "HSB"])

  R <- try(img@.Data[,,1], TRUE)
  G <- try(img@.Data[,,2], TRUE)
  B <- try(img@.Data[,,3], TRUE)

  if(any(index %in% hsb_ind)){
    hsb <- rgb_to_hsb(data.frame(R = c(R), G = c(G), B = c(B)))
    h <- matrix(hsb$h, nrow = nrow(img), ncol = ncol(img))
    s <- matrix(hsb$s, nrow = nrow(img), ncol = ncol(img))
    b <- matrix(hsb$b, nrow = nrow(img), ncol = ncol(img))
  }

  if(!is.null(re)|!is.null(nir)){
    if(index %in% nir_ind & is.null(nir)){
      stop(paste("Index ", index, " need NIR/RedEdge band to be calculated."), call. = FALSE)
    }
    if(!is.null(re)){
      RE <-  try(img@.Data[,,re], TRUE)
    }
    if(!is.null(nir)){
      NIR <- try(img@.Data[,,nir], TRUE)
    }
    test_nir_ne <- any(lapply(list(RE, NIR), class)  == "try-error" )
    if(isTRUE(test_nir_ne)){
      stop("RE and/or NIR is/are not available in your image.", call. = FALSE)
    }
  }
  if(isTRUE(has_white_bg)){
    R[which(R == 1 & G == 1 & B == 1)] <- NA
    G[which(R == 1 & G == 1 & B == 1)] <- NA
    B[which(R == 1 & G == 1 & B == 1)] <- NA
  }

  if(index %in% ind$Index){
    img_gray <- EBImage::Image(eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index]))))
  } else{
    img_gray <- EBImage::Image(eval(parse(text = as.character(index))))
  }
  invisible(img_gray)
}


#' Create an `Image` object
#'
#' This function is a simple wrapper around [EBImage::Image()].
#'
#' @param data A vector or array containing the pixel intensities of an image.
#'   If missing, the default 1x1 zero-filled array is used.
#' @param ... Additional arguments passed to [EBImage::Image()].
#' @return An `Image` object.
#' @export
#'
#' @examples
#' img <-
#' as_image(rnorm(150 * 150 * 3),
#'          dim = c(150, 150, 3),
#'          colormode = 'Color')
#' plot(img)
as_image <- function(data, ...){
  EBImage::Image(data, ...)
}


#' Prepare images to analyze_objects_shp()
#'
#' It is a simple wrapper around [image_align()] and [image_crop()]. In this case, only the option `viewer = "base"` is used. To use `viewer = "mapview"`, please, use such functions separately.
#'
#' @param img A `Image` object
#' @inheritParams image_align
#'
#' @return An aligned and cropped `Image` object.
#' @export
#'
#' @examples
#' if(interactive()){
#' img <- image_pliman("flax_leaves.jpg")
#' prepare_to_shp(img)
#' }
prepare_to_shp <- function(img,
                           align = "vertical"){
  aligned <- image_align(img, viewer = "base")
  cropped <- image_crop(aligned, viewer = "base", plot = TRUE)
  return(cropped)
}

