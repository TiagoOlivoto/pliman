#'Combines images to a grid
#'
#'Combines several images to a grid
#' @param ... a comma-separated name of image objects or a list containing image
#'   objects.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @import ggplot2
#' @importFrom stats reshape
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A grid with the images in `...`
#' @examples
#' library(pliman)
#'img1 <- image_import(image_pliman("sev_leaf.jpg"))
#'img2 <- image_import(image_pliman("sev_leaf_nb.jpg"))
#'image_combine(img1, img2)
image_combine <- function(..., nrow = NULL, ncol = NULL){
  if(is.list(c(...))){
    plots <- as.list(...)
  }else{
    plots <- list(...)
  }
  num_plots <- length(plots)
  if (is.null(nrow) && is.null(ncol)){
    ncol <- ceiling(sqrt(num_plots))
    nrow <- ceiling(num_plots/ncol)
  }
  if (is.null(ncol)){
    nrow <- ceiling(num_plots/ncol)
  }
  if (is.null(nrow)){
    nrow <- ceiling(num_plots/ncol)
  }
  op <- par(mfrow = c(nrow, ncol))
  on.exit(par(op))
  for(i in 1:length(plots)){
    plot(plots[[i]])
  }
}
#'Import, display and export images
#'
#'Import images from files and URLs, write images to files, and show images.
#'*
#' @name utils_image
#' @param image
#' * For `image_import()`, a character vector of file names or URLs.
#' * For `image_export()`, an Image object, an array or a list of images.
#' @param name An string specifying the name of the image.
#' @param img_pattern A pattern of file name used to identify images to be
#'   imported. For example, if `img_pattern = "im"` all images in the current
#'   working directory that the name matches the pattern (e.g., img1.-,
#'   image1.-, im2.-) will be imported as a list. Providing any number as
#'   pattern (e.g., `img_pattern = "1"`) will select images that are named as
#'   1.-, 2.-, and so on.
#' @param ... Alternative arguments passed to the corresponding functions from
#'   the `jpeg`, `png`, and `tiff` packages.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' * `image_import()` returns a new `Image` object.
#' * `image_export()` returns an invisible vector of file names.
#' * `image_pliman()` returns a character string with the path to the example
#' image installed with the package.
#' @examples
#' library(pliman)
#'img <- image_import(image_pliman("sev_leaf.jpg"))
image_import <- function(image, ..., img_pattern = NULL){
  if(!is.null(img_pattern)){
    if(img_pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      img_pattern <- "^[0-9].*$"
    }
    imgs <- list.files(pattern = img_pattern, getwd())
    extensions <- as.character(sapply(imgs, file_extension))
    names_plant <- as.character(sapply(imgs, file_name))
    if(length(grep(img_pattern, imgs)) == 0){
      stop(paste("'", img_pattern, "' pattern not found in '",
                 paste0(getwd())),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    list_img <-
      lapply(imgs, function(x){
        readImage(x, ...)
      })
    names(list_img) <- imgs
    return(list_img)
  }
  img_dir <- file_dir(image)
  all_files <- sapply(list.files(img_dir), file_name)
  img_name <- file_name(image)
  if(!grepl("http", img_dir, fixed = TRUE) & !img_name %in% all_files){
    stop(" '", img_name, "' not found in ", img_dir,  call. = FALSE)
  }
  readImage(image, ...)
}

#' @export
#' @name utils_image
image_export <- function(image, name, ...){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    lapply(seq_along(image), function(i){
      writeImage(x = image[[i]], files = names(image[i]), ...)
    })
  }
  writeImage(image, name, ...)
}
#' @export
#' @name utils_image
image_show <- function(image){
  if(any(class(image) != "Image")){
    grid.raster(image)
  } else{
    plot(image)
  }
}
#' @export
#' @name utils_image
image_pliman <- function(image){
  files <- list.files(system.file("tmp_images", package = "pliman"))
  if(!image %in% files){
    stop("Image not available in pliman.\nAvaliable images: ", paste(files, collapse = ", "), call. = FALSE)
  }
  system.file(paste0("tmp_images/", image), package = "pliman")
}

##### Spatial transformations
#'Spatial transformations
#'
#' Performs image rotation and reflection
#' * `image_dimension()` Gives the dimension (width and height) of an image.
#' * `image_rotate()` rotates the image clockwise by the given angle.
#' * `image_horizontal()` converts (if needed) an image to a horizontal image.
#' * `image_vertical()` converts (if needed) an image to a vertical image.
#' * `image_hreflect()` performs horizontal reflection of the `image`.
#' * `image_vreflect()` performs vertical reflection of the `image`.
#' * `image_resize()` resize the `image`.
#' @name utils_transform
#' @param image An image or a list of images of class `Image`.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param angle The rotation angle in degrees.
#' @param bg_col Color used to fill the background pixels, defaults to `"white"`.
#' @param rel_size The relative size of the resized image. Defaults to 100. For
#'   example, setting `rel_size = 50` to an image of width `1280 x 720`, the new
#'   image will have a size of `640 x 360`.
#' @param width,height Width and height of the resized image. These arguments
#'   can be missing. In this case, the image is resized according to the
#'   relative size informed in `rel_size`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @md
#' @importFrom parallel detectCores clusterExport makeCluster parLapply stopCluster
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A modified version of `image` depending on the function used.
#' @examples
#' library(pliman)
#'img <- image_import(image_pliman("sev_leaf.jpg"))
#'image_show(img)
#'img <- image_resize(img, 50)
#'img1 <- image_rotate(img, 45)
#'img2 <- image_hreflect(img)
#'img3 <- image_vreflect(img)
#'img4 <- image_vertical(img)
#'image_combine(img1, img2, img3, img4)
image_dimension <- function(image,
                            parallel = FALSE,
                            workers = NULL,
                            verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
      as.data.frame(
      do.call(rbind,
              lapply(image, image_dimension, verbose =  FALSE))
      )
      res <- transform(res, image = rownames(res))[,c(3, 1, 2)]
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
                         verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
    rotate(image, angle, bg.col = bg_col)
  }
}
#' @name utils_transform
#' @export
image_horizontal <- function(image,
                             parallel = FALSE,
                             workers = NULL,
                             verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
      image_rotate(image, 90)
    } else{
      image
    }
  }
}
#' @name utils_transform
#' @export
image_vertical <- function(image,
                           parallel = FALSE,
                           workers = NULL,
                           verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
      image_rotate(image, 90)
    } else{
      image
    }
  }
}
#' @name utils_transform
#' @export
image_hreflect <- function(image,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
    flop(image)
  }
}
#' @name utils_transform
#' @export
image_vreflect <- function(image,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
    flip(image)
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
                         verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
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
    resize(image, width, height)
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
#'   target mode for conversion to binary image. One of the following:  `"R"`,
#'   `"G"`, `"B"` `"GR"`, `"NR"`, `"NG"`, `"NB"`, `"BI"`, `"BIM"`, `"SCI"`,
#'   `"GLI"`, `"HI"`, `"NGRDI"`, `"SI"`, `"VARI"`, `"HUE"`, `"HUE2"`, `"BGI"`,
#'   `"BGI"`. See [image_index()] for more details.
#' @param my_index User can calculate a different index using the bands names,
#'   e.g. `my_index = "R+B/G"`.
#' @param resize Resize the image before processing? Defaults to `TRUE`. Resizes
#'   the image to 30% of the original size to speed up image processing.
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
#'   up to 90% of available cores.
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
#'img <- image_import(image_pliman("soybean_touch.jpg"))
#'image_binary(img, index = c("R, NR"), nrow = 1)
image_binary <- function(image,
                         index = NULL,
                         my_index = NULL,
                         resize = TRUE,
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
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_binary, index, my_index, resize, fill_hull, re, nir, invert, show_image, nrow, ncol)
    } else{
    lapply(image, image_binary, index, my_index, resize, fill_hull, re, nir, invert, show_image, nrow, ncol)
    }
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
    img2 <- image_index(image, indx, my_index, resize, re, nir, show_image = FALSE, nrow, ncol)[[1]]
    threshold <- otsu(img2, range = range(img2[!is.infinite(img2)], na.rm = TRUE))
    if(invert == FALSE){
      img2 <- combine(mapply(function(frame, th) frame < th, getFrames(img2), threshold, SIMPLIFY=FALSE))
    } else{
      img2 <- combine(mapply(function(frame, th) frame > th, getFrames(img2), threshold, SIMPLIFY=FALSE))
    }
    ifelse(fill_hull == TRUE,
           imgs[[i]] <- fillHull(img2),
           imgs[[i]] <- img2)
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
#' Builds image indexes using Red, Green, Blue, Red-Edge, and NIR bands.
#' @details
#' The following indexes are available in pliman.
#'
#'  | Index | Equation                   | Band |
#'  |-------|----------------------------|------|
#'  | R     | R                          | C    |
#'  | G     | G                          | C    |
#'  | B     | B                          | C    |
#'  | NR    | R/(R+G+B)                  | C    |
#'  | NG    | G/(R+G+B)                  | C    |
#'  | NB    | B/(R+G+B)                  | C    |
#'  | BI    | sqrt((R^2+G^2+B^2)/3)      | C    |
#'  | BIM   | sqrt((R*2+G*2+B*2)/3)      | C    |
#'  | SCI   | (R-G)/(R+G)                | C    |
#'  | GLI   | (2*G-R-B)/(2*G+R+B)        | C    |
#'  | HI    | (2*R-G-B)/(G-B)            | C    |
#'  | NGRDI | (G-R)/(G+R)                | C    |
#'  | SI    | (R-B)/(R+B)                | C    |
#'  | VARI  | (G-R)/(G+R-B)              | C    |
#'  | HUE   | atan(2*(B-G-R)/30.5*(G-R)) | C    |
#'  | HUE2  | atan(2*(R-G-R)/30.5*(G-B)) | C    |
#'  | BGI   | B/G                        | C    |
#' @name image_index
#' @param image An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. One of the following:  `"R"`,
#'   `"G"`, `"B"` `"GR"`, `"NR"`, `"NG"`, `"NB"`, `"BI"`, `"BIM"`, `"SCI"`,
#'   `"GLI"`, `"HI"`, `"NGRDI"`, `"SI"`, `"VARI"`, `"HUE"`, `"HUE2"`, `"BGI"`,
#'   `"BGI"`. Defaults to `NULL` ((normalized) Red, Green and Blue).  One can
#'   also use "RGB" for RGB only, "NRGB" for normalized RGB, or "all" for all
#'   indexes.
#' @param my_index User can calculate a different index using the bands names,
#'   e.g. `my_index = "R+B/G"`.
#' @param resize Resize the image before processing? Defaults to `TRUE`. Resizes
#'   the image to 30% of the original size to speed up image processing.
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
#'   up to 90% of available cores.
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
#'img <- image_import(image_pliman("soybean_touch.jpg"))
#'image_index(img, c("R, NR"), nrow = 1)
image_index <- function(image,
                        index = NULL,
                        my_index = NULL,
                        resize = TRUE,
                        re = NULL,
                        nir = NULL,
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
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_index, index, my_index, resize, re, nir, show_image, nrow, ncol)
    } else{
      lapply(image, image_index, index, my_index, resize, re, nir, show_image, nrow, ncol)
    }
  } else{
  if(resize == TRUE){
    image <- image_resize(image, 30)
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
                G = "blue",
                B = "green",
                GR = "gray"
        )
      img2 <- channel(image, indx)
    } else{
      if(is.null(my_index) & !indx %in% ind$Index){
        stop(paste("Index '",indx,"' is not available in pliman",sep = ""), call. = FALSE)
      }
      frames <- getFrames(image)
      num_band <- length(frames)
      if(num_band < 3){
        stop("At least 3 bands (RGB) are necessary to calculate indices available in pliman.", call. = FALSE)
      }
      R <- frames[[1]]
      G <- frames[[2]]
      B <- frames[[3]]
      if(indx %in% nir_ind & is.null(nir)){
        stop(paste("Index ", indx, " need NIR/RedEdge band to be calculated."), call. = FALSE)
      }
      if(!is.null(re)|!is.null(nir)){
        if(num_band < 4){
          stop("RE and/or NIR is/are not available in your image.", call. = FALSE)
        }
        if(!is.null(re)){
          RE <- frames[[re]]
        }
        if(!is.null(re)){
          NIR <- frames[[nir]]
        }
      }
      if(is.null(my_index)){
      img2 <- eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==indx])))
      } else{
      img2 <- eval(parse(text = as.character(my_index)))
      }
    }

    imgs[[i]] <- img2
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
#' Produces an histogram of an `image_index` object
#'
#' @name image_index
#' @param x An object of class `image_index`.
#' @param facet Shows RGB values as a facet plot? Defaults to `TRUE`.
#' @param ... Currently not used
#' @method plot image_index
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A `ggplot` object containing the distribution of the pixels for each
#'   index.
#' @examples
#' library(pliman)
#'img <- image_import(image_pliman("sev_leaf.jpg"))
# A half size of the original image
#'img2 <- image_resize(img, 50)
#'ind <- image_index(img2)
#'plot(ind)
plot.image_index <- function(x, facet = TRUE, ...){
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
  ggplot(a, aes(value, fill = Spectrum)) +
    geom_density(alpha = 0.6) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    scale_x_continuous(expand = expansion(c(0, 0))) +
    {if(facet)facet_wrap(~Spectrum, ncol = 3, scales = "free")} +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.ticks.length = unit(0.2, "cm"),
          panel.grid.minor = element_blank())
}






#' Image segmentation
#'
#' Reduce a color, color near-infrared, or grayscale images to a segmented image
#' using a given color channel (red, green blue) or even color indexes (See
#' [image_index()] for more details). The Otsu's thresholding method (Otsu,
#' 1979) is used to automatically perform clustering-based image thresholding.
#'
#' @param image An image object.
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to binary image. One of the following:  `"R"`,
#'   `"G"`, `"B"` `"GR"`, `"NR"`, `"NG"`, `"NB"`, `"BI"`, `"BIM"`, `"SCI"`,
#'   `"GLI"`, `"HI"`, `"NGRDI"`, `"SI"`, `"VARI"`, `"HUE"`, `"HUE2"`, `"BGI"`,
#'   `"BGI"`. See [image_index()] for more details.
#' @param my_index User can calculate a different index using the bands names,
#'   e.g. `my_index = "R+B/G"`.
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
#'   up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return A list containing `n` objects where `n` is the number of indexes
#'   used. Each objects contains:
#' * `image` an image with the RGB bands (layers) for the segmented object.
#' * `mask` A mask with logical values of 0 and 1 for the segmented image.
#' @examples
#' library(pliman)
#'img <- image_import(image_pliman("soybean_touch.jpg"))
#'image_show(img)
#'image_segment(img, index = c("R, G, B"))
image_segment <- function(image,
                         index = NULL,
                         my_index = NULL,
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
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_segment, index, my_index, fill_hull, re, nir, invert, show_image, nrow, ncol)
    } else{
      lapply(image, image_segment, index, my_index, fill_hull, re, nir, invert, show_image, nrow, ncol)
    }
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
      img2 <- image_index(image, indx, my_index, resize = FALSE, re, nir, show_image = FALSE, nrow, ncol)[[1]]
      threshold <- otsu(img2, range = range(img2[!is.infinite(img2)], na.rm = TRUE))
      if(invert == FALSE){
        img2 <- combine(mapply(function(frame, th) frame < th, getFrames(img2), threshold, SIMPLIFY=FALSE))
      } else{
        img2 <- combine(mapply(function(frame, th) frame > th, getFrames(img2), threshold, SIMPLIFY=FALSE))
      }
      ifelse(fill_hull == TRUE,
             img2 <- fillHull(img2),
             img2 <- img2)
      ID <- which(img2@.Data == FALSE)
      img <- image
      img@.Data[,,1][ID] <- 1
      img@.Data[,,2][ID] <- 1
      img@.Data[,,3][ID] <- 1
      mask <- image@.Data[,,1]
      mask[ID] <- 1
      mask[!ID] <- 0
      mask <- as.Image(mask)
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
    invisible(imgs)
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
#'   up to 90% of available cores.
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
#'img <- image_import(image_pliman("sev_leaf.jpg"))
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
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, "image")
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_to_mat)
    } else{
      lapply(image, image_to_mat)
    }
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


#' Create an image palette
#'
#' Creates image palettes by applying the k-means algorithm to the RGB values.
#' @param image An image object.
#' @param npal The number of color palettes.
#' @param nstart How many random sets from `npal` should be chosen?
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return A list with `npal` color palettes of class `Image`.
#' @export
#' @examples
#' \donttest{
#' library(pliman)
#'img <- image_import(image_pliman("sev_leaf.jpg"))
#'pal <- image_palette(img, 2)
#'image_show(pal[[1]])
#'image_show(pal[[2]])
#'}
image_palette <- function(image,
                          npal,
                          nstart = 25,
                          parallel = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, c("image", "image_to_mat"))
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      parLapply(clust, image, image_palette, npal, nstart)
    } else{
      lapply(image, image_palette, npal, nstart)
    }
  } else{
  df <- image_to_mat(image)$df_in
  df$CODE <- NULL
  rgbs <- list()
  b <- kmeans(df, centers = npal, nstart = nstart)
  for (i in 1:npal) {
    df1 <- df[which(b$cluster == i),]
    dim_mat <- trunc(sqrt(nrow(df1)))
    R <- df1[1:dim_mat^2, 1]
    G <- df1[1:dim_mat^2, 2]
    B <- df1[1:dim_mat^2, 3]
    rgbs[[i]] <- as.Image(array(c(R, G, B), dim = c(dim_mat, dim_mat, 3)))
  }
  return(rgbs)
  }
}


#' Utilities for image resolution
#'
#' Provides useful conversions between size (cm), number of pixels (px) and
#' resolution (dpi)
#'
#' @name utils_dpi
#' @param dpi The image resolution in dots per inch.
#' @param px The number of pixels.
#' @param cm The size in centimeters.
#' @return A numeric value or a vector of numeric values if the input data is a
#'   vector
#' @export
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
