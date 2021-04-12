#'Combines images to a grid
#'
#'Combines several images to a grid
#' @param ... a comma-separated list of image objects.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @importFrom EBImage readImage writeImage
#' @import ggplot2
#' @importFrom stats reshape
#' @export
#' @examples
#' library(pliman)
#'img1 <- image_import(system.file("tmp_images", "sev1.png", package = "pliman"))
#'img2 <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'image_combine(img1, img2)
image_combine <- function(..., nrow = NULL, ncol = NULL){
  plots <- list(...)
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
#'Import, export and display images
#'
#'Import images from files and URLs, write images to files and show images.
#' @name utils_image
#' @param image
#' * For `image_import()`, a character vector of file names or URLs.
#' * For `image_export()`, an Image object or an array.
#' @param name An string specifying the name of the image.
#' @param ... Alternative arguments passed to the corresponding functions from
#'   the `jpeg`, `png`, and `tiff` packages.
#' @md
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
image_import <- function(image, ...){
  img_dir <- file_dir(image)
  all_files <- sapply(list.files(img_dir), file_name)
  img_name <- file_name(image)
  if(!img_name %in% all_files){
    stop(" '", img_name, "' not found in ", img_dir,  call. = FALSE)
  }
  readImage(image, ...)
}
#' @export
#' @name utils_image
image_export <- function(image, name, ...){
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
#' Creates a binary image
#'
#' Reduce a color or grayscale image to a binary image using a given color
#' channel. Otsu's thresholding method (Otsu, 1979) is used to automatically
#' perform clustering-based image thresholding.
#' @param image An image object.
#' @param channel A character value specifying the target mode for conversion to
#'   binary image. One of `gray`, `grey`, `red`, `green`, or `blue`. Defaults to
#'   `"blue"`.
#' @param invert Inverts the binary image, if desired.
#' @references Nobuyuki Otsu, "A threshold selection method from gray-level
#'   histograms". IEEE Trans. Sys., Man., Cyber. 9 (1): 62-66. 1979.
#'   \doi{10.1109/TSMC.1979.4310076}
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'image_show(image_binary(img, "red"))
#'image_show(image_binary(img, "red"))
image_binary <- function(image, channel = "blue", invert = FALSE){
  img2 <- channel(image, channel)
  threshold <- otsu(img2)
  if(invert == FALSE){
  img2 <- combine(mapply(function(frame, th) frame < th, getFrames(img2), threshold, SIMPLIFY=FALSE))
  } else{
  img2 <- combine(mapply(function(frame, th) frame > th, getFrames(img2), threshold, SIMPLIFY=FALSE))
  }
  return(img2)
}
#' Convert an image to numerical matrices
#'
#' Given an object image, converts it into three matrices (RGB) and a data frame
#' where each column corresponds to the RGB values.
#' @param image An image object.
#'frame? Padrao `TRUE`.
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'dim(img)
#'mat <- image_to_mat(img)
#'dim(mat[[1]])
image_to_mat <- function(image){
  d <- match.call()
  ncols <- ncol(image@.Data[,,1])
  im <- cbind(c(image@.Data[,,1]), c(image@.Data[,,2]), c(image@.Data[,,3]))
  df_in <-
    data.frame(im) %>%
    transform(code = paste(d[["image"]])) %>%
    .[c(4, 1, 2, 3)]
  colnames(df_in) <-  c("CODE", "R", "G", "B")
  rbg <- list(R = matrix(im[, 1], ncol = ncols),
              G = matrix(im[, 2], ncol = ncols),
              B = matrix(im[, 3], ncol = ncols),
              df_in = df_in)
  return(rbg)
}
#' Produces an RGB plot of an image
#' @param image An image object.
#' @param facet Shows RGB values as a facet plot? Defaults to `TRUE`.
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'dim(img)
# A half size of the original image
#'img2 <- image_resize(img, 50)
#'image_rgb(img2)
image_rgb <- function(image, facet = TRUE){
  mat <- image_to_mat(image)$df_in
  mat$CODE <- NULL
  mat$id <- rownames(mat)
  a <-
    reshape(data.frame(mat),
            idvar = "id",
            varying = list(1:3),
            times = c("R", "G", "B"),
            timevar = "Spectrum",
            v.names = "value",
            direction = "long")
  a$Spectrum <- factor(a$Spectrum, levels = c("R", "G", "B"))
  ggplot(a, aes(value, fill = Spectrum)) +
    geom_density(alpha = 0.6) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    scale_x_continuous(expand = expansion(c(0, 0))) +
    {if(facet)facet_wrap(~Spectrum, ncol = 3, scales = "free_y")} +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.ticks.length = unit(0.2, "cm"),
          panel.grid.minor = element_blank())
}
#'Resize an image object
#'
#'Scales the `image` to the desired dimensions.
#'
#'Image can be resized either by a percent value of the original image or by
#'explicitly  declaring the desired dimension in the final image. For example,
#'setting `rel_size = 50` to an image of width `1280 x 720`, the new image will
#'have a size of `640 x 360`. Declare one of `width` or `height` to enable
#'proportional resizing by pixel values.
#' @param image An image object.
#' @param rel_size The relative size of the resized image. Defaults to 100.
#' @param width,height Width and height of the resized image. These arguments
#'   can be missing. In this case, the image is resized according to the
#'   relative size informed in `rel_size.`.
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'dim(img)
#'dim(image_resize(img, 50))
image_resize <- function(image,
                         rel_size = 100,
                         width,
                         height){
  nrow <- dim(image)[[1]]
  new_row <- nrow * rel_size / 100
  width <- ifelse(missing(width), new_row, width)
  resize(image, width, height)
}
#' Create an image palette
#'
#' Creates image palettes by applying the k-means algorithm to the RGB values.
#' @param image An image object.
#' @param npal The number of color palletes.
#' @param nstart How many random sets from `npal` should be chosen?
#' @export
#' @examples
#' library(pliman)
#'img <- image_import(system.file("tmp_images", "sev3.png", package = "pliman"))
#'pal <- image_pallete(img, 2)
#'image_show(pal[[1]])
#'image_show(pal[[2]])
image_pallete <- function(image, npal, nstart = 25){
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
    rgbs[[i]] <- array(c(R, G, B), dim = c(dim_mat, dim_mat, 3))
  }
  return(rgbs)
}
