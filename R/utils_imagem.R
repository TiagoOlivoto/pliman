#'Import, export and display images
#'
#'Basic tools to work with images .png, .jpeg, .jpg, .tiff (ou .tif, .TIFF,
#'.TIF)
#'* `import_image()` Imports an image.
#'* `show_image()` Show an image.
#'* `save_image()` Saves an image.
#'* `image_to_mat()` Convert an image to a list, containing three matrices (RGB)
#'and a data frame with three columns, one for each of the RGB Spectrum. (RGB).
#'* `image_rgb()` Produces a density plot with the RGB values of an image.
#'* `image_pallete()` Generate desired number of color palettes from an image
#'using the k-means algorithm.
#' @name utils_image
#' @param image An image.
#'* For `show_image()`, any R object that can be coerced to a raster object.
#'* For `import_image()`, a string specifying the path of an image.
#'* For `save_image()`, a two-dimensional or three-dimensional object (matrix,
#'data.frame or array).
#' @param name An string specifying the name of the image
#' @param randomize Randomize as lines of array conversation for data
#'frame? Padrao `TRUE`.
#' @param nrows Number of lines to be selected, no data frame created.
#' @param facet Shows RGB values as a facet plot? Defaults to `FALSE`.
#' @param npal The number of color palletes.
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param nstart How many random sets from `npal` should be chosen?
#' @param ... Depends on the function:
#'    * For `import_image()` and `save_image()`, alternative arguments passed to
#' the corresponding functions from the `jpeg`, `png`, and `tiff` packages.
#'    * For `combine_images()`, a comma-separated list of image objects.
#' @md
#' @importFrom EBImage readImage writeImage
#' @import ggplot2
#' @importFrom stats reshape
#' @export
#' @examples
#' library(pliman)
#'img <- import_image(system.file("tmp_images", "sev3.png", package = "pliman"))
#'image_rgb(img)
show_image <- function(image){
  if(any(class(image) != "Image")){
    grid.raster(image)
  } else{
    plot(image)
  }
}
#' @export
#' @name utils_image
combine_images <- function(..., nrow = NULL, ncol = NULL){
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
#' @export
#' @name utils_image
import_image <- function(image, ...){
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
save_image <- function(image, name, ...){
  writeImage(image, name, ...)
}
#' @export
#' @name utils_image
image_to_mat <- function(image, randomize =  TRUE, nrows = 5000){
  d <- match.call()
  ncols <- ncol(image@.Data[,,1])
  im <- cbind(c(image@.Data[,,1]), c(image@.Data[,,2]), c(image@.Data[,,3]))
  df_in <- data.frame(im) %>% transform(code = paste(d[["image"]]))
  df_in <- df_in[c(4, 1, 2, 3)]
  df_man <- df_in
  if(randomize == TRUE){
    df_man <- df_man[sample(1:nrow(df_man)),]
  }
  if(!missing(nrows)){
    df_man <- df_man[1:nrows, ]
  }
  colnames(df_man) <- colnames(df_in) <-  c("CODE", "R", "G", "B")
  rbg <- list(R = matrix(im[, 1], ncol = ncols),
              G = matrix(im[, 2], ncol = ncols),
              B = matrix(im[, 3], ncol = ncols),
              df_man = data.frame(df_man),
              df_in = df_in)
  return(rbg)
}
#' @export
#' @name utils_image
image_rgb <- function(image, facet = FALSE){
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
#' @export
#' @name utils_image
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
correct_image =function(image, perc){
  t <- image
  n=round(perc*min(c(ncol(t),nrow(t))),0)
  p1=function(t){
    t2=t
    for( i in 2:(nrow(t)-n-1)){
      for(j in 1:ncol(t)){
        if(t[i,j]==1){
          if(t[i-1,j]==0){
            a=0
            while(a<n){
              a=a+1

              if(sum(t[i:(i+a),j]==1)<a){t2[i:(i+a),j]=0;a=n}
            }
          }
        }
      }
    }
    return(t2)
  }
  Pp=p1(t)
  Pp=p1(t(Pp))
  return(t(Pp))
}
