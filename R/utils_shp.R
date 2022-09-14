#' Construct a shape file from an image
#'
#' Creates a list of object coordinates given the desired number of rows and
#' columns. It starts by selecting 4 points at the corners of objects of
#' interest in the plot space. Then, given `rows` and `cols`, a grid is drawn
#' and the objects' coordinates are returned.
#'
#' @param img An object of class `Image`
#' @param rows The number of desired rows in the grid. Defaults to `1`.
#' @param cols The number of desired columns in the grid. Defaults to `1`.
#' @param col_line,col_text The color of the line/text in the grid. Defaults to
#'   `"blue"`.
#' @param size_line,size_text The size of the line/text in the grid. Defaults to
#'   `2.5`.
#'
#' @return A list with `row * col` objects containing the plot coordinates.
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' shape <- image_shp(flax, rows = 3, cols = 5)
#' }
#'
image_shp <- function(img,
                      rows = 1,
                      cols = 1,
                      col_line = "blue",
                      size_line = 2,
                      col_text = "blue",
                      size_text = 1){
  message("Select 2 points drawing the diagonal that includes the objects of interest.")
  plot(img)
  cord <- locator(type = "p", n = 2, col = "red", pch = 19)
  c1 <- data.frame(do.call(rbind, cord)) |> t()
  xmin <- min(c1[,1])
  xmax <- max(c1[,1])
  ymin <- min(c1[,2])
  ymax <- max(c1[,2])
  bbox <-
    data.frame(x = c(xmin, xmax, xmax, xmin, xmin),
               y = c(ymin, ymin, ymax, ymax, ymin))
  xr <- xmax - xmin
  yr <- ymax - ymin
  intx <- xr / (cols )
  xvec <- xmin
  for (i in 1:(cols)) {
    xvec <- append(xvec, xvec[length(xvec)] + intx)
  }
  inty <- yr / (rows )
  yvec <- ymin
  for (i in 1:(rows)) {
    yvec <- append(yvec, yvec[length(yvec)] + inty)
  }
  coords <- list()
  con <- 0
  plot(img)
  for(i in 1:rows){
    for(j in 1:cols){
      con <- con + 1
      tmp <-
        data.frame(plot = con,
                   x = c(xvec[j], xvec[j + 1], xvec[j + 1], xvec[j], xvec[j]),
                   y = c(yvec[i], yvec[i], yvec[i + 1], yvec[i + 1], yvec[i]))
      lines(tmp[, -1], col = col_line, type = "l", lwd = size_line)
      text(min(tmp$x[-1]), min(tmp$y[-1]),
           label = con,
           col = col_text,
           cex = size_text,
           adj = c(-0.2, 1.2))
      coords[[paste0("plot_", con)]] <- tmp
    }
  }
  lst <- list(shapefiles = coords,
              bbox = bbox,
              rows = rows,
              cols = cols)
  return(structure(lst, class = "image_shp"))
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
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg")
#' shape <- image_shp(flax, rows = 3, cols = 5)
#'
#' # grid on the existing image
#' plot(flax)
#' plot(shape)
#' }
plot.image_shp <- function(x,
                           img = NULL,
                           col_line = "blue",
                           size_line = 2,
                           col_text = "blue",
                           size_text = 1,
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
#' @param ... Other arguments passed on to [image_shp()]
#' @return A list of `Image` objects
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' objects <- object_split_shp(flax, rows = 3, cols = 5)
#' image_combine(objects)
#' }
object_split_shp <- function(img,
                             rows = 1,
                             cols = 1,
                             ...){
  shapefile <- image_shp(img, rows, cols, ...)$shapefiles
  imgs <- list()
  get_borders <- function(x){
    min_x <- min(x[,1])
    max_x <- max(x[,1])
    min_y <- min(x[,2])
    max_y <- max(x[,2])
    return(list(min_x, max_x, min_y, max_y))
  }
  for (i in 1:length(shapefile)) {
    tmp <- shapefile[[i]][-1]
    borders <- get_borders(tmp)
    imgs[[paste0("shp", i)]] <-
      image_crop(img,
                 width = borders[[1]]:borders[[2]],
                 height = borders[[3]]:borders[[4]])
  }
  return(imgs)
}



#' Aligns an `Image` object by hand
#'
#' [image_rotate()] rotate an image given an informed angle. Here, users will
#' need to draw a line along the y axis that corresponds to the alignment of the
#' objects (e.g., field plots). By default, the aligment will be to the
#' vertical, which means that if the drawed line have an angle < 90º parallel to
#' the x axis, the rotation angle wil be negative (anticlocwise rotation). If
#' the drawed line have an angle > 90º along the x axis, the rotation angle wil
#' be positive (clocwise rotation). If the aligment is horizontal, the image
#' will be rotated to align the drawed line paralell to the x axis assuming the
#' shortest angle (left or right rotation if the angle is < 90º or > 90º,
#' respectively)
#'
#' @param img An `Image` object
#' @param align The desired alignment. Either `"vertical"` (default) or
#'   `"horizontal"`.
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
                        plot = TRUE){
  message("Select 2 points drawing a line of desired aligment along the y axis.")
  plot(img)
  cord <- locator(type = "p", n = 2, col = "red", pch = 19)
  c1 <- data.frame(do.call(rbind, cord)) |> t()
  lines(c1, col = "red", lty = 2, lwd = 2)
  a <- abs(cord$x[1] - cord$x[2])
  b <- abs(cord$y[1] - cord$y[2])
  angle <- (atan(b / a) * 180) / pi
  if(svd(var(c1))$u[2] >= 0){
    anglev <- angle - 90
    angleh <- angle
    if(align[[1]] == "vertical"){
      img2 <- image_rotate(img, angle = anglev, plot = plot)
    } else{
      img2 <- image_rotate(img, angle = angleh, plot = plot)
    }
    message(paste("Angle to align in the vertical: ", round(anglev, 3)))
    message(paste("Angle to align in the horizontal: ", round(angleh, 3)))
  } else{
    anglev <- 90 - angle
    angleh <- angle * -1
    if(align[[1]] == "vertical"){
      img2 <- image_rotate(img, angle = anglev, plot = plot)
    } else{
      img2 <- image_rotate(img, angle = angleh, plot = plot)
    }
    message(paste("Angle to align in the vertical: ", round(anglev, 3)))
    message(paste("Angle to align in the horizontal: ", round(angleh, 3)))
  }
  return(img2)
}


#' Analyzes objects using shapefiles
#'
#' This function calls [analyze_objects()] in each image polygon of a shapefile
#' object generated with generated with [image_shp()] and bind the results into
#' read-ready data frames.
#'
#' @param img An `Image` object
#' @param shapefile (Optional) An object created with [image_shp()]. If `NULL`
#'   (default), both `rows` and `cols` must be declared.
#' @param rows,cols The number of rows and columns to generate the shapefile
#'   when `shapefile` is not declared.
#' @param show_image Shows the processed images? Defaults to `TRUE`.
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
#' flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
#' res <- analyze_objects_shp(flax, rows = 3, cols = 5)
#' }
analyze_objects_shp <- function(img,
                                shapefile = NULL,
                                rows = NULL,
                                cols = NULL,
                                show_image = TRUE,
                                object_size = "large",
                                ...){
  if(is.null(shapefile) & is.null(rows)){
    stop("When 'shapefile' is not informed, 'rows' and 'cols must be declared", call. = FALSE)
  }
  if(is.null(shapefile)){
    spl <- object_split_shp(img, rows, cols)
  } else{
    shapes <- shapefile$shapefiles
    spl <- list()
    rows <- shapefile$rows
    cols <- shapefile$cols
    get_borders <- function(x){
      min_x <- min(x[,1])
      max_x <- max(x[,1])
      min_y <- min(x[,2])
      max_y <- max(x[,2])
      return(list(min_x, max_x, min_y, max_y))
    }
    for (i in 1:length(shapes)) {
      tmp <- shapes[[i]][-1]
      borders <- get_borders(tmp)
      spl[[paste0("obj", i)]] <-
        image_crop(img,
                   width = borders[[1]]:borders[[2]],
                   height = borders[[3]]:borders[[4]])
    }
  }
  if(isTRUE(show_image)){
    op <- par(mfrow = c(rows, cols))
    on.exit(par(op))
  }
  list <- lapply(spl,
                 analyze_objects,
                 img,
                 show_image = show_image,
                 object_size = object_size,
                 ...)
  results <-
    do.call(rbind,
            lapply(list, function(x){x$results}))
  vect <- rownames(results)
  results$img <-
    sapply(seq_along(vect),
           function(i){
             strsplit(vect[[i]], split = "\\.")[[1]][[1]]
           })
  results <- results[, c(ncol(results), 1:(ncol(results) - 1))]
  rownames(results) <- NULL

  statistics <-
    do.call(rbind,
            lapply(seq_along(list), function(x){
              transform(list[[x]][["statistics"]], img = names(list[x]))[,c(3, 1, 2)]
            }))

  if(!is.null(list[[1]][["object_rgb"]])){
    object_rgb <-
      do.call(rbind,
              lapply(seq_along(list), function(i){
                transform(list[[i]][["object_rgb"]], img = names(list[i]))
              }))
    object_rgb <- object_rgb[, c(ncol(object_rgb), 1:(ncol(object_rgb) - 1))]
  } else{
    object_rgb <- NULL
  }

  if(!is.null(list[[1]][["object_index"]])){
    object_index <-
      do.call(rbind,
              lapply(seq_along(list), function(i){
                transform(list[[i]][["object_index"]], img = names(list[i]))
              }))
    object_index <- object_index[, c(ncol(object_index), 1:(ncol(object_index) - 1))]
  } else{
    object_index <- NULL
  }

  return(structure(
    list(results = results,
         statistics = statistics,
         object_rgb = object_rgb,
         object_index = object_index),
  class = "anal_obj"))

}

# #
# library(pliman)
# x11()
# flax <- image_pliman("flax_leaves.jpg", plot = TRUE)
# shapefile <- image_shp(flax, 3, 5)
#
# res <- analyze_objects_shp(flax, shapefile, object_index = "DGCI", object_size = "medium")
# c <- get_measures(res)
