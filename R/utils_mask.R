#' Makes a brush
#'
#' Generates brushes of various sizes and shapes that can be used as structuring
#' elements. See [EBImage::makeBrush()].
#'
#' @param size A numeric containing the size of the brush in pixels. This should
#'   be an odd number; even numbers are rounded to the next odd one.
#' @param shape A character vector indicating the shape of the brush. Can be
#'   `"box"`, `"disc"`, `"diamond"`, `"Gaussian"` or `"line"` Defaults to
#'   `"disc"`.
#' @param ... Further arguments passed on to [EBImage::makeBrush()].
#' @return A 2D matrix of 0s and 1s containing the desired brush.
#' @export
#' @importFrom graphics image
#'
#' @examples
#'
#' make_brush(size = 51) |> image()
#' make_brush(size = 51, shape = "diamond") |> image()
make_brush <- function(size,
                       shape = "disc",
                       ...){
  if(size %% 2 == 0){
    size <- size + 1
    warning("'size' was rounded to the next odd number:", size, call. = FALSE)
  }
  EBImage::makeBrush(size, shape, ...)
}



#' Makes a mask in an image
#'
#' Make a mask using an `Image` object and a brush.
#'
#' It applies a brush to an Image, selecting the `Image` pixels that match the
#' brush values equal to 1. The position of the brush in the original image is
#' controlled by the relative positions x (`rel_pos_x`) and y (`rel_pos_y`)
#' arguments.  The size of the brush must be smaller or equal to the smaller
#' dimension of `image`.
#'
#' @param img A `Image` object
#' @param brush An object created with `make_brush()`
#' @param rel_pos_x,rel_pos_y A relative position to include the brush in the
#'   image. Defaults to 0.5. This means that the brush will be centered in the
#'   original image. Smaller values move the brush toward the left and top,
#'   respectively.
#' @param plot Plots the generated mask? Defaults to `TRUE`.
#'
#' @return A binary image with 0s and 1s.
#' @export
#'
#' @examples
#' img <- image_pliman("soybean_touch.jpg")
#' make_mask(img, brush = make_brush(size = 201))
#' make_mask(img,
#'           brush = make_brush(size = 401, shape = "diamond"),
#'           rel_pos_x = 0.1,
#'           rel_pos_y = 0.8)
make_mask <- function(img,
                      brush,
                      rel_pos_x = 0.5,
                      rel_pos_y = 0.5,
                      plot = TRUE){
  min_dim <- min(dim(img)[1:2])
  if(nrow(brush) > min_dim){
    stop("The size of the brush cannot be greater than the smaller dimension of `image` (", min_dim, ")", call. = FALSE)
  }
  nrim <- nrow(img)
  ncim <- ncol(img)
  nrbrush <- nrow(brush)
  # difference in number of rows
  difr <- nrim - nrbrush
  nrbelow <- trunc(difr * (1 - rel_pos_x))
  nrabove <- difr - nrbelow
  # difference in number of cols
  difcol <- ncim - nrbrush
  ncolleft <- trunc(difcol * (1 - rel_pos_y))
  ncolrigth <- difcol - ncolleft
  # add the rows
  br2 <- rbind(brush, matrix(rep(0, nrow(brush) * nrbelow), ncol = nrow(brush)))
  br2 <- rbind(matrix(rep(0, nrow(brush) * nrabove), ncol = nrow(brush)), br2)
  # add the columns
  br2 <- cbind(br2, matrix(rep(0, nrow(br2) * ncolleft), nrow = nrow(br2)))
  br2 <- cbind(matrix(rep(0, nrow(br2) * ncolrigth), nrow = nrow(br2)), br2)
  mask <- EBImage::Image(br2)
  if(isTRUE(plot)){
    plot(mask)
  }
  invisible(mask)
}



#' Segment an `Image` object using a brush mask
#'
#' It combines [make_mask()] and [make_brush()] to segment an `Image` object
#' using a brush of desired size, shape, and position.
#'
#' @inheritParams make_mask
#' @inheritParams make_brush
#' @param type Defines the type of the mask. By default, a binary mask is
#'   applied. This results in white pixels in the original image that matches
#'   the 0s pixels in the brush. If `type = "shadow"` is used, a shadow mask is produced
#' @param col_background Background color after image segmentation. Defaults to
#'   `"white"`.
#' @return A color `Image` object
#' @export
#'
#' @examples
#' img <- image_pliman("soybean_touch.jpg")
#' plot(img)
#' image_segment_mask(img, size = 601)
#' image_segment_mask(img,
#'                    size = 401,
#'                    shape = "diamond",
#'                    rel_pos_x = 0,
#'                    rel_pos_y = 0,
#'                    type = "shadow")
image_segment_mask <- function(img,
                               size,
                               shape = "disc",
                               rel_pos_x = 0.5,
                               rel_pos_y = 0.5,
                               type = c("binary", "shadow"),
                               col_background = "white",
                               plot = TRUE,
                               ...){
  mask <- make_mask(img,
                    brush = make_brush(size = size, shape = shape, ...),
                    rel_pos_x = rel_pos_x,
                    rel_pos_y = rel_pos_y,
                    plot = FALSE)
  col_background <- col2rgb(col_background) / 255
  ID <- which(mask == 0)
  if(type[[1]] == "binary"){
    img@.Data[, , 1][ID] <- col_background[1]
    img@.Data[, , 2][ID] <- col_background[2]
    img@.Data[, , 3][ID] <- col_background[3]
  } else{
    img@.Data[, , 1][ID] <- 1
  }
  if(isTRUE(plot)){
    plot(img)
  }
  invisible(img)
}

