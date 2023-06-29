#' Create an interactive map view of an image
#'
#' This function allows users to interactively edit and analyze an image using
#' mapview and mapedit packages.
#'
#' @param img An `Image` object.
#' @param title The title of the map view. Use to provide short orientations to
#'   the user.
#' @param show The display option for the map view. Options are "rgb" for RGB
#'   view and "index" for index view.
#' @param index The index to use for the index view. Defaults to "B".
#' @param max_pixels integer > 0. Maximum number of cells to use for the plot.
#'   If `max_pixels < npixels(img)`, regular sampling is used before plotting.
#' @param color_regions The color palette for displaying index values. Default
#'   is [custom_palette()].
#' @param ... Additional arguments to be passed to [mapview::mapview()].
#' @return An `sf` object, the same object returned by [mapedit::editMap()].
#'
#' @importFrom raster raster brick
#' @examples
#' if(interactive()){
#' # Example usage:
#' img <- image_pliman("sev_leaf.jpg")
#' image_view(img)
#' }
#'
#' @export
image_view <- function(img,
                       title = "Edit the image",
                       show = c("rgb", "index"),
                       index = "B",
                       max_pixels = 2500000,
                       color_regions = custom_palette(),
                       ...){
  check_mapview()
  viewopt <- c("rgb", "index")
  viewopt <- viewopt[pmatch(show[[1]], viewopt)]
  if(npixels(img) > max_pixels){
    img <- reduce_dimensions(img, target_pixels = max_pixels)
  }
  if(viewopt == "rgb"){
    ras <- rgb_to_raster(img)
    map <-
      leaflet::leaflet() |>
      leafem::addRasterRGB(ras,
                           r = 1,
                           g = 2,
                           b = 3,
                           maxBytes = 4 * 2048 * 2048) |>
      mapedit::editMap(editor = "leafpm",
                       title = title)
  } else{
    ind <- help_imageindex(img, index = index)
    ras2 <- raster::raster(t(ind@.Data))
    map <-
      suppressWarnings(
        suppressMessages(
          mapview::mapview(ras2,
                           map.types = "CartoDB.Positron",
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE,
                           ...) |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        )
      )
  }
  return(map[["finished"]][["geometry"]])
}

#' Generate Custom Color Palette
#'
#' This function generates a custom color palette using the specified colors and
#' number of colors.
#'
#' @param colors A vector of colors to create the color palette. Default is
#'   c("steelblue", "salmon", "forestgreen").
#' @param n The number of gradient colors in the color palette. Default is 100.
#'
#' @return A vector of colors representing the custom color palette.
#'
#' @examples
#' # Generate a custom color palette with default colors and 10 colors
#' custom_palette()
#'
#' # Generate a custom color palette with specified colors and 20 colors
#' custom_palette(colors = c("blue", "red"), n = 20)
#'
#' @importFrom grDevices colorRampPalette
#'
#' @export
#' @examples
#' # example code
#' library(pliman)
#' custom_palette(n = 5)
#'
custom_palette <- function(colors = c("#4B0055", "#00588B", "#009B95", "#53CC67", "yellow"), n = 100){
  grDevices::colorRampPalette(colors)(n)
}

#' Plot an image index
#'
# This function plots the specified index of an image either using base plotting
# or mapview package.
#'
#' @param img An optional `Image` object or an object computed with
#'   [image_index()]. If `object` is provided, then the input image is obtained
#'   internally.
#' @param object An object computed with [analyze_objects_shp()]. By using this
#'   object you can ignore `img`.
#' @param index The index to plot. Defaults to the index computed from the
#'   `object` if provided. Otherwise, the `B` index is computed. See
#'   [image_index()] for more details.
#' @param remove_bg Logical value indicating whether to remove the background
#'   when `object` is provided. Defaults to `TRUE`.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param layer The layer to plot when `img` is an object computed with
#'   [image_index()] and `viewer = "mapview"`. Defaults to the first layer
#'   (first index computed).
#' @param max_pixels integer > 0. Maximum number of cells to plot the index. If
#'   `max_pixels < npixels(img)`, regular sampling is used before plotting.
#'   Using a large number of pixels may slow down the plotting time.
#' @param color_regions The color palette for displaying index values. Default
#'   is [custom_palette()].
#' @param nrow,ncol The number of rows or columns in the plot grid. Defaults to
#'   `NULL`, i.e., a square grid is produced.
#' @param aspect_ratio Numeric, giving the aspect ratio y/x. Defaults to `NA`.
#'   See [graphics::plot.window()] for more details.
#' @return None
#' @export
#'
#' @examples
#' if(interactive()){
#' # Example usage:
#' library(pliman)
#' img <- image_pliman("sev_leaf.jpg")
#' plot_index(img, index = "B")
#' }
#'
plot_index <- function(img = NULL,
                       object = NULL,
                       index = NULL,
                       remove_bg = TRUE,
                       viewer = get_pliman_viewer(),
                       layer = 1,
                       max_pixels = 1000000,
                       color_regions = custom_palette(),
                       ncol = NULL,
                       nrow = NULL,
                       aspect_ratio = NA){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(!is.null(img) & inherits(img, "image_index")){
    rast <- lapply(img, function(x){
      if(npixels(x) > max_pixels){
        x <- reduce_dimensions(x, max_pixels)
      }
      raster::raster(t(x@.Data))
    })

    num_plots <-length(rast)
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
    if(vieweropt == "base"){
      rbrick <- raster::brick(rast)
      if(raster::nlayers(rbrick) > 16){
        warning("The number of layers is too large and plots may not fit well to the plot area. Consider reducing the number of indexes used.", call. = FALSE)
      }
      raster::plot(rbrick,
                   axes = FALSE,
                   nc = ncol,
                   nr = nrow,
                   asp = aspect_ratio,
                   maxnl = raster::nlayers(rbrick),
                   col = color_regions)
    } else{
      if(layer > length(rast)){
        warning("The layer number is greater than the total number of layers. Plotting the first layer.", call. = FALSE)
        layer <- 1
      }
      suppressWarnings(
        suppressMessages(
          mapview::mapview(rast[[layer]],
                           layer.name = names(rast[layer]),
                           map.types = "CartoDB.Positron",
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE)
        )
      )
    }
  } else{
    if(!is.null(object)){
      img <- object$final_image
      mask <- object$mask
    } else if(is.null(img)){
      stop("One of 'img' or 'object' must be informed.", call. = FALSE)
    }
    if(!is.null(index)){
      index <- index
    } else if(!missing(object) & !is.null(object$object_index_computed)){
      index <- object$object_index_computed[[1]]
    } else{
      index <- "B"
    }
    ind <- help_imageindex(img, index = index)
    if(!is.null(object)){
      if(isTRUE(remove_bg)){
        ind@.Data[which(mask@.Data == 0)] <- NA
        ras <- raster::raster(t(ind@.Data))
      } else{
        ras <- raster::raster(t(ind@.Data))
      }
    } else{
      ras <- raster::raster(t(ind@.Data))
    }
    if(vieweropt == "base"){
      raster::plot(ras, col = color_regions,  axes = FALSE, asp = aspect_ratio)
    } else{
      suppressWarnings(
        suppressMessages(
          mapview::mapview(ras,
                           layer.name = index,
                           map.types = "CartoDB.Positron",
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE)
        )
      )
    }
  }
}

#' Prepare an image (align and crop)
#'
#' This function aligns and crops the image using either base or mapview
#' visualization. This is useful to prepare the images to be analyzed with
#' [analyze_objects_shp()]
#'
#' @param img An optional `Image` object
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @importFrom raster extent<- projection crs
#' @return The alighed/cropped image for further visualization or analysis.
#'
#' @examples
#' # Example usage:
#' if(interactive()){
#' img <- image_pliman("mult_leaves.jpg")
#' image_prepare_mv(img, viewer = "mapview")
#'}
#' @export
image_prepare_mv <- function(img,
                             viewer = get_pliman_viewer()){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]

  align <- image_align(img, viewer = vieweropt)
  if(vieweropt != "base"){
    image_view(img |> reduce_dimensions(1000))
  }
  cropped <- image_crop(align, viewer = vieweropt)
  return(cropped)
}




# rgb_to_raster: Converts an RGB image object to a raster brick object.
rgb_to_raster <- function(img){
  r <- t(img@.Data[,,1])
  g <- t(img@.Data[,,2])
  b <- t(img@.Data[,,3])
  # Get the dimensions of the matrices
  nrows <- nrow(r)
  ncols <- ncol(r)
  # Create the RGB array
  rgb_array <- array(0, dim = c(nrows, ncols, 3))
  rgb_array[, , 1] <- r
  rgb_array[, , 2] <- g
  rgb_array[, , 3] <- b
  # Calculate the aspect ratio of the original image
  aspect_ratio <- nrows / ncols
  # Determine the extent based on the aspect ratio
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- aspect_ratio
  # Create a raster brick from the RGB array
  ras <- raster::brick(rgb_array)
  # Set the extent of the raster brick
  raster::extent(ras) <- c(xmin, xmax, ymin, ymax)
  # Set the spatial reference system (SRS) information for the raster brick
  raster::projection(ras) <- raster::crs("+proj=longlat +datum=WGS84")
  return(ras)
}

img_scale <- function(img){
  nr <- nrow(img)
  nc <- ncol(img)
  if(image_orientation(img) == "vertical"){
    scale <- min(nr, nc)
  } else{
    scale <- max(nr, nc)
  }
  return(scale)
}

# mv_two_points: Allows the user to select two points in an image and returns their coordinates.
mv_two_points <- function(img,
                          show = "rgb",
                          index = "NGRDI",
                          title = "Use the 'Draw Polyline' tool to Select 2 points"){
  e <- image_view(img, show = show, title = title, index = index)
  if(!inherits(e, "sfc_LINESTRING")){
    stop("The geometry used is not valid. Please, use 'Draw Polyline' tool to select two points.", call. = FALSE)
  }
  x1 <- e[[1]][1]
  x2 <- e[[1]][2]
  x3 <- e[[1]][3]
  x4 <- e[[1]][4]
  nc <- ncol(img)
  scale <- img_scale(img)
  return(list(x1 = x1 * scale,
              y1 = nc - (x3 * scale),
              x2 = x2 * scale,
              y2 = nc - (x4 * scale)))
}

# mv_rectangle: Enables the user to create a rectangle in an image and retrieves its coordinates.
mv_rectangle <- function(img,
                         show = "rgb",
                         index = "NGRDI",
                         title = "Use the 'Draw Rectangle' tool to create a rectangle in the image"){
  e <- image_view(img, show = show, title = title, index = index)
  if(!inherits(e, "sfc_POLYGON")){
    stop("The geometry used is not valid. Please, use 'Draw Rectangle' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  scale <- img_scale(img)
  coords <- e[[1]][[1]] * scale
  colnames(coords) <- c("x", "y")
  coords[, 2] <- nc - coords[, 2]
  return(coords)
}

# mv_polygon: Allows the user to create a polygon in an image and retrieves its coordinates.
mv_polygon <- function(img,
                       show = "rgb",
                       index = "NGRDI",
                       title = "Use the 'Draw Polygon' tool to create a polygon in the image"){
  e <- image_view(img, show = show, title = title, index = index)
  if(!inherits(e, "sfc_POLYGON")){
    stop("The geometry used is not valid. Please, use 'Draw Polygon' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  scale <- img_scale(img)
  coords <- e[[1]][[1]] * scale
  colnames(coords) <- c("x", "y")
  coords[, 2] <- nc - coords[, 2]
  return(coords)
}

# mv_points: Enables the user to select multiple points in an image and returns their coordinates.
mv_points <- function(img,
                      show = "rgb",
                      index = "NGRDI",
                      title = "Use the 'Draw Marker' tool to select points in the image"){
  e <- image_view(img, show = show, title = title, index = index)
  if(!inherits(e, "sfc_POINT")){
    stop("The geometry used is not valid. Please, use 'Draw Marker' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  scale <- img_scale(img)
  coords <-
    do.call(rbind, lapply(e, function(x){
      as.vector(x)
    })) * scale
  colnames(coords) <- c("x", "y")
  coords[, 2] <- nc - coords[, 2]
  return(coords)
}


# reduce_dimensions: Resizes an image by reducing its dimensions while maintaining the aspect ratio.
reduce_dimensions <- function(img, target_pixels = 2500000) {
  original_rows <- dim(img)[1]
  original_cols <- dim(img)[2]
  original_aspect_ratio <- original_rows / original_cols
  original_total_pixels <- original_rows * original_cols
  reduction_factor <- sqrt(target_pixels / original_total_pixels)
  new_rows <- round(original_rows * reduction_factor)
  EBImage::resize(img, new_rows)
}



# image_orientation: Determines the orientation of an image (vertical or horizontal).
image_orientation <- function(img){
  ifelse(ncol(img) > nrow(img), "vertical", "horizontal")
}
