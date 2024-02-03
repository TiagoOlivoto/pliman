#' Create an interactive map view of an image
#'
#' This function allows users to interactively edit and analyze an image using
#' mapview and mapedit packages.
#'
#' @inheritParams plot_index
#' @param img An `Image` object.
#' @param object (Optional). An object computed with [analyze_objects()]. If an
#'   object is informed, an additional layer is added to the plot, showing the
#'   contour of the analyzed objects, with a color gradient defined by
#'   `attribute`.
#' @param r,g,b The layer for the Red, Green and Blue band, respectively.
#'   Defaults to `1`, `2`, and `3`.
#' @param edit If `TRUE` enable editing options using [mapedit::editMap()].
#' @param alpha The transparency level of the rectangles' color (between 0 and 1).
#' @param attribute The name of the quantitative variable in the
#'   \code{object_index} to be used for coloring the rectangles.
#' @param title The title of the map view. Use to provide short orientations to
#'   the user.
#' @param show The display option for the map view. Options are "rgb" for RGB
#'   view and "index" for index view.
#' @param index The index to use for the index view. Defaults to "B".
#' @param max_pixels integer > 0. Maximum number of cells to use for the plot.
#'   If `max_pixels < npixels(img)`, regular sampling is used before plotting.
#' @param color_regions The color palette for displaying index values. Default
#'   is [custom_palette()].
#' @param quantiles the upper and lower quantiles used for color stretching. Set
#'   to `c(0, 1)`
#' @param ... Additional arguments to be passed to `downsample_fun`.
#' @return An `sf` object, the same object returned by [mapedit::editMap()].
#'
#' @examples
#' if(interactive()){
#' # Example usage:
#' img <- image_pliman("sev_leaf.jpg")
#' image_view(img)
#' }
#'
#' @export
#'
image_view <- function(img,
                       object = NULL,
                       r = 1,
                       g = 2,
                       b = 3,
                       edit = FALSE,
                       alpha = 0.7,
                       attribute = "area",
                       title = "Edit the image",
                       show = c("rgb", "index"),
                       index = "B",
                       max_pixels = 1000000,
                       downsample = NULL,
                       color_regions = custom_palette(),
                       quantiles = c(0, 1),
                       ...){
  # check_mapview()
  viewopt <- c("rgb", "index")
  viewopt <- viewopt[pmatch(show[[1]], viewopt)]
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

  ras <- terra::rast(EBImage::transpose(img)@.Data)
  nly <- terra::nlyr(ras)
  terra::crs(ras) <- terra::crs("EPSG:3857")
  dimsto <- dim(ras)[1:2]
  nr <- dimsto[1]
  nc <- dimsto[2]
  npix <- nc * nr
  if(npix > max_pixels){
    possible_downsamples <- 0:50
    possible_npix <- sapply(possible_downsamples, function(x){
      compute_downsample(nr, nc, x)
    })
    if(is.null(downsample)){
      downsample <- which.min(abs(possible_npix - max_pixels))
      downsample <- ifelse(downsample == 1, 0, downsample)
    }
    if(downsample > 0){
      message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
      ras <- mosaic_aggregate(ras, pct = round(100 / downsample))
    }
  }
  if(viewopt == "rgb"){
    if(nly >= 3){
      if(!is.null(object)){
        sf_df <- sf::st_sf(
          geometry = lapply(object$contours, function(x) {
            tmp <- x
            tmp[, 2] <- ncol(img) - tmp[, 2]
            sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(get_measures(object)),
          crs = sf::st_crs("EPSG:3857")
        )
        colnames(sf_df) <- gsub("data.", "", colnames(sf_df))

        mapview::viewRGB(
          as(ras, "Raster"),
          layer.name = "base",
          r = r,
          g = g,
          b = b,
          na.color = "#00000000",
          maxpixels = 60000000,
          quantiles = quantiles
        ) +
          mapview::mapview(sf_df,
                           map.types = "OpenStreetMap",
                           # col.regions = color_regions,
                           zcol = attribute,
                           legend = TRUE,
                           alpha.regions = alpha,
                           layer.name = attribute)
      } else{
        if(isTRUE(edit)){
          map <-
            mapview::viewRGB(
              as(ras, "Raster"),
              layer.name = "base",
              r = r,
              g = g,
              b = b,
              na.color = "#00000000",
              maxpixels = 60000000,
              quantiles = quantiles
            ) |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
          invisible(sf::st_transform(map[["finished"]][["geometry"]], sf::st_crs(ras)))
        } else{
          mapview::viewRGB(
            as(ras, "Raster"),
            layer.name = "base",
            r = r,
            g = g,
            b = b,
            na.color = "#00000000",
            maxpixels = 60000000,
            quantiles = quantiles
          )
        }
      }
    } else{
      if(isTRUE(edit)){
        map <-
          mapview::mapview(ras,
                           map.types = "OpenStreetMap",
                           layer.name = "layer",
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE) |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
        invisible(sf::st_transform(map[["finished"]][["geometry"]], sf::st_crs(ras)))
      } else{
        map <-
          mapview::mapview(ras,
                           map.types = "OpenStreetMap",
                           layer.name = "layer",
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
        map
      }
    }
  } else{
    ind <- mosaic_index(ras, index = index)
    terra::crs(ind) <- terra::crs("EPSG:3857")
    map <-
      suppressWarnings(
        suppressMessages(
          mapview::mapview(ind,
                           map.types = "CartoDB.Positron",
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE) |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        )
      )
    invisible(sf::st_transform(map[["finished"]][["geometry"]], sf::st_crs(ras)))
  }
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
custom_palette <- function(colors = c("yellow", "#53CC67", "#009B95", "#00588B","#4B0055"), n = 5){
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
#' @param object An object computed with [analyze_objects()] using the argument
#'   `return_mask = TRUE`.
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
#' @param all_layers Render all layers when `img` is an object computed with
#'   [image_index()] and `viewer = "mapview"`?.
#' @param layer The layer to plot when `img` is an object computed with
#'   [image_index()] and `viewer = "mapview"`. Defaults to the first layer
#'   (first index computed).
#' @param max_pixels integer > 0. Maximum number of cells to plot the index. If
#'   `max_pixels < npixels(img)`, downsampling is performed before plotting the
#'   index. Using a large number of pixels may slow down the plotting time.
#' @param downsample integer; for each dimension the number of
#'   pixels/lines/bands etc that will be skipped; Defaults to `NULL`, which will
#'   find the best downsampling factor to approximate the `max_pixels` value.
#' @param downsample_fun function; if given, downsampling will apply
#'   `downsample_fun`` ` to each of the the subtiles.
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
#' plot_index(img, index = c("R", "G"))
#' }
#'
plot_index <- function(img = NULL,
                       object = NULL,
                       index = NULL,
                       remove_bg = TRUE,
                       viewer = get_pliman_viewer(),
                       all_layers = TRUE,
                       layer = 1,
                       max_pixels = 1000000,
                       downsample = NULL,
                       downsample_fun = NULL,
                       color_regions = custom_palette(n = 100),
                       ncol = NULL,
                       nrow = NULL,
                       aspect_ratio = NA){
  # check_mapview()
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
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]

  if(!is.null(img) & inherits(img, c("SpatRaster", "image_index"))){
    if(inherits(img, "image_index")){
      for(x in 1:length(img)){
        img[[x]][is.infinite(img[[x]])] <- NA
      }
      sts <-  terra::rast(
        lapply(1:length(img), function(i){
          sto <-  terra::rast(t(img[[i]]@.Data))
          dimsto <- dim(sto)
          nr <- dimsto[1]
          nc <- dimsto[2]
          npix <- nc * nr
          if(npix > max_pixels){
              possible_downsamples <- 0:50
              possible_npix <- sapply(possible_downsamples, function(x){
                compute_downsample(nr, nc, x)
              })
              if(is.null(downsample)){
                downsample <- which.min(abs(possible_npix - max_pixels))
                downsample <- ifelse(downsample == 1, 0, downsample)
              }
              if(downsample > 0){
                sto <- mosaic_aggregate(sto, pct = round(100 / downsample))
              }
          }
          sto
        }
        )
      )
    } else{
      dimsto <- dim(img)
      nr <- dimsto[1]
      nc <- dimsto[2]
      npix <- nr * nc
      if(npix > max_pixels){
        possible_downsamples <- 0:50
        possible_npix <- sapply(possible_downsamples, function(x){
          compute_downsample(nr, nc, x)
        })
        if(is.null(downsample)){
          downsample <- which.min(abs(possible_npix - max_pixels))
          downsample <- ifelse(downsample == 1, 0, downsample)
        }
        if(downsample > 0){
          sts <- mosaic_aggregate(img, pct = round(100 / downsample))
        }
      } else{
        sts <- img
      }
    }

    names(sts) <- names(img)

    num_plots <- terra::nlyr(sts)
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
      if(terra::nlyr(sts) > 16){
        warning("The number of layers is too large and plots may not fit well to the plot area. Consider reducing the number of indexes used.", call. = FALSE)
      }
      terra::plot(sts,
                  col = color_regions,
                  axes = FALSE,
                  nc = ncol,
                  nr = nrow,
                  loc.main = "topleft",
                  cex.main = 1,
                  smooth = TRUE)
    } else{
      if(layer > terra::nlyr(sts)){
        warning("The layer number is greater than the total number of layers. Plotting the first layer.", call. = FALSE)
        layer <- 1
      }
      if(terra::nlyr(sts) == 1 & all_layers == TRUE){
        all_layers <- FALSE
      }
      if(isTRUE(all_layers)){
        mapbase <- sts[[1]]
        terra::crs(mapbase) <- terra::crs("EPSG:3857")
        mapbase <- mapview::mapview(mapbase,
                                    maxpixels = 50000000,
                                    layer.name = names(mapbase),
                                    map.types = "CartoDB.Positron",
                                    col.regions = color_regions,
                                    alpha.regions = 1,
                                    na.color = "transparent",
                                    verbose = FALSE)

        for (i in 2:terra::nlyr(sts)) {
          lyrtmp <- sts[[i]]
          terra::crs(lyrtmp) <- terra::crs("EPSG:3857")
          mapbase <- mapview::mapview(lyrtmp,
                                      maxpixels = 50000000,
                                      map = mapbase,
                                      hide = TRUE,
                                      layer.name = names(lyrtmp),
                                      map.types = "CartoDB.Positron",
                                      col.regions = color_regions,
                                      alpha.regions = 1,
                                      na.color = "transparent",
                                      verbose = FALSE)

        }
      } else{
        mapbase <- sts[[layer]]
        terra::crs(mapbase) <- terra::crs("EPSG:3857")
        mapbase <- mapview::mapview(mapbase,
                                    maxpixels = 50000000,
                                    layer.name = names(mapbase),
                                    map.types = "CartoDB.Positron",
                                    col.regions = color_regions,
                                    alpha.regions = 1,
                                    na.color = "transparent",
                                    verbose = FALSE)
      }
      mapbase
    }
  } else{
    if(!is.null(object)){
      if(is.null(object$mask)){
        stop("Use `return_mask = TRUE` in `analyze_objects()` to plot the image index.")
      }
      mask <- object$mask
    }
    if(!is.null(index)){
      index <- index
    } else if(!missing(object) & !is.null(object$parms$object_index)){
      index <- object$parms$object_index[[1]]
    } else{
      index <- "B"
    }
    ind <- image_index(img, index = index, plot = FALSE)[[1]]
    if(!is.null(object)){
      if(isTRUE(remove_bg)){
        ind@.Data[which(mask@.Data == 0)] <- NA
        ras <- terra::rast(EBImage::transpose(ind)@.Data)
      } else{
        ras <- terra::rast(EBImage::transpose(ind)@.Data)
      }
    } else{
      ras <-terra::rast(EBImage::transpose(ind)@.Data)
    }
    dimsto <- dim(ras)
    nr <- dimsto[1]
    nc <- dimsto[2]
    npix <- nc * nr
    if(npix > max_pixels){
      possible_downsamples <- 0:50
      possible_npix <- sapply(possible_downsamples, function(x){
        compute_downsample(nr, nc, x)
      })
      if(is.null(downsample)){
        downsample <- which.min(abs(possible_npix - max_pixels))
        downsample <- ifelse(downsample == 1, 0, downsample)
      }
      if(downsample > 0){
        ras <- mosaic_aggregate(ras, pct = round(100 / downsample))
      }
    }
    if(vieweropt == "base"){
      terra::plot(ras,
                  axes = FALSE,
                  loc.main = "topleft",
                  cex.main = 1,
                  smooth = TRUE)
    } else{
      terra::crs(ras) <- terra::crs("EPSG:3857")
      suppressWarnings(
        suppressMessages(
          mapview::mapview(ras,
                           maxpixels = 50000000,
                           layer.name = index,
                           map.types = "CartoDB.Positron",
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE)
        )
      )
    }
  }
}

#' Prepare an image
#'
#' This function aligns and crops the image using either base or mapview
#' visualization. This is useful to prepare the images to be analyzed with
#' [analyze_objects_shp()]
#' @inheritParams image_view
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
#' @return The alighed/cropped image for further visualization or analysis.
#'
#' @examples
#' # Example usage:
#' if(interactive()){
#' img <- image_pliman("mult_leaves.jpg")
#' image_prepare(img, viewer = "mapview")
#'}
#' @export
image_prepare <- function(img,
                          viewer = get_pliman_viewer(),
                          downsample = NULL,
                          max_pixels = 1000000){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]

  align <- image_align(img, viewer = vieweropt)
  if(vieweropt != "base"){
    image_view(img[1:5, 1:5, ], edit = TRUE)
  }
  cropped <- image_crop(align, viewer = vieweropt, downsample = downsample, max_pixels = max_pixels)
  invisible(cropped)
}




# mv_two_points: Allows the user to select two points in an image and returns their coordinates.
mv_two_points <- function(img,
                          show = "rgb",
                          index = "NGRDI",
                          title = "Use the 'Draw Polyline' tool to Select 2 points",
                          downsample = NULL,
                          max_pixels = 1000000){
  e <- image_view(img,
                  show = show,
                  title = title,
                  index = index,
                  downsample = downsample,
                  max_pixels = max_pixels,
                  edit = TRUE)
  if(!inherits(e, "sfc_LINESTRING")){
    stop("The geometry used is not valid. Please, use 'Draw Polyline' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  x1 <- e[[1]][1]
  x2 <- e[[1]][2]
  x3 <- e[[1]][3]
  x4 <- e[[1]][4]
  invisible(list(x1 = x1,
                 y1 = nc - x3,
                 x2 = x2,
                 y2 = nc - x4))
}

# mv_rectangle: Enables the user to create a rectangle in an image and retrieves its coordinates.
mv_rectangle <- function(img,
                         show = "rgb",
                         index = "NGRDI",
                         title = "Use the 'Draw Rectangle' tool to create a rectangle in the image",
                         downsample = NULL,
                         max_pixels = 1000000){
  e <- image_view(img,
                  show = show,
                  title = title,
                  index = index,
                  downsample = downsample,
                  max_pixels = max_pixels,
                  edit = TRUE)
  if(!inherits(e, "sfc_POLYGON")){
    stop("The geometry used is not valid. Please, use 'Draw Rectangle' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  coords <- e[[1]][[1]]
  coords[, 2] <- nc - coords[, 2]
  colnames(coords) <- c("x", "y")
  invisible(coords)
}

# mv_polygon: Allows the user to create a polygon in an image and retrieves its coordinates.
mv_polygon <- function(img,
                       show = "rgb",
                       index = "NGRDI",
                       title = "Use the 'Draw Polygon' tool to create a polygon in the image",
                       downsample = NULL,
                       max_pixels = 1000000){
  e <- image_view(img,
                  show = show,
                  title = title,
                  index = index,
                  downsample = downsample,
                  max_pixels = max_pixels,
                  edit = TRUE)
  if(!inherits(e, "sfc_POLYGON")){
    stop("The geometry used is not valid. Please, use 'Draw Polygon' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  coords <- e[[1]][[1]]
  colnames(coords) <- c("x", "y")
  coords[, 2] <- nc - coords[, 2]
  invisible(coords)
}

# mv_points: Enables the user to select multiple points in an image and returns their coordinates.
mv_points <- function(img,
                      show = "rgb",
                      index = "NGRDI",
                      title = "Use the 'Draw Marker' tool to select points in the image",
                      downsample = NULL,
                      max_pixels = 1000000){
  e <- image_view(img,
                  show = show,
                  title = title,
                  index = index,
                  downsample = downsample,
                  max_pixels = max_pixels,
                  edit = TRUE)
  if(!inherits(e, "sfc_POINT")){
    stop("The geometry used is not valid. Please, use 'Draw Marker' tool to select two points.", call. = FALSE)
  }
  nc <- ncol(img)
  coords <-
    do.call(rbind, lapply(e, function(x){
      as.vector(x)
    }))
  colnames(coords) <- c("x", "y")
  coords[, 2] <- nc - coords[, 2]
  invisible(coords)
}

