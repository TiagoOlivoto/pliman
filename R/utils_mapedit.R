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
#' @param quantiles the upper and lower quantiles used for color stretching. If
#'   set to `NULL`, stretching is performed basing on 'domain' argument.
#' @param domain the upper and lower values used for color stretching. This is
#'   used only if `'quantiles'` is `NULL`. If both '`domain'` and `'quantiles'`
#'   are set to `NULL`, stretching is applied based on min-max values.
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
                       alpha = 0.7,
                       attribute = "area",
                       title = "Edit the image",
                       show = c("rgb", "index"),
                       index = "B",
                       max_pixels = 1000000,
                       downsample = NULL,
                       color_regions = custom_palette(),
                       quantiles = c(0, 1),
                       domain = NULL,
                       ...){
  if(!is.null(domain)){
    quantiles <- NULL
  }
  check_mapview()
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

  ras <- terra::rast(EBImage::transpose(img)@.Data) |> stars::st_as_stars()
  dimsto <- dim(ras[,,,1])
  nr <- dimsto[1]
  nc <- dimsto[2]
  npix <- nc * nr
  if(npix > max_pixels){
    if(is.null(downsample)){
      possible_downsamples <- 0:100
      possible_npix <- sapply(possible_downsamples, function(x){
        compute_downsample(nr, nc, x)
      })
      downsample <- which.min(abs(possible_npix - max_pixels)) - 1
      message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
    }

    ras <- terra::rast(
      Map(c,
          lapply(1:3, function(x){
            stars::st_downsample(ras[,,,x], n = downsample) |> terra::rast()
          }))
    ) |>
      stars::st_as_stars()
  }
  if(viewopt == "rgb"){
    sf::st_crs(ras) <- sf::st_crs("EPSG:3857")
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
      mapview::mapview(sf_df,
                       map.types = "OpenStreetMap",
                       col.regions = color_regions,
                       zcol = attribute,
                       legend = TRUE,
                       alpha.regions = alpha,
                       layer.name = attribute) |>
        leafem::addRasterRGB(ras,
                             r = r,
                             g = g,
                             b = b,
                             maxBytes = 64 * 1024 * 1024,
                             domain = domain,
                             quantiles = quantiles)
    } else{
      map <-
        leaflet::leaflet() |>
        leafem::addRasterRGB(ras,
                             r = r,
                             g = g,
                             b = b,
                             maxBytes = 64 * 1024 * 1024,
                             domain = domain,
                             quantiles = quantiles) |>
        mapedit::editMap(editor = "leafpm",
                         title = "title")
      invisible(sf::st_transform(map[["finished"]][["geometry"]], sf::st_crs(ras)))
    }
  } else{
    ind <- mosaic_index(terra::rast(ras), index = index)
    sf::st_crs(ind) <- sf::st_crs("EPSG:3857")

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
#' plot_index(img, index = "B")
#' }
#'
plot_index <- function(img = NULL,
                       object = NULL,
                       index = NULL,
                       remove_bg = TRUE,
                       viewer = get_pliman_viewer(),
                       all_layers = TRUE,
                       layer = 1,
                       max_pixels = 500000,
                       downsample = NULL,
                       downsample_fun = NULL,
                       color_regions = custom_palette(),
                       ncol = NULL,
                       nrow = NULL,
                       aspect_ratio = NA){
  check_mapview()
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
          sto <-  terra::rast(t(img[[i]]@.Data)) |> stars::st_as_stars(proxy = FALSE)
          dimsto <- dim(sto)
          nr <- dimsto[1]
          nc <- dimsto[2]
          npix <- nc * nr
          if(npix > max_pixels){
            if(is.null(downsample)){
              possible_downsamples <- 0:100
              possible_npix <- sapply(possible_downsamples, function(x){
                compute_downsample(nr, nc, x)
              })
              downsample <- which.min(abs(possible_npix - max_pixels)) - 1
              if(i == 1){
                message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
              }
            }
            if(!is.null(downsample_fun)){
              sto <- stars::st_downsample(sto, n = downsample, FUN = downsample_fun) |> terra::rast()
            } else{
              sto <- stars::st_downsample(sto, n = downsample) |> terra::rast()
            }
          } else{
            sto <-  terra::rast(sto)
          }
        }
        )
      )
    } else{
      sts <-  terra::rast(
        lapply(1:terra::nlyr(img), function(i){
          sto <-  img[[i]] |> stars::st_as_stars(proxy = FALSE)
          dimsto <- dim(sto)
          nr <- dimsto[1]
          nc <- dimsto[2]
          npix <- nc * nr
          if(npix > max_pixels){
            if(is.null(downsample)){
              possible_downsamples <- 0:100
              possible_npix <- sapply(possible_downsamples, function(x){
                compute_downsample(nr, nc, x)
              })
              downsample <- which.min(abs(possible_npix - max_pixels)) - 1
              if(i == 1){
                message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
              }
            }
            if(!is.null(downsample_fun)){
              sto <- stars::st_downsample(sto, n = downsample, FUN = downsample_fun) |> terra::rast()
            } else{
              sto <- stars::st_downsample(sto, n = downsample) |> terra::rast()
            }
          } else{
            sto <-  terra::rast(sto)
          }
        }
        )
      )

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
        mapbase <- stars::st_as_stars(sts[[1]])
        sf::st_crs(mapbase) <- sf::st_crs("EPSG:3857")
        mapbase <- mapview::mapview(mapbase,
                                    maxpixels = 50000000,
                                    layer.name = names(mapbase),
                                    map.types = "CartoDB.Positron",
                                    col.regions = color_regions,
                                    alpha.regions = 1,
                                    na.color = "transparent",
                                    verbose = FALSE)

        for (i in 2:terra::nlyr(sts)) {
          lyrtmp <- stars::st_as_stars(sts[[i]])
          sf::st_crs(lyrtmp) <- sf::st_crs("EPSG:3857")
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
        mapbase <- stars::st_as_stars(sts[[layer]])
        sf::st_crs(mapbase) <- sf::st_crs("EPSG:3857")
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
    ind <- image_index(img, index = index, plot = FALSE)[[1]]
    if(!is.null(object)){
      if(isTRUE(remove_bg)){
        ind@.Data[which(mask@.Data == 0)] <- NA
        ras <- terra::rast(t(ind@.Data))
      } else{
        ras <- terra::rast(t(ind@.Data))
      }
    } else{
      ras <- terra::rast(t(ind@.Data))
    }
    sto <-  ras |> stars::st_as_stars(proxy = FALSE)
    dimsto <- dim(sto)
    nr <- dimsto[1]
    nc <- dimsto[2]
    npix <- nc * nr
    if(npix > max_pixels){
      if(is.null(downsample)){
        possible_downsamples <- 0:100
        possible_npix <- sapply(possible_downsamples, function(x){
          compute_downsample(nr, nc, x)
        })
        downsample <- which.min(abs(possible_npix - max_pixels)) - 1
        message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
      }
      if(!is.null(downsample_fun)){
        ras <- stars::st_downsample(sto, n = downsample, FUN = downsample_fun) |> terra::rast()
      } else{
        ras <- stars::st_downsample(sto, n = downsample) |> terra::rast()
      }
    } else{
      ras <-  terra::rast(sto)
    }
    if(vieweropt == "base"){
      terra::plot(ras,
                  axes = FALSE,
                  loc.main = "topleft",
                  cex.main = 1,
                  smooth = TRUE)
    } else{
      stob <- stars::st_as_stars(ras)
      sf::st_crs(stob) <- sf::st_crs("EPSG:3857")
      suppressWarnings(
        suppressMessages(
          mapview::mapview(stob,
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
#' @return The alighed/cropped image for further visualization or analysis.
#'
#' @examples
#' # Example usage:
#' if(interactive()){
#' img <- image_pliman("mult_leaves.jpg")
#' image_prepare(img, viewer = "mapview")
#'}
#' @export
image_prepare <- function(img, viewer = get_pliman_viewer()){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]

  align <- image_align(img, viewer = vieweropt)
  if(vieweropt != "base"){
    image_view(img[1:10, 1:10, ])
  }
  cropped <- image_crop(align, viewer = vieweropt)
  invisible(cropped)
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
                         title = "Use the 'Draw Rectangle' tool to create a rectangle in the image"){
  e <- image_view(img, show = show, title = title, index = index)
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
                       title = "Use the 'Draw Polygon' tool to create a polygon in the image"){
  e <- image_view(img, show = show, title = title, index = index)
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
                      title = "Use the 'Draw Marker' tool to select points in the image"){
  e <- image_view(img, show = show, title = title, index = index)
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



mosaic_index <- function(mosaic,
                         index = "R",
                         r = 1,
                         g = 2,
                         b = 3,
                         nir = 4,
                         re = 5){
  if(inherits(mosaic, "Image")){
    ras <- t(terra::rast(mosaic@.Data))
  } else{
    ras <- mosaic
  }
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
  if (!index %in% ind$Index) {
    message(paste("Index '", index, "' is not available. Trying to compute your own index.",
                  sep = ""))
  }
  R <- try(ras[[r]], TRUE)
  G <- try(ras[[g]], TRUE)
  B <- try(ras[[b]], TRUE)
  NIR <- try(ras[[nir]], TRUE)
  RE <- try(ras[[re]], TRUE)
  if(index %in% ind$Index){
    mosaic_gray <-
      eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index]))) |>
      stars::st_as_stars()
  } else{
    mosaic_gray <-
      eval(parse(text = as.character(index))) |>
      stars::st_as_stars()
  }
  names(mosaic_gray) <- index
  if(!is.na(sf::st_crs(mosaic))){
    suppressWarnings(sf::st_crs(mosaic_gray) <- sf::st_crs(mosaic))
  } else{
    suppressWarnings(sf::st_crs(mosaic_gray) <- "+proj=utm +zone=32 +datum=WGS84 +units=m")
  }
  invisible(mosaic_gray)
}

# image_orientation: Determines the orientation of an image (vertical or horizontal).
image_orientation <- function(img){
  ifelse(ncol(img) > nrow(img), "vertical", "horizontal")
}
