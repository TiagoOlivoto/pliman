#' Mosaic View
#'
#' @details The function can generate either an interactive map using the
#'   'mapview' package or a static plot using the 'base' package, depending on
#'   the `viewer` and `show` parameters. If show = "index" is used, the function
#'   first computes an image index that can be either an RGB-based index or a
#'   multispectral index, if a multispectral mosaic is provided.
#'
#' @param mosaic A mosaic of class `SpatRaster`, generally imported with
#'   [mosaic_input()].
#' @inheritParams image_view
#' @inheritParams image_align
#' @param r The layer for the Red band (default: 3).
#' @param g The layer for the Green band (default: 2).
#' @param b The layer for the Blue band (default: 1).
#' @param re The layer for the Red-edge band (default: 4).
#' @param nir The layer for the Near-infrared band(default: 5).
#' @param title A title for the generated map or plot (default: "").
#' @param max_pixels Maximum number of pixels to render in the map or plot
#'   (default: 500000).
#' @param downsample Downsampling factor to reduce the number of pixels
#'   (default: NULL). In this case, if the number of pixels in the image (width
#'   x height) is greater than `max_pixels` a downsampling factor will be
#'   automatically chosen so that the number of plotted pixels approximates the
#'   `max_pixels`.
#' @param alpha opacity of the fill color of the raster layer(s).
#' @param quantiles the upper and lower quantiles used for color stretching. If
#'   set to `NULL`, stretching is performed basing on 'domain' argument.
#' @param domain the upper and lower values used for color stretching. This is
#'   used only if `'quantiles'` is `NULL`. If both '`domain'` and `'quantiles'` are
#'   set to `NULL`, stretching is applied based on min-max values.
#' @param axes logical. Draw axes? Defaults to `FALSE`.
#' @param ... Additional arguments passed on to [terra::plot()] when `viewer =
#'   "base"`.
#' @return An sf object, the same object returned by [mapedit::editMap()].
#'
#' @importFrom terra rast crs nlyr
#' @importFrom sf st_crs st_transform st_make_grid
#' @importFrom stars st_downsample st_as_stars
#' @examples
#' if(interactive()){
#' library(pliman)
#' # Load a raster showing the elevation of Luxembourg
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # Generate an interactive map using 'mapview'
#' mosaic_view(mosaic)
#'
#' # Generate a static plot using 'base'
#' mosaic_view(mosaic, viewer = "base")
#' }
#'
#'
#' @export
mosaic_view <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        re = 4,
                        nir = 5,
                        title = "",
                        viewer = c("mapview", "base"),
                        show = c("rgb", "index"),
                        index = "B",
                        max_pixels = 500000,
                        downsample = NULL,
                        alpha = 1,
                        quantiles = c(0, 1),
                        domain = NULL,
                        color_regions = custom_palette(),
                        axes = FALSE,
                        ...){
  check_mapview()
  if(!is.null(domain)){
    quantiles <- NULL
  }
  mapview::mapviewOptions(layers.control.pos = "topright")
  on.exit(mapview::mapviewOptions(default = TRUE))
  viewopt <- c("rgb", "index")
  viewopt <- viewopt[pmatch(show[[1]], viewopt)]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]
  if(viewopt == "rgb" & vieweropt == "base" & terra::nlyr(mosaic) > 1){
    message("`viewer = 'base' can only be used with `show = 'index'`. Defaulting to viewer = 'mapview'")
    vieweropt <- "mapview"
  }
  if(inherits(mosaic, "Image")){
    mosaic <- t(terra::rast(mosaic@.Data))
  }
  sto <- suppressWarnings(stars::st_as_stars(mosaic, proxy = FALSE))
  dimsto <- dim(sto)
  nr <- dimsto[1]
  nc <- dimsto[2]
  npix <- nc * nr
  if(max_pixels > 500000){
    message("The number of pixels is too high, which might slow the rendering process.")
  }
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
  possible_downsamples <- 0:100
  possible_npix <- sapply(possible_downsamples, function(x){
    compute_downsample(nr, nc, x)
  })
  if(is.null(downsample)){
    downsample <- which.min(abs(possible_npix - max_pixels)) - 1
  }
  if(downsample > 0){
    message(paste0("Using downsample = ", downsample, " so that the number of rendered pixels approximates the `max_pixels`"))
    if(terra::nlyr(mosaic) == 1){
      sto <- stars::st_downsample(sto, n = downsample)
    } else{
      downsampled <-
        lapply(1:terra::nlyr(mosaic), function(x){
          stars::st_downsample(sto[,,,x], n = downsample) |> terra::rast()
        })

      sto <- terra::rast(Map(c, downsampled)) |> stars::st_as_stars()
    }
  }

  sto[sto == 65535] <- NA

  if(viewopt == "rgb"){
    sto <- sto
  } else{
    if(terra::nlyr(mosaic) > 2){
      sto <- mosaic_index(terra::rast(sto), index = index)
    } else{
      sto <- sto
    }
  }
  if(!is.na(sf::st_crs(mosaic))){
    suppressWarnings(sf::st_crs(sto) <- sf::st_crs(mosaic))
  } else{
    suppressWarnings(sf::st_crs(sto) <- "+proj=utm +zone=32 +datum=WGS84 +units=m")
  }
  if(viewopt == "rgb"){
    if(terra::nlyr(mosaic) > 2){
      if(!is.na(sf::st_crs(mosaic))){

        message("Using `show = 'rgb' may not produce accurate cropping coordinates.\n Please, consider using `show = 'index'`instead.")
        map <-
          leaflet::leaflet() |>
          leaflet::addScaleBar(position = "bottomleft") |>
          leaflet::addTiles(options = leaflet::providerTileOptions(minZoom = 3, maxZoom = 30)) |>
          leafem::addStarsRGB(sto,
                              r = r,
                              g = g,
                              b = b,
                              quantiles = quantiles,
                              domain = domain,
                              maxBytes = 64 * 1024 * 1024,
                              na.color = "#00000000") |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
      } else{

        map <-
          leaflet::leaflet() |>
          leafem::addStarsRGB(sto,
                              r = r,
                              g = g,
                              b = b,
                              quantiles = quantiles,
                              domain = domain,
                              maxBytes = 64 * 1024 * 1024,
                              na.color = "#00000000") |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
      }
    } else{
      if(vieweropt == "base"){
        terra::plot(terra::rast(sto),
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        index <- gsub("[/\\\\]", "_", index, perl = TRUE)
        map <-
          mapview::mapview(sto,
                           layer.name = index,
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE) |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
        invisible(map)
      }

    }
  } else{
    if(terra::nlyr(mosaic) > 2){
      index <- gsub("[/\\\\]", "_", index, perl = TRUE)
      if(vieweropt == "base"){
        terra::plot(terra::rast(sto),
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(sto,
                           layer.name = index,
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE) |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
        invisible(map)
      }
    } else{
      if(vieweropt == "base"){
        terra::plot(terra::rast(sto),
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(sto,
                           layer.name = names(mosaic),
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE) |>
          mapedit::editMap(editor = "leafpm",
                           title = title)
        invisible(map)
      }
    }
  }
}


#' Create and Export mosaics
#' @details
#' * `mosaic_input()` is a simply wrapper around [terra::rast()]. It creates a
#' `SpatRaster` object from scratch, from a filename, or from another object.
#' * `mosaic_export()` is a simply wrapper around [terra::writeRaster()]. It write
#' a `SpatRaster` object to a file.
#'
#' @name mosaic_input
#' @param mosaic
#'  * For `mosaic_input()`, a file path to the raster to imported, a matrix,
#'    array or a list of `SpatRaster` objects.
#'  * For `mosaic_export()`, an `SpatRaster` object.
#' @param filename character. The Output filename.
#' @param overwrite logical. If `TRUE`, filename is overwritten.
#' @param ... Additional arguments passed to [terra::rast()] (`mosaic_input()`)
#'   or  [terra::writeRaster()] (`mosaic_output()`)
#'
#' @return
#' * `mosaic_input()` returns an `SpatRaster` object.
#' * `mosaic_export()` do not return an object.
#' @export
#' @examples
#' library(pliman)
#'
#' # create an SpatRaster object based on a matrix
#' x <- matrix(1:20, nrow = 4, ncol = 5)
#' rast <- mosaic_input(x)
#' mosaic_view(rast, viewer = "base", axes = TRUE)
#'
#' # create a temporary filename for the example
#' f <- file.path(tempdir(), "test.tif")
#' mosaic_export(rast, f, overwrite=TRUE)
#' list.files(tempdir())
#'
mosaic_input <- function(mosaic, ...){
  terra::rast(mosaic, ...)
}
#' @export
#' @name mosaic_input
mosaic_export <- function(mosaic,
                          filename,
                          overwrite = FALSE,
                          ...){
  terra::writeRaster(mosaic,
                     filename = filename,
                     overwrite = overwrite,
                     ...)
}



#' Crop a mosaic
#'
#' Crop a `SpatRaster` object based on user-defined selection using an
#' interactive map or plot.
#'
#' @details This function uses the `mosaic_view` function to display an
#'   interactive map or plot of the mosaic raster, allowing users to draw a
#'   rectangle to select the cropping area. The selected area is then cropped
#'   from the input mosaic and returned as a new `SpatRaster` object.
#'
#' @inheritParams mosaic_view
#' @param ... Additional arguments passed to [mosaic_view()].
#'
#' @return A cropped version of `mosaic` based on the user-defined selection.
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' # Load a raster showing the elevation of Luxembourg
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # Generate an interactive map using 'mapview' (works only in an interactive section)
#' cropped <- mosaic_crop(mosaic)
#' mosaic_view(cropped)
#' }
#'
mosaic_crop <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        re = 4,
                        nir = 5,
                        show = c("rgb", "index"),
                        index = "R",
                        max_pixels = 500000,
                        downsample = NULL,
                        ...){
  showopt <- c("rgb", "index")
  showopt <- showopt[pmatch(show[[1]], showopt)]
  controls <- mosaic_view(mosaic,
                          show = showopt,
                          index = index,
                          r = r,
                          g = g,
                          b = b,
                          nir = nir,
                          re = re,
                          max_pixels = max_pixels,
                          downsample = downsample,
                          title = "Use the 'Draw rectangle' tool to select the cropping area.",
                          ...)
  if(!is.na(sf::st_crs(mosaic))){
    grids <-
      sf::st_make_grid(controls$finished, n = c(1, 1)) |>
      sf::st_transform(sf::st_crs(mosaic))
  } else{
    sf::st_crs(mosaic) <- sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
    grids <-
      sf::st_make_grid(controls$finished, n = c(1, 1)) |>
      sf::st_transform(sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m"))
  }
  cropped <- terra::crop(mosaic, grids)
  invisible(cropped)
}


#' Mosaic Index
#'
#' Compute or extract an index layer from a multi-band mosaic raster.
#' @inheritParams mosaic_view
#' @inheritParams image_index
#'
#' @return An index layer extracted/computed from the mosaic raster.
#'
#' @details This function computes or extracts an index layer from the input
#'   mosaic raster based on the specified index name. If the index is not found
#'   in the package's predefined index list (see [image_index()] for more
#'   details), it attempts to compute the index using the specified band
#'   indices. The resulting index layer is returned as an `SpatRaster` object.
#' @export
#' @examples
#' library(pliman)
#' library(terra)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' confusion <-
#'      matrix(rnorm(90*95, 100, 30),
#'             nrow = nrow(mosaic),
#'             ncol = ncol(mosaic))
#' confusion <- mosaic_input(confusion)
#' terra::ext(confusion) <- terra::ext(mosaic)
#' terra::crs(confusion) <- terra::crs(mosaic)
#' names(confusion) <- "confusion"
#' two_layers <- c(mosaic, confusion)
#' final <- mosaic_index(two_layers, "B+R", b = 1, r = 2)
#' mosaic_view(mosaic_input(final), viewer = "base")
#'
mosaic_index <- function(mosaic,
                         index = "R",
                         r = 3,
                         g = 2,
                         b = 1,
                         re = 4,
                         nir = 5){
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




#' Mosaic to pliman
#'
#' Convert an `SpatRaster` object to a `Image` object with optional scaling.
#' @inheritParams mosaic_view
#' @param rescale Rescale the final values? If `TRUE` the final values are
#'   rescaled so that the maximum value is 1.
#' @param coef An addition coefficient applied to the resulting object. This is
#'   useful to adjust the brightness of the final image. Defaults to 0.
#'
#' @return An `Image` object with the same number of layers as `mosaic`.
#'
#' @details This function converts `SpatRaster` into an `Image` object, which
#'   can be used for image analysis in `pliman`. Note that if a large
#'   `SpatRaster` is loaded, the resulting object may increase considerably the
#'   memory usage.
#' @importFrom stars st_dimensions
#' @export
#' @examples
#' library(pliman)
#' # Convert a mosaic raster to an Image object
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' pliman_image <- mosaic_to_pliman(mosaic)
#' plot(pliman_image)
#'
mosaic_to_pliman <- function(mosaic,
                             r = 3,
                             g = 2,
                             b = 1,
                             re = 4,
                             nir = 5,
                             rescale =  TRUE,
                             coef = 0){
  if(class(mosaic) %in% c("RasterStack","RasterLayer","RasterBrick")){
    mosaic <- terra::rast(mosaic)
  }
  nlr <- terra::nlyr(mosaic)
  mosaic <- stars::st_as_stars(mosaic, proxy = FALSE)
  mosaic[mosaic == 65535] <- NA
  if(nlr == 5){
    mosaic <- EBImage::Image(mosaic[[1]])[,, c(r, g, b, re, nir)]
  } else if(nlr == 3){
    mosaic <- EBImage::Image(mosaic[[1]])[,, c(r, g, b)]
  } else{
    mosaic <- EBImage::Image(mosaic[[1]])
  }
  if(isTRUE(rescale)){
    mosaic <- mosaic / max(mosaic, na.rm = TRUE)
  }
  if(nlr == 3){
    EBImage::colorMode(mosaic) <- "color"
  }
  return(mosaic + coef)
}


#' Mosaic to RGB
#'
#' Convert an `SpatRaster` to a three-band RGB image of class `Image`.
#'
#' @inheritParams mosaic_to_pliman
#' @param plot Logical, whether to display the resulting RGB image (default:
#'   TRUE).
#'
#' @return A three-band RGB image represented as a pliman (EBImage) object.
#'
#' @details This function converts `SpatRaster` that contains the RGB bands into
#'   a three-band RGB image using pliman (EBImage). It allows you to specify the
#'   band indices for the red, green, and blue channels, as well as apply a
#'   scaling coefficient to the final image. By default, the resulting RGB image
#'   is displayed, but this behavior can be controlled using the `plot`
#'   parameter.
#'
#' @export
#' @examples
#' library(pliman)
#' # Convert a mosaic raster to an RGB image and display it
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # Convert a mosaic raster to an RGB image without displaying it
#' rgb_image <- mosaic_to_rgb(c(mosaic * 2, mosaic - 0.3, mosaic * 0.8))
#' plot(rgb_image)
#'
#'
mosaic_to_rgb <- function(mosaic,
                          r = 3,
                          g = 2,
                          b = 1,
                          coef = 0,
                          plot = TRUE){
  ebim <- mosaic_to_pliman(mosaic,
                           r = r,
                           g = g,
                           b = b,
                           coef = coef)[,,c(r, g, b)]
  EBImage::colorMode(ebim) <- "color"
  invisible(ebim)
}


#' Prepare a mosaic
#'
#' Prepare an `SpatRaster` object to be analyzed in pliman. This includes
#' cropping the original mosaic, aligning it, and cropping the aligned object.
#' The resulting object is an object of class `Image` that can be further
#' analyzed.
#' @inheritParams mosaic_view
#' @inheritParams mosaic_to_pliman
#' @param crop_mosaic Logical, whether to crop the mosaic interactively before
#'   aligning it (default: FALSE).
#' @param align Logical, whether to align the mosaic interactively (default:
#'   TRUE).
#' @param crop_aligned Logical, whether to crop the aligned mosaic interactively
#'   (default: TRUE).
#'
#' @return A prepared object of class `Image`.
#'
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_prepare(mosaic)
#' }
#'
mosaic_prepare <- function(mosaic,
                           r = 3,
                           g = 2,
                           b = 1,
                           re = 4,
                           nir = 5,
                           crop_mosaic = FALSE,
                           align = TRUE,
                           crop_aligned = TRUE,
                           rescale =  TRUE,
                           coef = 0,
                           viewer = "mapview",
                           max_pixels = 500000,
                           show = "rgb",
                           index = "R"){
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[1], vieweropt)]
  if(isTRUE(crop_mosaic)){
    cropped <- mosaic_crop(mosaic,
                           show = show,
                           index = index,
                           max_pixels = max_pixels,
                           r = r,
                           g = g,
                           b = b,
                           nir = nir,
                           re = re)
    ebimg <- mosaic_to_pliman(cropped,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              rescale = rescale,
                              coef = coef)
    if(vieweropt != "base"){
      image_view(ebimg[1:5, 1:5,])
    }
  } else{
    ebimg <- mosaic_to_pliman(mosaic,
                              r = r,
                              g = g,
                              b = b,
                              re = re,
                              nir = nir,
                              rescale = rescale,
                              coef = coef)
  }
  if(isTRUE(align)){
    aligned <- image_align(ebimg, viewer = viewer)
    if(vieweropt != "base"){
      image_view(aligned[1:5, 1:5,])
    }
  } else{
    aligned <- ebimg
  }

  if(isTRUE(crop_aligned)){
    cropped <- image_crop(aligned, viewer = vieweropt)
  } else{
    cropped <- aligned
  }
  if(dim(cropped)[3] == 3){
    EBImage::colorMode(cropped) <- "color"
  }
  invisible(cropped)
}
