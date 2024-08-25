
add_width_height <- function(grid, width, height, mosaic, points_align) {
  gridl <-lapply(sf::st_geometry(grid), sf::st_coordinates)
  gridadj <- add_width_height_cpp(gridl, height, width, points_align)
  grd <- lapply(gridadj, function(x){sf::st_polygon(list(x))}) |> sf::st_sfc()
  grd <-  sf::st_sf(geometry = grd)
  sf::st_crs(grd) <- sf::st_crs(grid)
  return(grd)
}
create_buffer <- function(coords, buffer_col, buffer_row) {
  # Calculate the new x-min, x-max, y-min, and y-max after adjustment
  coords <- sf::st_coordinates(coords)
  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  y_min <- min(coords[, 2])
  y_max <- max(coords[, 2])
  new_x_min <- x_min - buffer_col * (x_max - x_min)
  new_x_max <- x_max + buffer_col * (x_max - x_min)
  new_y_min <- y_min - buffer_row * (y_max - y_min)
  new_y_max <- y_max + buffer_row * (y_max - y_min)

  # Calculate the scaling factors for x and y
  x_scale_factor <- (new_x_max - new_x_min) / (x_max - x_min)
  y_scale_factor <- (new_y_max - new_y_min) / (y_max - y_min)

  # Apply the scaling to the coordinates
  resized_coords <- coords
  resized_coords[, 1] <- (resized_coords[, 1] - x_min) * x_scale_factor + new_x_min
  resized_coords[, 2] <- (resized_coords[, 2] - y_min) * y_scale_factor + new_y_min
  sf::st_polygon(list(resized_coords[, 1:2]))
  # return(resized_coords)
  # return(data.frame(resized_coords) |> sf::st_as_sf(coords = c("X", "Y")))
}

make_grid <- function(points, nrow, ncol, mosaic, buffer_col = 0, buffer_row = 0, plot_width = NULL, plot_height = NULL) {
  points_align <-
    sf::st_transform(points, sf::st_crs(mosaic)) |>
    sf::st_coordinates()

  grids <-
    sf::st_make_grid(points, n = c(nrow, ncol)) |>
    sf::st_transform(sf::st_crs(mosaic))

  sxy <-
    points |>
    sf::st_make_grid(n = c(1, 1)) |>
    sf::st_cast("POINT") |>
    rev() |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_coordinates()
  txy <-
    points |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_coordinates()
  txy <- txy[1:4, 1:2]

  cvm <- lm(txy ~ sxy[1:4, ])
  parms <- cvm$coefficients[2:3, ]
  intercept <- cvm$coefficients[1, ]
  geometry <- grids * parms + intercept

  if (buffer_row != 0 | buffer_col != 0) {
    geometry <- lapply(geometry, function(g) {
      create_buffer(g, buffer_col, buffer_row)
    }) |>
      sf::st_sfc(crs = sf::st_crs(mosaic))
  }
  if (!is.null(plot_width) & !is.null(plot_height)) {
    geometry <-
      add_width_height(grid = geometry, width = plot_width, height = plot_height, points_align = points_align[2:3, 1:2]) |>
      sf::st_as_sfc(crs = sf::st_crs(mosaic))
  }
  return(sf::st_sf(geometry = geometry, crs = sf::st_crs(mosaic)))
}

#' Generate plot IDs with different layouts
#'
#' Based on a shapefile, number of columns and rows, generate plot IDs with
#' different layouts.
#'
#' @param shapefile An object computed with [shapefile_build()]
#' @param nrow The number of columns
#' @param ncol The number of rows
#' @param layout Character: one of
#'  * `'tblr'` for top/bottom left/right orientation
#'  * `'tbrl'` for top/bottom right/left orientation
#'  * `'btlr'` for bottom/top left/right orientation
#'  * `'btrl'` for bottom/top right/left orientation
#'  * `'lrtb'` for left/right top/bottom orientation
#'  * `'lrbt'` for left/right bottom/top orientation
#'  * `'rltb'` for right/left top/bottom orientation
#'  * `'rlbt'` for right/left bottom/top orientation
#' @param plot_prefix The plot_id prefix. Defaults to `'P'`.
#' @param serpentine Create a serpentine-based layout? Defaults to `FALSE`.
#' @return A vector of plot IDs with specified layout
#' @export
#'
plot_id <- function(shapefile,
                    nrow,
                    ncol,
                    layout = c("tblr", "tbrl", "btlr", "btrl", "lrtb", "lrbt", "rltb", "rlbt"),
                    plot_prefix = "P",
                    serpentine = FALSE) {
  # Ensure the specified layout is valid
  allowed <- c("tblr", "tbrl", "btlr", "btrl", "lrtb", "lrbt", "rltb", "rlbt")
  layout <- layout[[1]]
  if (!layout %in% allowed) {
    stop(paste0("`layout` must be one of the following: ", paste0(allowed, collapse = ", ")))
  }

  # Ensure that the number of rows in the shapefile matches expected dimensions
  expected_rows <- nrow * ncol
  if (nrow(shapefile) != expected_rows) {
    stop(paste("Expected", expected_rows, "rows, but shapefile has", nrow(shapefile), "rows."))
  }

  # Helper function for generating plot names
  leading_zeros <- function(x, n) {
    sprintf(paste0("%0", n, "d"), x)
  }

  plots_tblr <- paste0(plot_prefix, leading_zeros(1:nrow(shapefile), 4))

  # Define layout functions
  make_tblr <- function() {
    plots_tblr
  }

  make_tbrl <- function() {
    plots_tblr_rev <- rev(plots_tblr)
    plots_tbrl <- NULL
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_tbrl <- c(plots_tbrl, rev(plots_tblr_rev[start:end]))
    }
    plots_tbrl
  }

  make_btrl<- function() {
    plots_rev <- rev(plots_tblr)
    plots_btlr <- NULL
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_btlr <- c(plots_btlr, plots_rev[start:end])
    }
    plots_btlr
  }

  make_btlr <- function() {
    plots_btlr_rev <- rev(make_btrl())
    plots_btrl <- NULL
    for (i in seq_len(ncol)) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_btrl <- c(plots_btrl, rev(plots_btlr_rev[start:end]))
    }
    plots_btrl
  }

  make_lrtb <- function() {
    plots_lrtb <- NULL
    for (i in 1:ncol) {
      plots_lrtb <- c(plots_lrtb, plots_tblr[seq(i, length(plots_tblr), by = ncol)])
    }
    plots_lrtb
  }

  make_lrbt <- function() {
    plots_lrbt <- NULL
    plots_lrtb <- make_lrtb()
    for (i in 1:ncol) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_lrbt <- c(plots_lrbt, rev(plots_lrtb[start:end]))
    }
    plots_lrbt
  }

  make_rltb <- function() {
    plots_rltb <- NULL
    for (i in 1:ncol) {
      # Columns from right to left
      plots_rltb <- c(plots_rltb, plots_tblr[seq(ncol - i + 1, length(plots_tblr), by = ncol)])
    }
    plots_rltb
  }

  make_rlbt <- function() {
    plots_rltb <- make_rltb()
    plots_rlbt <- NULL
    for (i in seq_len(ncol)) {
      start <- (i - 1) * nrow + 1
      end <- start + nrow - 1
      plots_rlbt <- c(plots_rlbt, rev(plots_rltb[start:end]))
    }
    plots_rlbt
  }

  # Return the appropriate layout
  plots <-  switch(layout,
                   "tblr" = make_tblr(),
                   "tbrl" = make_tbrl(),
                   "btlr" = make_btlr(),
                   "btrl" = make_btrl(),
                   "lrtb" = make_lrtb(),
                   "lrbt" = make_lrbt(),
                   "rltb" = make_rltb(),
                   "rlbt" = make_rlbt())
  mat <- matrix(plots, ncol = ncol, nrow = nrow)
  if(serpentine){
    # column serpentine
    if(layout %in% c("tblr", "btlr")){
      mat2 <- mat
      for (j in 1:ncol(mat)) {
        if(j %% 2 == 0){
          mat2[, j] <- rev(mat[, j])
        } else{
          mat2[, j]
        }
      }
    }
    if(layout %in% c( "tbrl", "btrl")){
      mat2 <- mat
      cols <- ncol(mat):1
      for (j in cols[seq(2, ncol, by = 2)]) {
        mat2[, j] <- rev(mat[, j])
      }
    }
    # row serpentine
    if(layout %in% c("lrtb", "rltb")){
      mat2 <- mat
      for (j in 1:nrow(mat)) {
        if(j %% 2 == 0){
          mat2[j, ] <- rev(mat[j, ])
        } else{
          mat2[j, ]
        }
      }
    }
    if(layout %in% c("lrbt", "rlbt")){
      mat2 <- mat
      rows <- nrow(mat):1
      for (j in rows[seq(2, nrow, by = 2)]) {
        mat2[j, ] <- rev(mat[j, ])
      }
    }
  } else{
    mat2 <- mat
  }
  return(as.vector(mat2))
}

#' Build a shapefile from a mosaic raster
#'
#' This function takes a mosaic raster to create a shapefile containing polygons
#' for the specified regions. Users can drawn Areas of Interest (AOIs) that can
#' be either a polygon with n sides, or a grid, defined by `nrow`, and `ncol`
#' arguments.
#' @details
#' Since multiple blocks can be created, the length of arguments `grid`, `nrow`,
#' `ncol`, `buffer_edge`, `buffer_col`, and `buffer_row` can be either an scalar
#' (the same argument applied to all the drawn blocks), or a vector with the
#' same length as the number of drawn blocks. In the last, shapefiles in each
#' block can be created with different dimensions.
#' @param sf_to_polygon Convert sf geometry like POINTS and LINES to POLYGONS?
#'   Defaults to `FALSE`. Using `TRUE` allows using POINTS to extract values
#'   from a raster using `exactextractr::exact_extract()`.
#' @param basemap An optional `mapview` object.
#' @param controlpoints An `sf` object created with [mapedit::editMap()],
#'   containing the polygon that defines the region of interest to be analyzed.
#' @inheritParams mosaic_analyze
#' @inheritParams mosaic_index
#' @inheritParams mosaic_view
#' @inheritParams utils_shapefile
#' @inheritParams plot_id
#' @return A list with the built shapefile. Each element is an `sf` object with
#'   the coordinates of the drawn polygons.
#' @export
#' @examples
#' if(interactive()){
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' shps <-
#'       shapefile_build(mosaic,
#'                       nrow = 6,
#'                       ncol = 3,
#'                       buffer_row = -0.05,
#'                       buffer_col = -0.25,
#'                       check_shapefile = FALSE,
#'                       build_shapefile = FALSE) ## Use TRUE to interactively build the plots
#' mosaic_plot(mosaic)
#' shapefile_plot(shps[[1]], add = TRUE)
#' }
#'

shapefile_build <- function(mosaic,
                            basemap = NULL,
                            controlpoints = NULL,
                            r = 3,
                            g = 2,
                            b = 1,
                            crop_to_shape_ext = TRUE,
                            grid = TRUE,
                            nrow = 1,
                            ncol = 1,
                            plot_width = NULL,
                            plot_height = NULL,
                            layout = "lrtb",
                            serpentine = TRUE,
                            build_shapefile = TRUE,
                            check_shapefile = FALSE,
                            sf_to_polygon = FALSE,
                            buffer_edge = 1,
                            buffer_col = 0,
                            buffer_row = 0,
                            as_sf = TRUE,
                            verbose = TRUE,
                            max_pixels = 1000000,
                            downsample = NULL,
                            quantiles =  c(0, 1)){
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:4326")
  }
  ress <- terra::res(mosaic)
  nlyrs <- terra::nlyr(mosaic)
  if(build_shapefile){
    if(verbose){
      message("\014","\nBuilding the mosaic...\n")
    }
    if(is.null(basemap)){
      basemap <- mosaic_view(mosaic,
                             r = r,
                             g = g,
                             b = b,
                             max_pixels = max_pixels,
                             downsample = downsample,
                             quantiles = quantiles)
    }
    if(is.null(controlpoints)){
      points <- mapedit::editMap(basemap, editor = "leafpm")
      cpoints <- points$finished
    } else{
      cpoints <- controlpoints
    }
    if(sf_to_polygon){
      cpoints <- cpoints |> sf_to_polygon()
    }
  } else{
    extm <- terra::ext(mosaic)
    xmin <- extm[1]
    xmax <- extm[2]
    ymin <- extm[3]
    ymax <- extm[4]
    coords <- matrix(c(xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin, xmin, ymax), ncol = 2, byrow = TRUE)
    # Create a Polygon object
    polygon <- sf::st_polygon(list(coords))
    # Create an sf object with a data frame that includes the 'geometry' column
    if(sum(ress) == 2){
      crop_to_shape_ext <- FALSE

    }
    cpoints <- sf::st_sf(data.frame(id = 1),
                         geometry = sf::st_sfc(polygon),
                         crs = sf::st_crs(mosaic))
  }

  # crop to the analyzed area
  if(crop_to_shape_ext){
    if(sum(ress) != 2){
      cpoints <- cpoints |> sf::st_transform(crs = sf::st_crs(terra::crs(mosaic)))
    }
    poly_ext <-
      cpoints |>
      terra::vect() |>
      terra::buffer(buffer_edge) |>
      terra::ext()
    mosaiccr <- terra::crop(mosaic, poly_ext)

  } else{
    mosaiccr <- mosaic
  }
  # check the parameters
  nrow <- validate_and_replicate2(nrow, cpoints)
  ncol <- validate_and_replicate2(ncol, cpoints)
  layout <- validate_and_replicate2(layout, cpoints)
  buffer_col <- validate_and_replicate2(buffer_col, cpoints)
  buffer_row <- validate_and_replicate2(buffer_row, cpoints)
  plot_width <- validate_and_replicate2(plot_width, cpoints)
  plot_height <- validate_and_replicate2(plot_height, cpoints)
  serpentine <- validate_and_replicate2(serpentine, cpoints)
  grid <- validate_and_replicate2(grid, cpoints)

  # check the created shapes?
  if(verbose){
    message("\014","\nCreating the shapes...\n")
  }
  created_shapes <- list()
  for(k in 1:nrow(cpoints)){
    if(inherits(cpoints[k, ]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(cpoints[k, ])) == 5 & grid[[k]]){
      pg <-
        make_grid(cpoints[k, ],
                  nrow = nrow[k],
                  ncol = ncol[k],
                  mosaic = mosaic,
                  buffer_col = buffer_col[k],
                  buffer_row = buffer_row[k],
                  plot_width = plot_width[k],
                  plot_height = plot_height[k]) |>
        dplyr::mutate(row = rep(1:nrow[k], ncol[k]),
                      column = rep(1:ncol[k], each = nrow[k]),
                      .before = geometry)
      pg <-
        pg |>
        dplyr::mutate(unique_id = dplyr::row_number(),
                      block = paste0("B", leading_zeros(k, 2)),
                      plot_id = plot_id(pg, nrow = nrow[k], ncol = ncol[k], layout = layout[k], serpentine = serpentine[k]),
                      .before = 1)
    } else{
      pg <-
        cpoints[k, ] |>
        sf::st_transform(sf::st_crs(mosaic)) |>
        dplyr::select(geometry) |>
        dplyr::mutate(unique_id = dplyr::row_number(),
                      block = paste0("B", leading_zeros(k, 2)),
                      plot_id = "P0001",
                      row = 1,
                      column = 1,
                      .before = 1)
    }

    created_shapes[[k]] <- pg

  }
  if(check_shapefile){
    if(verbose){
      message("\014","\nChecking the built shapefile...\n")
    }
    lengths <- sapply(created_shapes, nrow)
    pg_edit <-
      do.call(rbind, lapply(seq_along(created_shapes), function(i){
        created_shapes[[i]] |>
          dplyr::mutate(`_leaflet_id` = 1:nrow(created_shapes[[i]]),
                        feature_type = "polygon") |>
          dplyr::relocate(geometry, .after = 4) |>
          sf::st_transform(crs = 4326)
      }))
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- mosaic_aggregate(mosaiccr, pct = round(100 / downsample))
    }
    if(build_shapefile){
      mapview::mapview() |> mapedit::editMap()
    }
    edited <-
      mapedit::editFeatures(pg_edit, basemap) |>
      dplyr::select(geometry, block, plot_id, row, column) |>
      dplyr::mutate(unique_id = dplyr::row_number(), .before = 1) |>
      sf::st_transform(sf::st_crs(mosaic))
    sfeat <- sf::st_as_sf(edited)
    sf::st_geometry(sfeat) <- "geometry"
    created_shapes <- split(sfeat, edited$block)
  }
  if(verbose){
    message("\014","\nShapefile finished...\n")
  }
  if(!as_sf){
    return(shapefile_input(created_shapes, info = FALSE, as_sf = FALSE))
  } else{
    return(created_shapes)
  }
}


#' A wrapper around terra::plot()
#'
#' Plot the values of a SpatVector
#'
#' @param shapefile An SpatVector of sf object.
#' @param ... Further arguments passed on to [terra::plot()].
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' library(pliman)
#' r <- shapefile_input(system.file("ex/lux.shp", package="terra"))
#' shapefile_plot(r)
shapefile_plot <- function(shapefile, ...){
  if(!inherits(shapefile, "SpatVector") & !inherits(shapefile, "sf") ){
    stop("'mosaic' must be an object of class 'SpatVector' of 'sf'")
  }
  if(inherits(shapefile, "sf")){
    shapefile <- terra::vect(shapefile)
  }
  terra::plot(shapefile, ...)
}


#' Import/export shapefiles.
#' @description
#'
#' * `shapefile_input()` creates or importes a shapefile and optionally convert
#'  it to an `sf` object.
#' * `shapefile_export()` exports an object (`sf` or `SpatVector`) to a file
#' * `shapefile_view()` is a simple wrapper around mapview() to plot a shapefile.
#' @name utils_shapefile
#'
#' @param shapefile
#'
#'   For `shapefile_input()`, character (filename),  or an object that can be
#'   coerced to a SpatVector, such as an `sf` (simple features) object. See
#'   [terra::vect()] for more details.
#'
#'   For `shapefile_export()`, `SpatVector` or an `sf` object to be exported as
#'   a shapefile.
#'
#' @param info Logical value indicating whether to print information about the
#'   imported shapefile (default is `TRUE`).
#' @param as_sf Logical value indicating whether to convert the imported
#'   shapefile to an `sf` object (default is `TRUE`).
#' @param filename The path to the output shapefile.
#' @param attribute The attribute to be shown in the color key. It must be a
#'   variable present in `shapefile`.
#' @param color_regions The color palette to represent `attribute`.
#' @param ... Additional arguments to be passed to [terra::vect()]
#'   (`shapefile_input()`), [terra::writeVector()] (`shapefile_export()`) or
#'   [mapview::mapview()] (`shapefile_view()`).
#'
#' @return
#'  * `shapefile_input()` returns an object of class `sf` (default) representing
#'  the imported shapefile.
#'
#'  * `shapefile_export()` returns a `NULL` object.
#'
#'  * `shapefile_view()` returns an object of class `mapview`.
#'
#' @examples
#' library(pliman)
#' shp <- system.file("ex/lux.shp", package="terra")
#' shp_file <- shapefile_input(shp, as_sf = FALSE)
#' shapefile_view(shp_file)
#'
#' @export
shapefile_input <- function(shapefile,
                            info = TRUE,
                            as_sf = TRUE,
                            ...) {
  create_shp <- function(shapefile, info, as_sf, ...){
    shp <- terra::vect(shapefile, ...)
    if (terra::crs(shp) == "") {
      message("Missing Coordinate Reference System. Setting to EPSG:3857")
      terra::crs(shp) <- terra::crs("EPSG:3857")
    }
    if (info) {
      print(shp)
    }
    if (as_sf) {
      shp <- sf::st_as_sf(shp)
    }
    return(shp)
  }

  if(inherits(shapefile, "list")){
    shapes <- do.call(rbind, lapply(shapefile, function(x){x}))
    create_shp(shapes, info, as_sf, ...)
  } else{
    create_shp(shapefile, info, as_sf, ...)
  }
}
#' @name utils_shapefile
#' @export
shapefile_export <- function(shapefile, filename, ...) {
  if (!inherits(shapefile, "SpatVector")) {
    shapefile <- try(terra::vect(shapefile))
  }
  if (inherits(shapefile, "try-error")) {
    stop("Unable to coerce the input shapefile to a SpatVector object")
  }
  terra::writeVector(shapefile, filename, ...)
}

#' @name utils_shapefile
#' @export
shapefile_view <- function(shapefile,
                           attribute = NULL,
                           color_regions = custom_palette(c("red", "yellow", "forestgreen")),
                           ...){
  suppressWarnings(
    mapview::mapview(shapefile,
                     zcol = attribute,
                     col.regions = color_regions,
                     ...)
  )
}

#' Edit Features in a Shapefile
#'
#' This function allows you to interactively edit features in a shapefile using
#' the mapedit package.
#'
#' @param shapefile A shapefile (`sf` object) that can be created with
#'   [shapefile_input()].
#' @param mosaic Optionally, a mosaic (SpatRaster) to be displayed as a
#'   background.
#' @param basemap An optional `mapview` object.
#' @param r Red band index for RGB display (default is 3).
#' @param g Green band index for RGB display (default is 2).
#' @param b Blue band index for RGB display (default is 1).
#' @param max_pixels Maximum number of pixels for down-sampling the mosaic
#'   (default is 3e6).
#' @export
#' @return A modified shapefile with user-edited features.
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' shp <- shapefile_input(system.file("ex/lux.shp", package="terra"))
#' edited <- shapefile_edit(shp)
#' }
shapefile_edit <- function(shapefile,
                           mosaic = NULL,
                           basemap = NULL,
                           r = 3,
                           g = 2,
                           b = 1,
                           max_pixels = 3e6){
  shapefile <- shapefile_input(shapefile, info = FALSE)
  if(!is.null(mosaic)){
    if(is.null(basemap)){
      downsample <- find_aggrfact(mosaic, max_pixels = max_pixels)
      if(downsample > 0){
        mosaic <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
      }
      nlyrs <- terra::nlyr(mosaic)
      if(nlyrs > 2){
        map <-
          mapview::viewRGB(
            x = as(mosaic, "Raster"),
            layer.name = "base",
            r = r,
            g = g,
            b = b,
            na.color = "#00000000",
            maxpixels = 5e6
          )
      } else{
        map <-
          mapview::mapview() %>%
          leafem::addGeoRaster(x = as(mosaic[[1]], "Raster"),
                               colorOptions = leafem::colorOptions(palette = custom_palette(),
                                                                   na.color = "transparent"))
      }
    } else{
      map <- basemap
    }
    edited <- mapedit::editFeatures(shapefile |> sf::st_transform(crs = 4326), map)
  } else{
    edited <- mapedit::editFeatures(shapefile)
  }
  return(edited)
}
