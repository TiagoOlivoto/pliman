# gathered from https://github.com/r-spatial/leafem/blob/6d6831352038b8f7462ff7afa698050c4e46fb5e/R/addRasterRGB.R#L143

rscl = function(x,
                from = range(x, na.rm = TRUE, finite = TRUE),
                to = c(0, 1),
                ...) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}
add_rgb <- function(
    map,
    x,
    r = 3, g = 2, b = 1,
    group = NULL,
    layerId = NULL,
    resolution = 96,
    opacity = 0.8,
    options = leaflet::tileOptions(),
    colorOptions = NULL,
    project = TRUE,
    pixelValuesToColorFn = NULL,
    ...
) {

  if (inherits(x, "Raster")) {
    x = stars::st_as_stars(x)
  }

  if (project & !sf::st_is_longlat(x)) {
    x = stars::st_warp(x, crs = 4326)
  }

  if (is.null(colorOptions)) {
    colorOptions = leafem::colorOptions()
  }

  fl = tempfile(fileext = ".tif")

  if (inherits(x, "stars_proxy")) {
    # file.copy(x[[1]], fl)
    fl = x[[1]]
  }

  if (!inherits(x, "stars_proxy")) {
    stars::write_stars(x, dsn = fl)
  }

  minband = min(r, g, b)

  rgbPixelfun = htmlwidgets::JS(
    sprintf(
      "
        pixelValuesToColorFn = values => {
        // debugger;
          if (isNaN(values[0])) return '%s';
          return rgbToHex(
            Math.ceil(values[%s])
            , Math.ceil(values[%s])
            , Math.ceil(values[%s])
          );
        };
      "
      , colorOptions[["naColor"]]
      , r - minband
      , g - minband
      , b - minband
    )
  )

  # todo: streching via quantiles and domain...

  leafem::addGeotiff(
    map
    , file = fl
    , url = NULL
    , group = group
    , layerId = layerId
    , resolution = resolution
    , bands = c(r, g, b)
    , arith = NULL
    , opacity = opacity
    , options = options
    , colorOptions = colorOptions
    , rgb = TRUE
    , pixelValuesToColorFn = rgbPixelfun
  )
}

create_buffer <- function(coords, buffer_col, buffer_row) {
  # Calculate the new x-min, x-max, y-min, and y-max after adjustment
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
  return(resized_coords)
}
make_grid <- function(points,
                      nrow,
                      ncol,
                      buffer_col = 0,
                      buffer_row = 0){
  grids <-
    sf::st_make_grid(points, n = c(nrow, ncol)) |>
    sf::st_transform(sf::st_crs(mosaic))

  for(k in 1:length(grids)){
    coords <- grids[[k]][[1]]
    grids[[k]][[1]] <- create_buffer(coords,
                                     buffer_col = buffer_row,
                                     buffer_row = buffer_col)
  }

  sxy <-
    points |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_make_grid(n = c(1, 1)) |>
    sf::st_cast("POINT") |>
    rev() |>
    sf::st_coordinates()

  txy <-
    points |>
    sf::st_transform(sf::st_crs(mosaic)) |>
    sf::st_coordinates()
  txy <- txy[1:4, 1:2]
  cvm <- lm(formula = txy ~ sxy[1:4, ])
  parms <- cvm$coefficients[2:3, ]
  intercept <- cvm$coefficients[1, ]
  geometry <- grids * parms + intercept
  gshp <-
    geometry |>
    sf::st_sf(crs = sf::st_crs(mosaic))
  return(gshp)
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


compute_measures_mosaic <- function(contour){
  lw <- help_lw(contour)
  cdist <- help_centdist(contour)
  data.frame(area = help_area(contour),
             perimeter = sum(help_distpts(contour)),
             length = lw[[1]],
             width = lw[[2]],
             diam_min = min(cdist),
             diam_mean = mean(cdist),
             diam_max = max(cdist))

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
#'
#' @inheritParams mosaic_analyze
#' @inheritParams utils_shapefile
#' @export
#'
shapefile_build <- function(mosaic,
                            r = 1,
                            g = 2,
                            b = 3,
                            re = 4,
                            nir = 5,
                            grid = TRUE,
                            nrow = 1,
                            ncol = 1,
                            build_shapefile = TRUE,
                            check_shapefile = FALSE,
                            buffer_edge = 5,
                            buffer_col = 0,
                            buffer_row = 0,
                            as_sf = TRUE,
                            verbose = TRUE,
                            max_pixels = 1000000,
                            downsample = NULL,
                            quantiles =  c(0, 1)){
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  nlyrs <- terra::nlyr(mosaic)
  if(verbose){
    cat("\014","\nBuilding the mosaic...\n")
  }
  if(build_shapefile){
    points <- mosaic_view(mosaic,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          max_pixels = max_pixels,
                          downsample = downsample,
                          quantiles = quantiles,
                          edit = TRUE)
    cpoints <- points$finished
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
    cpoints <- sf::st_sf(data.frame(id = 1),
                         geometry = sf::st_sfc(polygon),
                         crs = sf::st_crs(mosaic))
  }

  # crop to the analyzed area
  poly_ext <-
    cpoints |>
    sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
    terra::vect() |>
    terra::buffer(buffer_edge) |>
    terra::ext()
  mosaiccr <- terra::crop(mosaic, poly_ext)
  mosaiccr[mosaiccr == 65535] <- NA

  # check the parameters
  if(length(nrow) == 1 & nrow(cpoints) != 1){
    nrow <- rep(nrow, nrow(cpoints))
  }
  if(length(nrow) != nrow(cpoints)){
    warning(paste0("`nrow` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
  }
  if(length(ncol) == 1 & nrow(cpoints) != 1){
    ncol <- rep(ncol, nrow(cpoints))
  }
  if(length(ncol) != nrow(cpoints)){
    warning(paste0("`ncol` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
  }
  if(length(buffer_col) == 1 & nrow(cpoints) != 1){
    buffer_col <- rep(buffer_col, nrow(cpoints))
  }
  if(length(buffer_col) != nrow(cpoints)){
    warning(paste0("`buffer_col` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
  }
  if(length(buffer_row) == 1 & nrow(cpoints) != 1){
    buffer_row <- rep(buffer_row, nrow(cpoints))
  }
  if(length(buffer_row) != nrow(cpoints)){
    warning(paste0("`buffer_row` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
  }
  if(length(grid) == 1 & nrow(cpoints) != 1){
    grid <- rep(grid, nrow(cpoints))
  }
  if(length(grid) != nrow(cpoints)){
    warning(paste0("`grid` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
  }

  # check the created shapes?
  if(verbose){
    cat("\014","\nCreating the shapes...\n")
  }
  # cpoints <- re$shapefile$finished
  created_shapes <- list()
  for(k in 1:nrow(cpoints)){
    if(inherits(cpoints[k, ]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(cpoints[k, ])) == 5 & grid[[k]]){
      plot_grid <-
        make_grid(cpoints[k, ],
                  nrow = nrow[k],
                  ncol = ncol[k],
                  buffer_col = buffer_col[k],
                  buffer_row = buffer_row[k])
    } else{
      plot_grid <-  cpoints[k, ] |> sf::st_transform(sf::st_crs(mosaic)) |> poorman::select(geometry)
    }
    created_shapes[[k]] <- plot_grid

  }
  if(check_shapefile){
    if(verbose){
      cat("\014","\nChecking the build shapefile...\n")
    }
    lengths <- sapply(created_shapes, nrow)
    pg_edit <-
      do.call(rbind, lapply(seq_along(created_shapes), function(i){
        created_shapes[[i]] |>
          poorman::mutate(`_leaflet_id` = 1:nrow(created_shapes[[i]]),
                          feature_type = "polygon",
                          block = i) |>
          poorman::relocate(geometry, .after = 3) |>
          sf::st_transform(crs = 4326)
      }))

    if(nlyrs > 2){
      map <-
        mapview::mapview() %>%
        add_rgb(x = as(mosaiccr, "Raster"),
                r = r,
                g = g,
                b = b,
                na.color = "#00000000",
                quantiles = quantiles)
    } else{
      map <-
        mapview::mapview() %>%
        leafem::addGeoRaster(x = as(mosaic[[1]], "Raster"),
                             colorOptions = leafem::colorOptions(palette = custom_palette(),
                                                                  na.color = "transparent"))
    }
    if(build_shapefile){
      mapview::mapview() |> mapedit::editMap()
    }
    edited <-
      mapedit::editFeatures(pg_edit, map) |>
      poorman::select(geometry, block) |>
      sf::st_transform(sf::st_crs(mosaic))
    sfeat <- sf::st_as_sf(edited$geometry)
    sf::st_geometry(sfeat) <- "geometry"
    created_shapes <- split(sfeat, edited$block)
  }
  if(verbose){
    cat("\014","\nShapefile finished...\n")
  }
  if(!as_sf){
    return(shapefile_input(created_shapes, info = FALSE, as_sf = FALSE))
  } else{
    return(created_shapes)
  }
}


#' Analyze a mosaic of remote sensing data
#'
#' This function analyzes a mosaic of remote sensing data, extracting information
#' from specified regions of interest (ROIs) defined in a shapefile or interactively
#' drawn on the mosaic. It computes various vegetation indices and statistical
#' summaries for segmentation-based analysis of the ROIs.
#'
#' @details
#' Since multiple blocks can be analyzed, the length of arguments `grid`,
#' `nrow`, `ncol`, `buffer_edge`, , `buffer_col`, `buffer_row`, `segment_plot`,
#' `segment_individuals`, `includ_if`, `threshold`, `segment_index`, `invert`,
#' and `filter`, can be either an scalar (the same argument applied to all the
#' drawn blocks), or a vector with the same length as the number of drawn. In
#' the last, each block can be analyzed with different arguments.
#'
#' @inheritParams mosaic_view
#' @param grid Logical, indicating whether to use a grid for segmentation
#'   (default: TRUE).
#' @param nrow Number of rows for the grid (default: 1).
#' @param ncol Number of columns for the grid (default: 1).
#' @param shapefile An optional shapefile containing regions of interest (ROIs)
#'   for analysis.
#' @param build_shapefile Logical, indicating whether to interactively draw ROIs
#'   if shapefile is `NULL` (default: TRUE).
#' @param check_shapefile Logical, indicating whether to validate the shapefile
#'   with an interactive map view (default: FALSE). If `TRUE`, one can edit the
#'   drawn shapefile by deleting or changing the drawn grids.
#' @param buffer_edge Width of the buffer around the shapefile (default: 5).
#' @param buffer_col,buffer_row Buffering factor for the columns and rows,
#'   respectively, of each individual plot's side. A value between 0 and 0.5
#'   where 0 means no buffering and 0.5 means complete buffering (default: 0). A
#'   value of 0.25 will buffer the plot by 25% on each side.
#' @param segment_plot Logical, indicating whether to segment plots (default:
#'   FALSE). If `TRUE`, the `segment_index` will be computed and pixels with
#'   values below the `threshold` will be selected.
#' @param segment_individuals Logical, indicating whether to segment individuals
#'   within plots (default: FALSE). If `TRUE`, the `segment_index` will be
#'   computed and pixels with values below the `threshold` will be selected and
#'   a watershed-based segmentation will be performed.
#' @param watershed 	If `TRUE` (default) performs watershed-based object
#'   detection. This will detect objects even when they are touching one other.
#'   If FALSE, all pixels for each connected set of foreground pixels are set to
#'   a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Higher value smooths out small objects.
#' @param include_if Character vector specifying the type of intersection
#'   defaults to "centroid" (individuals in with the centroid is included within
#'   the drawn plot will be included in such plot). Other possible values
#'   includes `"covered"`, `"overlap"`, and `"intersect"`. See Details to a
#'   detailed explanation about these intersecting controls.
#' @param plot_index The index(es) to be computed for the drawn plots. Either a
#'   single vegetation index (eg., `"GLAI"`), a vector of indexes (eg.,
#'   `c("GLAI", "NGRDI", "HUE")`, or an own index based on the available bands
#'   (eg., `"(R-B)/(R+B)"`. See [pliman_indexes()] and [image_index()] for more
#'   details.
#' @param segment_index The index used for segmentation. The same rule as
#'   `plot_index`.
#' @param threshold By default (threshold = "Otsu"), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is informed, this value will be used as a threshold.
#' @param filter Performs median filtering in the binary image? See more at
#'   image_filter(). Defaults to FALSE. Use a positive integer to define the
#'   size of the median filtering. Larger values are effective at removing
#'   noise, but adversely affect edges.
#' @param summarize_fun The function to compute summaries for the pixel values.
#'   Defaults to "mean", ie., the mean value of the pixels (either at a plot- or
#'   individual-level) is returned.
#' @param attribute The attribute to be shown at the plot when `plot` is `TRUE`
#'   (default: "GLAI").
#' @param invert Logical, indicating whether to invert the mask (default: FALSE).
#' @param color_regions The color palette for regions (default:
#'   rev(grDevices::terrain.colors(50)).
#' @param plot Logical, indicating whether to generate plots (default: TRUE).
#' @param verbose Logical, indicating whether to display verbose output
#'   (default: TRUE).
#'
#' @return A list containing the following objects:
#' * `result_plot`: The results at a plot level
#' * `result_plot_summ`: The summary of results at a plot level. When `segment_individuals = TRUE`, the number of individuals, canopy coverage, and mean values of some shape statistics such as perimeter, length, width, and diameter are computed.
#' * `result_individ`: The results at an individual level.
#' * `map_plot`: An object of class `mapview` showing the plot-level results.
#' * `map_individual`: An object of class `mapview` showing the individual-level results.
#' * `shapefile`: The generated shapefile, with the drawn grids/blocks
#' @export
#'
#' @examples
#' library(pliman)
mosaic_analyze <- function(mosaic,
                           r = 1,
                           g = 2,
                           b = 3,
                           re = 4,
                           nir = 5,
                           grid = TRUE,
                           nrow = 1,
                           ncol = 1,
                           shapefile = NULL,
                           build_shapefile = TRUE,
                           check_shapefile = FALSE,
                           buffer_edge = 5,
                           buffer_col = 0,
                           buffer_row = 0,
                           segment_plot = FALSE,
                           segment_individuals = FALSE,
                           watershed = TRUE,
                           tolerance = 1,
                           extension = 1,
                           include_if = "centroid",
                           plot_index = "GLAI",
                           segment_index = "GLAI",
                           threshold = "Otsu",
                           filter = FALSE,
                           summarize_fun = "mean",
                           attribute = "mean.GLAI",
                           invert = FALSE,
                           color_regions = rev(grDevices::terrain.colors(50)),
                           alpha = 1,
                           max_pixels = 2000000,
                           downsample = NULL,
                           quantiles = c(0, 1),
                           plot = TRUE,
                           verbose = TRUE){
  includeopt <- c("intersect", "covered", "overlap", "centroid")
  includeopt <- includeopt[sapply(include_if, function(x){pmatch(x, includeopt)})]
  if(any(segment_individuals) & !segment_index %in% plot_index){
    plot_index <- unique(append(plot_index, segment_index))
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  nlyrs <- terra::nlyr(mosaic)
  if(verbose){
    cat("\014","\nBuilding the mosaic...\n")
  }
  if(is.null(shapefile)){

    created_shapes <- shapefile_build(mosaic,
                                      r = r,
                                      g = g,
                                      b = b,
                                      re = re,
                                      nir = nir,
                                      grid = grid,
                                      nrow = nrow,
                                      ncol = ncol,
                                      build_shapefile = build_shapefile,
                                      check_shapefile = check_shapefile,
                                      buffer_edge = buffer_edge,
                                      buffer_col = buffer_col,
                                      buffer_row = buffer_row,
                                      max_pixels = max_pixels,
                                      verbose = verbose,
                                      downsample = downsample,
                                      quantiles = quantiles)


    # if(build_shapefile){
    #   points <- mosaic_view(mosaic,
    #                         r = r,
    #                         g = g,
    #                         b = b,
    #                         re = re,
    #                         nir = nir,
    #                         max_pixels = max_pixels,
    #                         downsample = downsample,
    #                         quantiles = quantiles,
    #                         edit = TRUE)
    #   cpoints <- points$finished
    # } else{
    #   extm <- terra::ext(mosaic)
    #   xmin <- extm[1]
    #   xmax <- extm[2]
    #   ymin <- extm[3]
    #   ymax <- extm[4]
    #   coords <- matrix(c(xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin, xmin, ymax), ncol = 2, byrow = TRUE)
    #   # Create a Polygon object
    #   polygon <- sf::st_polygon(list(coords))
    #   # Create an sf object with a data frame that includes the 'geometry' column
    #   cpoints <- sf::st_sf(data.frame(id = 1),
    #                        geometry = sf::st_sfc(polygon),
    #                        crs = sf::st_crs(mosaic))
    # }


    # crop to the analyzed area
    poly_ext <-
      do.call(rbind, lapply(created_shapes, function(x){
        x
      })) |>
      sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
      terra::vect() |>
      terra::buffer(buffer_edge) |>
      terra::ext()

    mosaiccr <- terra::crop(mosaic, poly_ext)
    mosaiccr[mosaiccr == 65535] <- NA

    # # check the parameters
    # if(length(nrow) == 1 & nrow(cpoints) != 1){
    #   nrow <- rep(nrow, nrow(cpoints))
    # }
    # if(length(nrow) != nrow(cpoints)){
    #   warning(paste0("`nrow` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
    # }
    # if(length(ncol) == 1 & nrow(cpoints) != 1){
    #   ncol <- rep(ncol, nrow(cpoints))
    # }
    # if(length(ncol) != nrow(cpoints)){
    #   warning(paste0("`ncol` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
    # }
    # if(length(buffer_col) == 1 & nrow(cpoints) != 1){
    #   buffer_col <- rep(buffer_col, nrow(cpoints))
    # }
    # if(length(buffer_col) != nrow(cpoints)){
    #   warning(paste0("`buffer_col` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
    # }
    # if(length(buffer_row) == 1 & nrow(cpoints) != 1){
    #   buffer_row <- rep(buffer_row, nrow(cpoints))
    # }
    # if(length(buffer_row) != nrow(cpoints)){
    #   warning(paste0("`buffer_row` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
    # }
    # if(length(grid) == 1 & nrow(cpoints) != 1){
    #   grid <- rep(grid, nrow(cpoints))
    # }
    # if(length(grid) != nrow(cpoints)){
    #   warning(paste0("`grid` must have length 1 or ", nrow(cpoints), " (the number of drawn polygons)."))
    # }
    if(length(segment_plot) == 1 & length(created_shapes) != 1){
      segment_plot <- rep(segment_plot, length(created_shapes))
    }
    if(length(segment_plot) != length(created_shapes)){
      warning(paste0("`segment_plot` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(segment_individuals) == 1 & length(created_shapes) != 1){
      segment_individuals <- rep(segment_individuals, length(created_shapes))
    }
    if(length(segment_individuals) != length(created_shapes)){
      warning(paste0("`segment_individuals` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(threshold) == 1 & length(created_shapes) != 1){
      threshold <- rep(threshold, length(created_shapes))
    }
    if(length(threshold) != length(created_shapes)){
      warning(paste0("`threshold` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(segment_index) == 1 & length(created_shapes) != 1){
      segment_index <- rep(segment_index, length(created_shapes))
    }
    if(length(segment_index) != length(created_shapes)){
      warning(paste0("`segment_index` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(invert) == 1 & length(created_shapes) != 1){
      invert <- rep(invert, length(created_shapes))
    }
    if(length(invert) != length(created_shapes)){
      warning(paste0("`invert` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(includeopt) == 1 & length(created_shapes) != 1){
      includeopt <- rep(includeopt, length(created_shapes))
    }
    if(length(includeopt) != length(created_shapes)){
      warning(paste0("`includeopt` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(filter) == 1 & length(created_shapes) != 1){
      filter <- rep(filter, length(created_shapes))
    }
    if(length(filter) != length(created_shapes)){
      warning(paste0("`filter` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(length(grid) == 1 & length(created_shapes) != 1){
      grid <- rep(grid, length(created_shapes))
    }
    if(length(grid) != length(created_shapes)){
      warning(paste0("`grid` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    # # check the created shapes?
    # if(verbose){
    #   cat("\014","\nCreating the shapes...\n")
    # }
    # cpoints <- re$shapefile$finished
    # created_shapes <- list()
    # for(k in 1:nrow(cpoints)){
    #   if(inherits(cpoints[k, ]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(cpoints[k, ])) == 5 & grid[[k]]){
    #     plot_grid <-
    #       make_grid(cpoints[k, ],
    #                 nrow = nrow[k],
    #                 ncol = ncol[k],
    #                 buffer_col = buffer_col[k],
    #                 buffer_row = buffer_row[k])
    #   } else{
    #     plot_grid <-  cpoints[k, ] |> sf::st_transform(sf::st_crs(mosaic)) |> poorman::select(geometry)
    #   }
    #   if(check_shapefile){
    #     pg_edit <-
    #       plot_grid |>
    #       poorman::mutate(`_leaflet_id` = 1:nrow(plot_grid), feature_type = "polygon") |>
    #       poorman::relocate(geometry, .after = 3) |>
    #       sf::st_transform(crs = 4326)
    #     if(nlyrs > 2){
    #       map <-
    #         mapview::mapview() %>%
    #         add_rgb(x = as(mosaiccr, "Raster"),
    #                 r = r,
    #                 g = g,
    #                 b = b,
    #                 na.color = "#00000000",
    #                 quantiles = quantiles)
    #     } else{
    #       map <-
    #         mapview::mapview() %>%
    #         leafem::addGeoRaster(x = as(mosaic[[1]], "Raster"),
    #                              colorOptions = leafem:::colorOptions(palette = color_regions,
    #                                                                   na.color = "transparent"))
    #     }
    #     if(build_shapefile){
    #       mapview::mapview() |> mapedit::editMap()
    #     }
    #     edited <-
    #       mapedit::editFeatures(pg_edit, map) |>
    #       poorman::select(geometry) |>
    #       sf::st_transform(sf::st_crs(mosaic))
    #     created_shapes[k] <- edited
    #   } else{
    #     created_shapes[k] <- plot_grid
    #   }
    # }

  } else{
    extm <- terra::ext(shapefile)
    xmin <- extm[1]
    xmax <- extm[2]
    ymin <- extm[3]
    ymax <- extm[4]
    coords <- matrix(c(xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin, xmin, ymax), ncol = 2, byrow = TRUE)
    # Create a Polygon object
    # polygon <- sf::st_polygon(list(coords))
    # Create an sf object with a data frame that includes the 'geometry' column
    # cpoints <- sf::st_sf(data.frame(id = 1),
    #                      geometry = sf::st_sfc(polygon),
    #                      crs = sf::st_crs(mosaic))
    geoms <- sf::st_as_sf(shapefile$geometry)
    sf::st_geometry(geoms) <- "geometry"
    created_shapes <- list(geoms)

    # crop to the analyzed area
    poly_ext <-
      created_shapes[[1]] |>
      sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
      terra::vect() |>
      terra::buffer(buffer_edge) |>
      terra::ext()
    mosaiccr <- terra::crop(mosaic, poly_ext)
    mosaiccr[mosaiccr == 65535] <- NA
  }

  # return(created_shapes)
  # compute the indexes
  if(verbose){
    cat("\014","\nComputing the indexes...\n")
  }
  if(nlyrs > 1){
    mind <- terra::rast(
      Map(c,
          lapply(seq_along(plot_index), function(i){
            mosaic_index(mosaiccr,
                         index = plot_index[[i]],
                         r = r,
                         g = g,
                         b = b,
                         re = re,
                         nir = nir)
          })
      )
    )
  } else{
    plot_index <- names(mosaiccr)
    mind <- mosaiccr
  }

  results <- list()
  result_indiv <- list()
  extends <- terra::ext(mosaiccr)
  # return(created_shapes)
  # return(created_shapes)
  for(j in seq_along(created_shapes)){
    if(segment_plot[j] & segment_individuals[j]){
      stop("Only `segment_plot` OR `segment_individuals` can be used", call. = FALSE)
    }
    if(verbose){
      cat("\014","\nExtracting data from block", j, "\n")
    }
    if(inherits(created_shapes[[j]]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(created_shapes[[j]][[1]][[1]])) == 5 & grid[[j]]){
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      ext_anal <-
        plot_grid |>
        sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
        terra::vect() |>
        terra::buffer(buffer_edge) |>
        terra::ext()
      mind_temp <- terra::crop(mind, terra::ext(ext_anal))
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(!segment_index[j] %in% names(mind_temp)){
          stop("`segment_index` must be one of used in `plot_index`.")
        }
        thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
        print(thresh)
        if(invert[j]){
          mask <- mind_temp[[segment_index[j]]] > thresh
        } else{
          mask <- mind_temp[[segment_index[j]]] < thresh
        }
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
      }
      # check if segmentation is performed (analyze individuals)
      if(segment_individuals[j]){
        thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
        if(invert[j]){
          mask <- mind_temp[[segment_index[j]]] > thresh
        } else{
          mask <- mind_temp[[segment_index[j]]] < thresh
        }
        dmask <- EBImage::Image(matrix(mask, ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(watershed){
          dmask <- EBImage::watershed(EBImage::distmap(dmask), tolerance = tolerance, ext = extension)
        } else{
          dmask <- EBImage::bwlabel(dmask)
        }
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        conts <- EBImage::ocontour(matrix(dmask, ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        conts <- conts[sapply(conts, nrow) > 2]
        sf_df <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extends[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extends[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        )
        centroids <- suppressWarnings(sf::st_centroid(sf_df))
        intersects <-
          switch (includeopt[j],
                  "intersect" = sf::st_intersects(sf_df, plot_grid),
                  "centroid" =  sf::st_within(centroids, plot_grid),
                  "covered" = sf::st_covered_by(sf_df, plot_grid),
                  "overlap" = sf::st_overlaps(sf_df, plot_grid),

          )
        plot_id <- data.frame(plot_id = paste0(intersects))
        valid_rows <- plot_id$plot_id != "integer(0)"
        sf_df <- sf_df[valid_rows, ]
        plot_id <- plot_id[valid_rows, ]
        addmeasures <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  }))
        gridindiv <- cbind(sf_df, plot_id, addmeasures)[c(2, 1, 3:10)]

        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = gridindiv,
                                       fun = summarize_fun,
                                       progress = FALSE,
                                       force_df = TRUE)

        if(inherits(valindiv, "list")){
          valindiv <-
            do.call(rbind, lapply(1:length(valindiv), function(i){
              tmp <- transform(valindiv[[i]],
                               individual = paste0(i),
                               block = paste0("B", j))
              tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

            }
            ))
          if("coverage_fraction" %in% colnames(valindiv)){
            valindiv$coverage_fraction <- NULL
          }
        }
        if(ncol(valindiv) == 1 & length(plot_index) == 1){
          colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
        }
        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", j), gridindiv, valindiv)
          # valindiv <- valindiv[valindiv$plot_id != "integer(0)", ]
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        } else{
          result_indiv[[j]] <- valindiv
        }
      } else{
        dmask <- NULL
        result_indiv[[j]] <- NULL
      }

      # extract the values for the individual plots
      vals <-
        exactextractr::exact_extract(x = mind_temp,
                                     y = plot_grid,
                                     fun = summarize_fun,
                                     progress = FALSE,
                                     force_df = TRUE)

    } else{
      ####### ANY TYPE OF POLYGON ########
      # check if segmentation is performed
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      ext_anal <-
        plot_grid |>
        terra::vect() |>
        terra::buffer(20) |>
        terra::ext()
      mind_temp <- terra::crop(mind, terra::ext(ext_anal))
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(!segment_index[j] %in% names(mind_temp)){
          stop("`segment_index` must be one of used in `plot_index`.")
        }
        thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
        if(invert[j]){
          mask <- mind_temp[[segment_index[j]]] > thresh
        } else{
          mask <- mind_temp[[segment_index[j]]] < thresh
        }
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
      }

      if(segment_individuals[j]){
        thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
        extends <- terra::ext(mind_temp)
        if(invert[j]){
          mask <- mind_temp[[segment_index[j]]] > thresh
        } else{
          mask <- mind_temp[[segment_index[j]]] < thresh
        }
        dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(watershed){
          dmask <- EBImage::watershed(EBImage::distmap(dmask), tolerance = tolerance, ext = extension)
        } else{
          dmask <- EBImage::bwlabel(dmask)
        }
        conts <- EBImage::ocontour(dmask)
        conts <- conts[sapply(conts, nrow) > 2]
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        sf_df <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extends[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extends[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        )
        centroids <- suppressWarnings(sf::st_centroid(sf_df))
        # intersect_individ <- sf::st_overlaps(sf_df, plot_grid, sparse = FALSE)[,1]
        intersect_individ <-
          switch (includeopt[j],
                  "intersect" = sf::st_intersects(sf_df, plot_grid, sparse = FALSE)[,1],
                  "centroid" = sf::st_within(centroids, plot_grid, sparse = FALSE)[,1],
                  "covered" = sf::st_covered_by(sf_df, plot_grid, sparse = FALSE)[,1],
                  "overlap" = sf::st_overlaps(sf_df, plot_grid, sparse = FALSE)[,1],

          )
        sf_df <- sf_df[intersect_individ, ]
        addmeasures <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  }))
        gridindiv <- cbind(sf_df, addmeasures)

        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = gridindiv,
                                       fun = summarize_fun,
                                       progress = FALSE,
                                       force_df = TRUE)

        if(inherits(valindiv, "list")){
          valindiv <-
            do.call(rbind, lapply(1:length(valindiv), function(i){
              tmp <- transform(valindiv[[i]],
                               individual = paste0(i),
                               block = paste0("B", j))
              tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]
            }
            ))
          if("coverage_fraction" %in% colnames(valindiv)){
            valindiv$coverage_fraction <- NULL
          }
        }
        if(ncol(valindiv) == 1 & length(plot_index) == 1){
          colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
        }

        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", j), plot_id = 1, gridindiv, valindiv)
          result_indiv[[j]] <- valindiv
        } else{
          result_indiv[[j]] <- valindiv
        }
      } else{
        result_indiv[[j]] <- NULL
      }

      vals <-
        exactextractr::exact_extract(x = mind_temp,
                                     y = plot_grid,
                                     fun = summarize_fun,
                                     progress = FALSE,
                                     force_df = TRUE)
    }

    # bind the results
    if(inherits(vals, "list")){
      vals <-
        do.call(rbind, lapply(1:length(vals), function(i){
          tmp <- transform(vals[[i]],
                           plot_id = paste0(i),
                           block = paste0("B", j))
          tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]
        }
        ))
      if("coverage_fraction" %in% colnames(vals)){
        vals$coverage_fraction <- NULL
      }
      if(ncol(vals) == 1 & length(plot_index) == 1){
        colnames(vals) <- paste0(colnames(vals), ".", plot_index)
      }
    } else{
      if(ncol(vals) == 1 & length(plot_index) == 1){
        colnames(vals) <- paste0(colnames(vals), ".", plot_index)
      }
      vals <- transform(vals,
                        plot_id = paste0(1:nrow(vals)),
                        block = paste0("B", j))
      vals <- vals[, c(ncol(vals), ncol(vals) - 1, 1:(ncol(vals) - 2))]
    }
    if(!is.null(summarize_fun)){
      results[[j]] <- cbind(plot_grid, vals)
    } else{
      results[[j]] <- vals
    }
    # end
  }



  # bind the results  ## at a level plot
  results <- do.call(rbind, lapply(results, function(x){x})) |> sf::st_sf()

  if(any(segment_individuals)){

    result_indiv <- do.call(rbind, lapply(result_indiv, function(x){x}))
    blockid <- unique(result_indiv$block)
    summres <-
      lapply(1:length(blockid), function(i){
        # tmp <-
        result_indiv |>
          poorman::filter(block == blockid[i]) |>
          as.data.frame() |>
          poorman::group_by(plot_id) |>
          poorman::summarise(poorman::across(poorman::where(is.numeric), \(x){mean(x, na.rm = TRUE)}),
                             n = length(area),
                             area_sum = sum(area, na.rm = TRUE)) |>
          poorman::mutate(block = blockid[i], .before = 1) |>
          poorman::ungroup() |>
          poorman::relocate(n, .after = plot_id)
      })
    names(summres) <- blockid
    # compute plot area
    plot_area <-
      results |>
      poorman::mutate(plot_area = sf::st_area(geometry)) |>
      as.data.frame() |>
      poorman::select(block, plot_id, plot_area)
    # compute coverage area
    result_plot_summ <-
      do.call(rbind, lapply(summres, function(x){x})) |>
      poorman::left_join(plot_area, by = c("block", "plot_id")) |>
      poorman::mutate(coverage = as.numeric(area_sum / plot_area), .after = area) |>
      poorman::left_join(results |>   poorman::select(block, plot_id, geometry), by = c("block", "plot_id")) |>
      sf::st_as_sf()

  } else{
    result_plot_summ <- NULL
    result_indiv <- NULL
  }

  ext_plot <-
    poly_ext |>
    terra::vect() |>
    terra::buffer(buffer_edge) |>
    terra::ext()

  mosaic_extent <- terra::ext(mosaic)
  # Ensure the extent of ext_plot is constrained within the mosaic_extent
  ext_plot[1] <- max(ext_plot[1], mosaic_extent[1])
  ext_plot[2] <- min(ext_plot[2], mosaic_extent[2])
  ext_plot[3] <- max(ext_plot[3], mosaic_extent[3])
  ext_plot[4] <- min(ext_plot[4], mosaic_extent[4])
  mosaicplot <- terra::crop(mosaic, ext_plot)
  if(!is.null(summarize_fun) & isTRUE(plot)){
    if(verbose){
      cat("\014","\nPreparing to plot...\n")
    }

    possible_downsamples <- 0:15
    possible_npix <- sapply(possible_downsamples, function(x){
      compute_downsample(nrow(mosaicplot), ncol(mosaicplot), x)
    })
    downsample <- which.min(abs(possible_npix - max_pixels))
    downsample <- ifelse(downsample == 1, 0, downsample)
    if(downsample > 0){
      mosaicplot <- terra::aggregate(mosaicplot, fact = downsample)
    }
    map <-
      suppressWarnings(
        mapview::mapview(results,
                         zcol = attribute,
                         layer.name = attribute,
                         col.regions = custom_palette(c("red", "yellow", "darkgreen"), n = 3),
                         alpha.regions = 0.8,
                         na.color = "#00000000",
                         maxBytes = 64 * 1024 * 1024,
                         verbose = FALSE)
      )

    map <-
      add_rgb(map,
              as(mosaicplot, "Raster"),
              r = r,
              g = g,
              b = b,
              na.color = "#00000000",
              quantiles = quantiles)

    if(any(segment_individuals)){
      mapindivid <-
        suppressWarnings(
          mapview::mapview(result_indiv,
                           zcol = attribute,
                           layer.name = attribute,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE) +
            mapview::mapview(results,
                             legend = FALSE,
                             alpha.regions = 0.4,
                             zcol = "block",
                             map.types = "OpenStreetMap") )|>
        add_rgb(as(mosaicplot, "Raster"),
                r = r,
                g = g,
                b = b,
                na.color = "#00000000",
                quantiles = quantiles)
    } else{
      mapindivid <- NULL
    }
  } else{
    map <- NULL
    mapindivid <- NULL
  }

  if(verbose){
    cat("\014","Done!\n")
  }
  return(list(result_plot = results,
              result_plot_summ = result_plot_summ,
              result_indiv = result_indiv,
              map_plot = map,
              map_indiv = mapindivid,
              shapefile = created_shapes))
}



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
#' @param edit If `TRUE` enable editing options using [mapedit::editMap()].
#' @param title A title for the generated map or plot (default: "").
#' @param max_pixels Maximum number of pixels to render in the map or plot
#'   (default: 500000).
#' @param downsample Downsampling factor to reduce the number of pixels
#'   (default: NULL). In this case, if the number of pixels in the image (width
#'   x height) is greater than `max_pixels` a downsampling factor will be
#'   automatically chosen so that the number of plotted pixels approximates the
#'   `max_pixels`.
#' @param alpha opacity of the fill color of the raster layer(s).
#' @param quantiles the upper and lower quantiles used for color stretching.
#' @param axes logical. Draw axes? Defaults to `FALSE`.
#' @param ... Additional arguments passed on to [terra::plot()] when `viewer =
#'   "base"`.
#' @return An sf object, the same object returned by [mapedit::editMap()].
#'
#' @importFrom terra rast crs nlyr terraOptions
#' @importFrom methods as
#' @importFrom sf st_crs st_transform st_make_grid
#' @importFrom stars st_downsample st_as_stars
#' @importFrom poorman summarise across mutate arrange left_join bind_cols
#'   bind_rows contains ends_with everything between pivot_longer pivot_wider
#'   where select filter relocate rename
#' @importFrom htmlwidgets JS
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
                        edit = FALSE,
                        title = "",
                        viewer = c("mapview", "base"),
                        show = c("rgb", "index"),
                        index = "B",
                        max_pixels = 1000000,
                        downsample = NULL,
                        alpha = 1,
                        quantiles = c(0, 1),
                        color_regions = custom_palette(),
                        axes = FALSE,
                        ...){
  terra::terraOptions(progress = 0)
  on.exit(terra::terraOptions(progress = 1))
  # mosaic[mosaic == 65535] <- NA
  # mosaic <- terra::subst(mosaic, 65535, NA)
  # check_mapview()
  mapview::mapviewOptions(layers.control.pos = "topright", raster.size = 64 * 1024 * 1024)
  on.exit(mapview::mapviewOptions(default = TRUE))
  viewopt <- c("rgb", "index")
  viewopt <- viewopt[pmatch(show[[1]], viewopt)]
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]

  if(inherits(mosaic, "Image")){
    mosaic <- terra::rast(EBImage::transpose(mosaic)@.Data)
  }
  if(viewopt == "rgb" & vieweropt == "base" & terra::nlyr(mosaic) > 1){
    message("`viewer = 'base' can only be used with `show = 'index'`. Defaulting to viewer = 'mapview'")
    vieweropt <- "mapview"
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("+proj=utm +zone=05 +datum=WGS84 +units=m")
  }
  dimsto <- dim(mosaic)
  nr <- dimsto[1]
  nc <- dimsto[2]

  if(max_pixels > 2000000){
    message("The number of pixels is too high, which might slow the rendering process.")
  }
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
    mosaic <- terra::aggregate(mosaic, fact = downsample)
  }
  if(viewopt == "index" & terra::nlyr(mosaic) > 2){
    mosaic <- mosaic_index(mosaic, index = index)
  }
  if(viewopt == "rgb"){
    if(terra::nlyr(mosaic) > 2){
      map <-
        mapview::viewRGB(as(mosaic, "Raster"),
                         na.color = "#00000000",
                         layer.name = "T",
                         r = r,
                         g = g,
                         b = b,
                         maxpixels = 60000000,
                         quantiles = quantiles)
      if(edit){
        map <-
          mapedit::editMap(map,
                           editor = "leafpm",
                           title = title)
      }
      map

    } else{
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        index <- gsub("[/\\\\]", "_", index, perl = TRUE)
        map <-
          mapview::mapview(mosaic,
                           map.types = mapview::mapviewGetOption("basemaps"),
                           layer.name = index,
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)

        if(edit){
          map <-
            map |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        }
        map
      }
    }
  } else{
    if(terra::nlyr(mosaic) > 2){
      index <- gsub("[/\\\\]", "_", index, perl = TRUE)
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(mosaic,
                           layer.name = index,
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
        if(edit){
          map <-
            map |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        }
        map
      }
    } else{
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        map <-
          mapview::mapview(mosaic,
                           layer.name = names(mosaic),
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = 1,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
        if(edit){
          map <-
            map |>
            mapedit::editMap(editor = "leafpm",
                             title = title)
        }
        map
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
#' @param info Print the mosaic informations (eg., CRS, extend). Defaults to `TRUE`
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
mosaic_input <- function(mosaic, info = TRUE, ...){
  mosaic <- suppressWarnings(terra::rast(mosaic, ...))
  if(terra::crs(mosaic) == ""){
    message("Missing Coordinate Reference System. Setting to EPSG:3857")
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  if(info){
    print(mosaic)
  }
  return(mosaic)
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
    do.call(rbind, lapply(shapefile, create_shp, info, as_sf, ...))
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
shapefile_view <- function(shapefile, ...){
  mapview::mapview(shapefile, ...)
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
#' @importFrom terra crs
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
                          edit = TRUE,
                          title = "Use the 'Draw rectangle' tool to select the cropping area.",
                          ...)
  if(!is.na(sf::st_crs(mosaic))){
    grids <-
      sf::st_make_grid(controls$finished, n = c(1, 1)) |>
      sf::st_transform(sf::st_crs(mosaic))
  } else{
    terra::crs(mosaic) <- terra::crs("+proj=utm +zone=32 +datum=WGS84 +units=m")
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
#' mosaic_view(final, viewer = "base")
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
      eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index])))
  } else{
    mosaic_gray <-
      eval(parse(text = as.character(index)))
  }
  names(mosaic_gray) <- index
  if(!is.na(terra::crs(mosaic))){
    suppressWarnings(terra::crs(mosaic_gray) <- terra::crs(mosaic))
  } else{
    suppressWarnings(terra::crs(mosaic_gray) <- "+proj=utm +zone=32 +datum=WGS84 +units=m")
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
#' @export
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
                           crop_mosaic = TRUE,
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
      image_view(ebimg[1:5, 1:5,], edit = TRUE)
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
    aligned <- image_align(ebimg, viewer = vieweropt)
    if(vieweropt != "base"){
      image_view(aligned[1:5, 1:5,], edit = TRUE)
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




#' Drawing Lines or Polygons with Raster Information
#'
#' @details
#' The `mosaic_draw` function enables you to create mosaic drawings from
#' remote sensing data and compute vegetation indices.
#'
#'  * If a line is drawn using the "Draw Polyline" tool, the profile of `index` is
#'  displayed on the y-axis along the line's distance, represented in meter
#'  units. It is important to ensure that the Coordinate Reference System (CRS)
#'  of `mosaic` has latitude/longitude units for accurate distance
#'  representation.
#'
#'  * If a rectangle or polygon is drawn using the "Draw Rectangle" or "Draw Polygon"
#'  tools, the `index` values are calculated for each object. By default, the
#'  raw data is returned. You can set the `summarize_fun` to compute a summary
#'  statistic for each object.
#'
#' @inheritParams mosaic_view
#' @inheritParams analyze_objects
#' @param color_regions The color palette for displaying index values. Defaults
#'   to `rev(grDevices::terrain.colors(50))`.
#' @param threshold By default (threshold = "Otsu"), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is informed, this value will be used as a threshold.
#' @param invert Inverts the mask if desired. Defaults to `FALSE`.
#' @param segment Should the raster object be segmented? If set to `TRUE`,
#'   pixels within each polygon/rectangle will be segmented based on the
#'   `threshold` argument.
#' @param summarize_fun An optional function or character vector. When
#'   `summarize_fun = "mean"`, the mean values of `index` are calculated within
#'   each object. For more details on available functions, refer to
#'  [exactextractr::exact_extract()].
#' @param buffer Adds a buffer around the geometries of the SpatVector created.
#'   Note that the distance unit of `buffer` will vary according to the CRS of
#'   `mosaic`.
#' @param plot Plots the draw line/rectangle? Defaults to `TRUE`.
#' @param plot_layout The de plot layout. Defaults to `plot_layout = c(1, 2, 3,
#'   3)`. Ie., the first row has two plots, and the second row has one plot.
#' @importFrom terra crop vect extract
#' @importFrom exactextractr exact_extract
#' @importFrom graphics layout
#' @importFrom stats smooth
#' @return An invisible list containing the mosaic, draw_data, distance,
#'   distance_profile, geometry, and map.
#' @export
#' @examples
#'
#' if(interactive()){
#' library(pliman)
#' # Load a raster showing the elevation of Luxembourg
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'
#' # draw a polyline to see the elevation profile along the line
#' mosaic_draw(mosaic, buffer = 1500)
#' }
mosaic_draw <- function(mosaic,
                        r = 3,
                        g = 2,
                        b = 1,
                        re = 4,
                        nir = 5,
                        index = "NGRDI",
                        show = "rgb",
                        segment = FALSE,
                        viewer = c("mapview", "base"),
                        threshold = "Otsu",
                        invert = FALSE,
                        summarize_fun = NULL,
                        buffer = 2,
                        color_regions = rev(grDevices::terrain.colors(50)),
                        alpha = 1,
                        max_pixels = 1000000,
                        downsample = NULL,
                        quantiles = c(0, 1),
                        plot = TRUE,
                        plot_layout = c(1, 2, 3, 3)){
  if(is.null(terra::crs(mosaic))){
    terra::crs(mosaic) <- "+proj=utm +zone=32 +datum=WGS84 +units=m"
  }
  vieweropt <- c("base", "mapview")
  vieweropt <- vieweropt[pmatch(viewer[[1]], vieweropt)]

  points <- mosaic_view(mosaic,
                        r = r,
                        g = g,
                        b = b,
                        re = re,
                        nir = nir,
                        max_pixels = max_pixels,
                        downsample = downsample,
                        alpha = alpha,
                        quantiles = quantiles,
                        index = index[[1]],
                        show = show,
                        edit = TRUE)

  nlyrs <- terra::nlyr(mosaic)
  polygons <- points$finished$geometry
  polygons_spv <- sf::st_transform(polygons, crs = sf::st_crs(mosaic))
  polygons_ext <- terra::vect(polygons_spv)
  ext <- terra::buffer(polygons_ext, buffer) |> terra::ext()
  mosaiccr <- terra::crop(mosaic, ext)
  mosaiccr[mosaiccr == 65535] <- NA

  # Compute the image indexes
  if(nlyrs > 1){
    mind <- terra::rast(
      Map(c,
          lapply(seq_along(index), function(i){
            mosaic_index(mosaiccr,
                         index = index[[i]],
                         r = r,
                         g = g,
                         b = b,
                         re = re,
                         nir = nir)
          })
      )
    )
  } else{
    index <- names(mosaiccr)
    mind <- mosaiccr
  }

  if(inherits(polygons, "sfc_LINESTRING")){
    vals <-
      terra::extract(x = mind,
                     y = polygons_ext,
                     fun = summarize_fun)
    coords <- as.matrix(polygons_spv[[1]])
    n <- nrow(coords)
    distances <- NULL
    for (j in 1:(n - 1)) {
      x1 <- coords[j, 1]
      y1 <- coords[j, 2]
      x2 <- coords[j + 1, 1]
      y2 <- coords[j + 1, 2]
      distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
      distances[j] <- distance
    }
    # distances
    dists <- cumsum(distances)
    dist <- max(dists)

    if(plot){
      if(nlyrs > 2){
        layout(
          matrix(plot_layout, nrow = 2, byrow = TRUE),
          heights = c(3, 3)
        )
        on.exit(layout(1))
        scl <- max(terra::minmax(mosaiccr))

        terra::plotRGB(mosaiccr,
                       r = r,
                       g = g,
                       b = b,
                       scale = ifelse(scl < 255, 255, scl * 0.7),
                       colNA = "#00000000",
                       mar = c(2, 2, 2, 2),
                       axes = FALSE)
        lines(coords,
              col = "red",
              lwd = 3)

        terra::plot(mind[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)
        lines(coords,
              col = "red",
              lwd = 3)

        plot(x = seq(0, dist, length.out = nrow(vals)),
             y = smooth(vals[, index[[1]]]),
             xlab = "Distance",
             ylab = index[[1]],
             # mar = c(0, 0, 0, 0),
             type = "l",
             xlim = c(0, dist),
             col = "red")
      } else{
        layout(
          matrix(plot_layout, nrow = 2, byrow = TRUE),
          heights = c(3, 3)
        )
        on.exit(layout(1))
        ext2 <- terra::buffer(polygons_ext, buffer * 3) |> terra::ext()
        mosaiccr2 <- terra::crop(mosaic[[1]], ext2)
        terra::plot(mosaiccr2[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)


        terra::plot(mind[[1]],
                    axes = FALSE,
                    maxcell=5000000,
                    mar = c(2, 2, 2, 2),
                    smooth=TRUE)
        lines(coords,
              col = "red",
              lwd = 3)

        plot(x = seq(0, dist, length.out = nrow(vals)),
             y = smooth(vals[, index[[1]]]),
             xlab = "Distance",
             ylab = index[[1]],
             type = "l",
             xlim = c(0, dist),
             col = "red")

      }
    }
    map <- NULL
  } else{
    if(segment){
      mask <- mind[[1]] < otsu(na.omit(terra::values(mind)[, index[1]]))
      mind <- terra::mask(mind, mask, inverse = invert, maskvalues = TRUE)
    }
    mind <- terra::mask(mind, polygons_ext)
    vals <-
      suppressWarnings(
        exactextractr::exact_extract(x = mind,
                                     y = polygons,
                                     fun = "mean",
                                     progress = FALSE,
                                     force_df = TRUE)
      )
    if(inherits(vals, "list")){
      vals <-
        do.call(rbind, lapply(1:length(vals), function(i){
          tmp <- transform(vals[[i]], id = i)
          tmp[, c(ncol(tmp), 1:(ncol(tmp) - 2))]
        }
        ))
    } else{
      vals <- transform(vals, id = paste0(1:nrow(vals)))
      vals <- vals[, c(ncol(vals), 1:(ncol(vals) - 2))]
    }
    colnames(vals) <- c("id", index)
    vals <- na.omit(vals)

    if(nlyrs > 2){
      if(vieweropt == "base"){
        terra::plot(mind, axes = FALSE)
        map <- NULL
      } else{

        map <-
          mapview::viewRGB(as(mosaiccr, "Raster"),
                           na.color = "#00000000",
                           layer.name = "",
                           r = r,
                           g = g,
                           b = b,
                           maxpixels = 10000000,
                           quantiles = quantiles)

        for (i in 1:length(index)) {
          map <-
            mapview::mapview(mind[[i]],
                             hide = TRUE,
                             map = map,
                             layer.name = index[[i]],
                             map.types = mapview::mapviewGetOption("basemaps"),
                             maxpixels =  max_pixels,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE)
        }
        map

      }
    } else{
      if(vieweropt == "base"){
        terra::plot(mind, axes = FALSE)
        map <- NULL
      } else{
        map <-
          mapview::mapview(mind,
                           layer.name = names(mind),
                           map.types = mapview::mapviewGetOption("basemaps"),
                           maxpixels =  max_pixels,
                           col.regions = color_regions,
                           alpha.regions = alpha,
                           na.color = "#00000000",
                           maxBytes = 64 * 1024 * 1024,
                           verbose = FALSE)
      }
    }

    dists <- NULL
    dist <- NULL
  }

  invisible(list(
    mosaic = mosaic,
    draw_data = vals,
    distance = dist,
    distance_profile = dists,
    geometry = points,
    map = map

  ))
}


