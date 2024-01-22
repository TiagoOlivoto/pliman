sf_to_polygon <- function(shps) {
  if(inherits(shps, "list")){
    shps <- do.call(rbind, shps)
  }
  classes <- sapply(lapply(sf::st_geometry(shps$geometry), class), function(x){x[2]})
  shps[classes %in% c("POINT", "LINESTRING"), ] <-
    shps[classes %in% c("POINT", "LINESTRING"), ] |>
    sf::st_buffer(0.0000001) |>
    sf::st_cast("POLYGON") |>
    sf::st_simplify(preserveTopology = TRUE)
  return(shps)
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
                      mosaic,
                      buffer_col = 0,
                      buffer_row = 0){

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
  cvm <- lm(formula = txy ~ sxy[1:4, ])
  parms <- cvm$coefficients[2:3, ]
  intercept <- cvm$coefficients[1, ]
  geometry <- grids * parms + intercept

  if(buffer_row != 0 | buffer_col != 0){
    for(k in 1:length(geometry)){
      coords <- geometry[[k]][[1]]
      geometry[[k]][[1]] <- create_buffer(coords,
                                      buffer_col = buffer_row,
                                      buffer_row = buffer_col)
    }
  }
  gshp <-
    geometry |>
    sf::st_sf(crs = sf::st_crs(mosaic))
  return(gshp)
}

find_aggrfact <- function(mosaic, max_pixels = 1000000){
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
  nr <- nrow(mosaic)
  nc <- ncol(mosaic)
  npixel <- nr * nc
  possible_downsamples <- 0:20
  possible_npix <- sapply(possible_downsamples, function(x){
    compute_downsample(nr, nc, x)
  })
  downsample <- which.min(abs(possible_npix - max_pixels))
  downsample <- ifelse(downsample == 1, 0, downsample)
  return(downsample)
}
compute_measures_mosaic <- function(contour){
  lw <- help_lw(contour)
  cdist <- help_centdist(contour)
  data.frame(area = help_area(contour),
             perimeter = sum(help_distpts(contour)),
             length = lw[[1]],
             width = lw[[2]],
             diam_min = min(cdist) * 2,
             diam_mean = mean(cdist) * 2,
             diam_max = max(cdist) * 2)

}
map_individuals <- function(object,
                            by_column = "plot_id",
                            direction = c("horizontal", "vertical")) {
  optdirec <- c("horizontal", "vertical")
  optdirec <- pmatch(direction[[1]], optdirec)
  object <- as.data.frame(object)
  unique_values <- unique(object[, by_column] )
  distances <- vector("list", length(unique_values))
  for (i in 1:length(unique_values)) {
    subset_coords <- object[object[, by_column] == unique_values[i], 4:5]
    n <- nrow(subset_coords)
    nearest <- order(subset_coords[, optdirec])
    subset_distances <- numeric(n - 1)
    for (j in 1:(n - 1)) {
      x1 <- subset_coords[nearest[j], 1]
      y1 <- subset_coords[nearest[j], 2]
      x2 <- subset_coords[nearest[j+1], 1]
      y2 <- subset_coords[nearest[j+1], 2]
      distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
      subset_distances[j] <- distance
    }
    distances[[i]] <- subset_distances
  }
  if(optdirec == 1){
    names(distances) <- paste0("row", 1:length(distances))
  } else{
    names(distances) <- paste0("column", 1:length(distances))
  }
  cvs <- sapply(distances, function(x){
    (sd(x) / mean(x)) * 100
  })
  means <- sapply(distances, mean)
  invisible(list(distances = distances, cvs = cvs, means = means))
}
linear_iterpolation <- function(mosaic, points, method = "loess"){
  if(inherits(points, "list")){
    points <- do.call(rbind, points)
  }
  xy <- sf::st_coordinates(points)[, 1:2]
  vals <- terra::values(mosaic)[terra::cellFromXY(mosaic, xy), ]
  vals <- data.frame(cbind(xy, vals))
  names(vals) <- c("x", "y", "z")
  newdata <- as.data.frame(terra::xyFromCell(mosaic, 1:terra::ncell(mosaic)))
  new_ras <-
    terra::rast(
      lapply(3:ncol(vals), function(i){
        if(method == "loess"){
          mod <- loess(vals[, i] ~ x + y, data = vals)
        } else{
          mod <- lm(vals[, i] ~ x + y, data = vals)
        }
        terra::rast(matrix(predict(mod, newdata = newdata),
                           nrow = nrow(mosaic),
                           ncol = ncol(mosaic),
                           byrow = TRUE))
      })
    )
  terra::crs(new_ras) <- terra::crs(mosaic)
  terra::ext(new_ras) <- terra::ext(mosaic)
  terra::resample(new_ras, mosaic)
}
idw_interpolation <- function(mosaic, points){
  downsample <- find_aggrfact(mosaic, max_pixels = 200000)
  if(downsample > 0){
    magg <- terra::aggregate(mosaic, downsample)
  } else{
    magg <- mosaic
  }
  if(inherits(points, "list")){
    points <- do.call(rbind, points)
  }
  xy <- sf::st_coordinates(points)[, 1:2]
  vals <- terra::values(magg)[terra::cellFromXY(magg, xy), ]
  vals <- data.frame(cbind(xy, vals))

  xy_grid <- terra::xyFromCell(magg, 1:terra::ncell(magg))
  newx <- seq(min(xy_grid[,1]), max(xy_grid[,1]), length.out = 1000)
  newy <- seq(min(xy_grid[,2]), max(xy_grid[,2]), length.out = 1000)

  new_ras <-
    terra::rast(
      lapply(3:ncol(vals), function(i){
        interp <- idw_interpolation_cpp(vals[, 1], vals[, 2], vals[, i], xy_grid[, 1], xy_grid[, 2])
        ra3 <-
          terra::rast(matrix(interp,
                             nrow = nrow(magg),
                             ncol = ncol(magg),
                             byrow = TRUE))
      })
    )

  terra::crs(new_ras) <- terra::crs(mosaic)
  terra::ext(new_ras) <- terra::ext(mosaic)
  terra::resample(new_ras, mosaic)
}


#' Mosaic interpolation
#'
#' Performs the interpolation of points from a raster object.
#'
#' @param mosaic An `SpatRaster` object
#' @param points An `sf` object with the points for x and y coordinates, usually
#'   obtained with [shapefile_build()]. Alternatively, an external shapefile
#'   imported with [shapefile_input()] containing the x and y coordinates can be
#'   used. The function will handle most used shapefile formats (eg.,
#'   .shp, .rds) and convert the imported shapefile to an sf object.
#' @param method One of "bilinear" (default), "loess" (local regression) or
#'   "idw" (Inverse Distance Weighting).
#' @importFrom stats loess
#'
#' @return An `SpatRaster` object with the same extend and crs from `mosaic`
#' @export
#'
mosaic_interpolate <- function(mosaic, points, method = c("bilinear", "loess", "idw")){
  if(terra::crs(points) != terra::crs(mosaic)){
    terra::crs(points) <- terra::crs(mosaic)
  }
  if(!method[[1]] %in% c("bilinear", "idw", "loess")){
    stop("'method' must be one of 'bilinear', 'loess', or 'idw'")
  }
  if(method[[1]]  %in%  c("bilinear", "loess")){
    linear_iterpolation(mosaic, points, method = method[[1]])
  } else{
    idw_interpolation(mosaic, points)
  }
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
#' @inheritParams utils_shapefile
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
#'
shapefile_build <- function(mosaic,
                            basemap = NULL,
                            controlpoints = NULL,
                            r = 3,
                            g = 2,
                            b = 1,
                            re = 4,
                            nir = 5,
                            crop_to_shape_ext = TRUE,
                            grid = TRUE,
                            nrow = 1,
                            ncol = 1,
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
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  nlyrs <- terra::nlyr(mosaic)
  if(build_shapefile){
    if(verbose){
      cat("\014","\nBuilding the mosaic...\n")
    }
    if(is.null(basemap)){
      basemap <- mosaic_view(mosaic,
                             r = r,
                             g = g,
                             b = b,
                             re = re,
                             nir = nir,
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
    cpoints <- sf::st_sf(data.frame(id = 1),
                         geometry = sf::st_sfc(polygon),
                         crs = sf::st_crs(mosaic))
  }

  # crop to the analyzed area
  if(crop_to_shape_ext){
    poly_ext <-
      cpoints |>
      sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
      terra::vect() |>
      terra::buffer(buffer_edge) |>
      terra::ext()
    mosaiccr <- terra::crop(mosaic, poly_ext)
  } else{
    mosaiccr <- mosaic
  }

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
  created_shapes <- list()
  for(k in 1:nrow(cpoints)){
    if(inherits(cpoints[k, ]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(cpoints[k, ])) == 5 & grid[[k]]){
      plot_grid <-
        make_grid(cpoints[k, ],
                  nrow = nrow[k],
                  ncol = ncol[k],
                  mosaic = mosaic,
                  buffer_col = buffer_col[k],
                  buffer_row = buffer_row[k])
    } else{
      plot_grid <-  cpoints[k, ] |> sf::st_transform(sf::st_crs(mosaic)) |> poorman::select(geometry)
    }
    created_shapes[[k]] <- plot_grid

  }
  if(check_shapefile){
    if(verbose){
      cat("\014","\nChecking the built shapefile...\n")
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
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- terra::aggregate(mosaiccr, fact = downsample)
    }
    if(build_shapefile){
      mapview::mapview() |> mapedit::editMap()
    }
    edited <-
      mapedit::editFeatures(pg_edit, basemap) |>
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
#' This function analyzes a mosaic of remote sensing data (UVAs or satellite
#' imagery), extracting information from specified regions of interest (ROIs)
#' defined in a shapefile or interactively drawn on the mosaic. It allows
#' counting and measuring individuals (eg., plants), computing canopy coverage,
#' and statistical summaries (eg., mean, coefficient of variation) for
#' vegetation indices (eg, NDVI) at a block, plot, individual levels or even
#' extract the raw results at pixel level.
#'
#' @details
#' Since multiple blocks can be analyzed, the length of arguments `grid`,
#' `nrow`, `ncol`, `buffer_edge`, , `buffer_col`, `buffer_row`, `segment_plot`,
#' `segment_i, ndividuals`, `includ_if`, `threshold`, `segment_index`, `invert`,
#' `filter`, `threshold`, `lower_size`, `upper_size`, `watershed`, and
#' `lower_noise`, can be either an scalar (the same argument applied to all the
#' drawn blocks), or a vector with the same length as the number of drawn. In
#' the last, each block can be analyzed with different arguments.
#'
#' When `segment_individuals = TRUE` is enabled, individuals are included within
#' each plot based on the `include_if` argument. The default value
#' (`'centroid'`) includes an object in a given plot if the centroid of that
#' object is within the plot. This makes the inclusion mutually exclusive (i.e.,
#' an individual is included in only one plot). If `'covered'` is selected,
#' objects are included only if their entire area is covered by the plot. On the
#' other hand, selecting `overlap` is the complement of `covered`; in other
#' words, objects that overlap the plot boundary are included. Finally, when
#' `intersect` is chosen, objects that intersect the plot boundary are included.
#' This makes the inclusion ambiguous (i.e., an object can be included in more
#' than one plot).

#'
#' @inheritParams mosaic_view
#' @inheritParams analyze_objects
#' @param crop_to_shape_ext Crop the mosaic to the extension of shapefile?
#'   Defaults to `TRUE`. This allows for a faster index computation when the
#'   region of the built shapefile is much smaller than the entire mosaic
#'   extension.
#' @param grid Logical, indicating whether to use a grid for segmentation
#'   (default: TRUE).
#' @param nrow Number of rows for the grid (default: 1).
#' @param ncol Number of columns for the grid (default: 1).
#' @param indexes An optional `SpatRaster` object with the image indexes,
#'   computed with [mosaic_index()].
#' @param shapefile An optional shapefile containing regions of interest (ROIs)
#'   for analysis.
#' @param basemap An optional basemap generated with [mosaic_view()].
#' @param build_shapefile Logical, indicating whether to interactively draw ROIs
#'   if the shapefile is `NULL` (default: TRUE).
#' @param check_shapefile Logical, indicating whether to validate the shapefile
#'   with an interactive map view (default: TRUE). This enables live editing of
#'   the drawn shapefile by deleting or changing the drawn grids.
#' @param buffer_edge Width of the buffer around the shapefile (default: 5).
#' @param buffer_col,buffer_row Buffering factor for the columns and rows,
#'   respectively, of each individual plot's side. A value between 0 and 0.5
#'   where 0 means no buffering and 0.5 means complete buffering (default: 0). A
#'   value of 0.25 will buffer the plot by 25% on each side.
#' @param segment_plot Logical, indicating whether to segment plots (default:
#'   FALSE). If `TRUE`, the `segment_index` will be computed, and pixels with
#'   values below the `threshold` will be selected.
#' @param segment_individuals Logical, indicating whether to segment individuals
#'   within plots (default: FALSE). If `TRUE`, the `segment_index` will be
#'   computed, and pixels with values below the `threshold` will be selected, and
#'   a watershed-based segmentation will be performed.
#' @param segment_pick When `segment_plot` or `segment_individuals` are `TRUE`,
#'   `segment_pick` allows segmenting background (eg., soil) and foreground
#'   (eg., plants) interactively by picking samples from background and
#'   foreground using [mosaic_segment_pick()]
#' @param map_individuals If `TRUE`, the distance between objects within plots
#'   is computed. The distance can be mapped either in the horizontal or vertical
#'   direction. The distances, coefficient of variation (CV), and mean of
#'   distances are then returned.
#' @param map_direction The direction for mapping individuals within plots.
#'   Should be one of `"horizontal"` or `"vertical"` (default).
#' @param watershed If `TRUE` (default), performs watershed-based object
#'   detection. This will detect objects even when they are touching one another.
#'   If FALSE, all pixels for each connected set of foreground pixels are set to
#'   a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. A higher value smooths out small objects.
#' @param include_if Character vector specifying the type of intersection.
#'   Defaults to "centroid" (individuals in which the centroid is included within
#'   the drawn plot will be included in that plot). Other possible values include
#'   `"covered"`, `"overlap"`, and `"intersect"`. See Details for a detailed
#'   explanation of these intersecting controls.
#' @param plot_index The index(es) to be computed for the drawn plots. Either a
#'   single vegetation index (e.g., `"GLAI"`), a vector of indexes (e.g.,
#'   `c("GLAI", "NGRDI", "HUE")`), or a custom index based on the available
#'   bands (e.g., `"(R-B)/(R+B)"`). See [pliman_indexes()] and [image_index()]
#'   for more details.
#' @param segment_index The index used for segmentation. The same rule as
#'   `plot_index`. Defaults to `NULL`
#' @param threshold By default (threshold = "Otsu"), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is provided, this value will be used as a threshold.
#' @param filter Performs median filtering in the binary image? See more at
#'   image_filter(). Defaults to FALSE. Use a positive integer to define the
#'   size of the median filtering. Larger values are effective at removing
#'   noise but adversely affect edges.
#' @param summarize_fun The function to compute summaries for the pixel values.
#'   Defaults to "mean," i.e., the mean value of the pixels (either at a plot- or
#'   individual-level) is returned.
#' @param summarize_quantiles quantiles to be computed when 'quantile' is on `summarize_fun`.
#' @param attribute The attribute to be shown at the plot when `plot` is `TRUE`. Defaults to the first `summary_fun` and first `segment_index`.
#' @param invert Logical, indicating whether to invert the mask. Defaults to
#'   `FALSE`, i.e., pixels with intensity greater than the threshold values are
#'   selected.
#' @param color_regions The color palette for regions (default:
#'   rev(grDevices::terrain.colors(50))).
#' @param plot Logical, indicating whether to generate plots (default: TRUE).
#' @param verbose Logical, indicating whether to display verbose output
#'   (default: TRUE).
#'
#' @return A list containing the following objects:
#' * `result_plot`: The results at a plot level.
#' *  `result_plot_summ`: The summary of results at a plot level. When
#'  `segment_individuals = TRUE`, the number of individuals, canopy coverage,
#'  and mean values of some shape statistics such as perimeter, length, width,
#'  and diameter are computed.
#' * `result_individ`: The results at an individual level.
#' * `map_plot`: An object of class `mapview` showing the plot-level results.
#' * `map_individual`: An object of class `mapview` showing the individual-level
#'   results.
#' * `shapefile`: The generated shapefile, with the drawn grids/blocks.
#' @export
#'
#' @examples
#' if(interactive()){
#' library(pliman)
#' url <- "https://github.com/TiagoOlivoto/images/raw/master/pliman/rice_field/rice_ex.tif"
#' mosaic <- mosaic_input(url)
#' # Draw a polygon (top left, top right, bottom right, bottom left, top left)
#' # include 8 rice lines and one column
#'res <-
#'  mosaic_analyze(mosaic,
#'                 r = 1, g = 2, b = 3,
#'                 segment_individuals = TRUE,     # segment the individuals
#'                 segment_index = "(G-B)/(G+B-R)",# index for segmentation
#'                 filter = 4,
#'                 nrow = 8,
#'                 map_individuals = TRUE)
#'# map with individual results
#'res$map_indiv
#' }
mosaic_analyze <- function(mosaic,
                           r = 3,
                           g = 2,
                           b = 1,
                           re = 4,
                           nir = 5,
                           crop_to_shape_ext = TRUE,
                           grid = TRUE,
                           nrow = 1,
                           ncol = 1,
                           indexes = NULL,
                           shapefile = NULL,
                           basemap = NULL,
                           build_shapefile = TRUE,
                           check_shapefile = TRUE,
                           buffer_edge = 1,
                           buffer_col = 0,
                           buffer_row = 0,
                           segment_plot = FALSE,
                           segment_individuals = FALSE,
                           segment_pick = FALSE,
                           map_individuals = FALSE,
                           map_direction = c("horizontal", "vertical"),
                           watershed = TRUE,
                           tolerance = 1,
                           extension = 1,
                           include_if = "centroid",
                           plot_index = "GLI",
                           segment_index = NULL,
                           threshold = "Otsu",
                           filter = FALSE,
                           lower_noise = 0.15,
                           lower_size = NULL,
                           upper_size = NULL,
                           topn_lower = NULL,
                           topn_upper = NULL,
                           summarize_fun = "mean",
                           summarize_quantiles = NULL,
                           attribute = NULL,
                           invert = FALSE,
                           color_regions = rev(grDevices::terrain.colors(50)),
                           alpha = 1,
                           max_pixels = 2e6,
                           downsample = NULL,
                           quantiles = c(0, 1),
                           plot = TRUE,
                           verbose = TRUE){
  includeopt <- c("intersect", "covered", "overlap", "centroid")
  includeopt <- includeopt[sapply(include_if, function(x){pmatch(x, includeopt)})]
  if(is.null(plot_index) & !is.null(segment_index)){
    plot_index <- segment_index
  }
  if(!is.null(indexes)){
    plot_index <- names(indexes)
  }
  if(!is.null(plot_index) & is.null(segment_index)){
    segment_index <- plot_index[[1]]
  }
  if(any(segment_individuals) | any(segment_plot) & !is.null(plot_index) & !segment_index %in% plot_index){
    plot_index <- unique(append(plot_index, segment_index))
  }
  if(is.null(attribute)){
    attribute <- paste(summarize_fun[[1]], segment_index[[1]], sep = ".")
  }
  if(terra::crs(mosaic) == ""){
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  nlyrs <- terra::nlyr(mosaic)
  if(verbose){
    cat("\014","\nBuilding the mosaic...\n")
  }
  if(is.null(basemap)){
    basemap <-
      suppressWarnings(
        mosaic_view(mosaic,
                    r = r,
                    g = g,
                    b = b,
                    re = re,
                    nir = nir,
                    max_pixels = max_pixels,
                    verbose = verbose,
                    downsample = downsample,
                    quantiles = quantiles,
                    edit = FALSE)
      )
  }
  if(is.null(shapefile)){
    created_shapes <-
      suppressWarnings(
        shapefile_build(mosaic,
                        basemap = basemap,
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
                        sf_to_polygon =  TRUE,
                        buffer_edge = buffer_edge,
                        buffer_col = buffer_col,
                        buffer_row = buffer_row,
                        max_pixels = max_pixels,
                        verbose = verbose,
                        downsample = downsample,
                        quantiles = quantiles)
      )

    # crop to the analyzed area
    if(crop_to_shape_ext){
      poly_ext <-
        do.call(rbind, lapply(created_shapes, function(x){
          x
        })) |>
        sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
        terra::vect() |>
        terra::buffer(buffer_edge) |>
        terra::ext()

      mosaiccr <- terra::crop(mosaic, poly_ext)
    } else{
      mosaiccr <- mosaic
    }

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
    if(length(watershed) == 1 & length(created_shapes) != 1){
      watershed <- rep(watershed, length(created_shapes))
    }
    if(length(watershed) != length(created_shapes)){
      warning(paste0("`watershed` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
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
    if(length(lower_noise) == 1 & length(created_shapes) != 1){
      lower_noise <- rep(lower_noise, length(created_shapes))
    }
    if(length(lower_noise) != length(created_shapes)){
      warning(paste0("`lower_noise` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(is.null(lower_size) | length(lower_size) == 1 & length(created_shapes) != 1){
      lower_size <- rep(lower_size, length(created_shapes))
    }
    if(!is.null(lower_size) & length(lower_size) != length(created_shapes)){
      warning(paste0("`lower_size` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(is.null(upper_size) | length(upper_size) == 1 & length(created_shapes) != 1){
      upper_size <- rep(upper_size, length(created_shapes))
    }
    if(!is.null(upper_size) & length(upper_size) != length(created_shapes)){
      warning(paste0("`upper_size` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(is.null(topn_lower) | length(topn_lower) == 1 & length(created_shapes) != 1){
      topn_lower <- rep(topn_lower, length(created_shapes))
    }
    if(!is.null(topn_lower) & length(topn_lower) != length(created_shapes)){
      warning(paste0("`topn_lower` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
    if(is.null(topn_upper) | length(topn_upper) == 1 & length(created_shapes) != 1){
      topn_upper <- rep(topn_upper, length(created_shapes))
    }
    if(!is.null(topn_upper) & length(topn_upper) != length(created_shapes)){
      warning(paste0("`topn_upper` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
    }
  } else{
    if(inherits(shapefile, "list")){
      shapefile <- shapefile_input(shapefile |> sf_to_polygon(), info = FALSE)
    } else{
      if(inherits(shapefile, "SpatVector")){
        shapefile <- sf::st_as_sf(shapefile)
      }
      shapefile <- shapefile |> sf_to_polygon()
    }
    extm <- terra::ext(shapefile)
    xmin <- extm[1]
    xmax <- extm[2]
    ymin <- extm[3]
    ymax <- extm[4]
    coords <- matrix(c(xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin, xmin, ymax), ncol = 2, byrow = TRUE)

    geoms <- sf::st_as_sf(shapefile$geometry)
    sf::st_geometry(geoms) <- "geometry"
    created_shapes <- list(geoms)

    # crop to the analyzed area
    if(crop_to_shape_ext){
      poly_ext <-
        created_shapes[[1]] |>
        sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
        terra::vect() |>
        terra::buffer(buffer_edge) |>
        terra::ext()
      mosaiccr <- terra::crop(mosaic, poly_ext)
    } else{
      mosaiccr <- mosaic
    }
  }

  # return(created_shapes)
  # compute the indexes

  if(is.null(indexes)){
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
                           nir = nir,
                           plot = FALSE)
            })
        )
      )
    } else{
      plot_index <- names(mosaiccr)
      mind <- mosaiccr
    }
  } else{
    mind <- indexes
    if(!segment_index %in% names(mind)){
      stop("`segment_index` must be present in `indexes`")
    }
  }

  results <- list()
  result_indiv <- list()
  extends <- terra::ext(mosaiccr)
  usepickmask <- segment_pick & (segment_individuals | segment_plot)
  if(usepickmask){
    if(build_shapefile & is.null(shapefile)){
      mapview::mapview() |> mapedit::editMap()
    }
    mask <- suppressWarnings(
      mosaic_segment_pick(mosaic,
                          basemap = basemap,
                          r = r,
                          g = g,
                          b = b,
                          re = re,
                          nir = nir,
                          max_pixels = max_pixels,
                          return = "mask")
    )
  }
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
      if(crop_to_shape_ext){
        ext_anal <-
          plot_grid |>
          sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
        mind_temp <- terra::crop(mind, terra::ext(ext_anal))
      } else{
        mind_temp <- mind
      }
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(usepickmask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            stop("`segment_index` must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        # compute plot coverage
        dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        dmask <- EBImage::bwlabel(dmask)
        conts <- EBImage::ocontour(dmask)
        conts <- conts[sapply(conts, nrow) > 2]
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        sf_plt <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extends[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extends[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        ) |>
          sf::st_make_valid()
        covered_area <-
          suppressWarnings(
            sapply(1:nrow(plot_grid), function(i){
              plot_grid[i, ] |>
                sf::st_intersection(sf_plt) |>
                sf::st_area() |>
                sum()
            })
          )
        plot_grid <-
          plot_grid |>
          poorman::mutate(covered_area = as.numeric(covered_area),
                          plot_area = as.numeric(sf::st_area(geometry)),
                          coverage = covered_area / plot_area)

        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE, inverse = !invert[j])
      }
      # check if segmentation is performed (analyze individuals)
      if(segment_individuals[j]){
        if(usepickmask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }

        } else{
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        dmask <- EBImage::Image(matrix(mask, ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(watershed[j]){
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
        plot_id <- leading_zeros(as.numeric(plot_id[valid_rows, ]), n = 4)
        addmeasures <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  }))
        gridindiv <- cbind(sf_df, plot_id, addmeasures)[c(2, 1, 3:10)]
        # control noise removing
        if(!is.null(lower_size[j]) & !is.null(topn_lower[j]) | !is.null(upper_size[j]) & !is.null(topn_upper[j])){
          stop("Only one of 'lower_*' or 'topn_*' can be used.")
        }
        ifelse(!is.null(lower_size[j]),
               gridindiv <- gridindiv[gridindiv$area > lower_size[j], ],
               gridindiv <- gridindiv[gridindiv$area > mean(gridindiv$area) * lower_noise[j], ])

        if(!is.null(upper_size[j])){
          gridindiv <- gridindiv[gridindiv$area < upper_size[j], ]
        }
        if(!is.null(topn_lower[j])){
          gridindiv <- gridindiv[order(gridindiv$area),][1:topn_lower[j],]
        }
        if(!is.null(topn_upper[j])){
          gridindiv <- gridindiv[order(gridindiv$area, decreasing = TRUE),][1:topn_upper[j],]
        }

        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = gridindiv,
                                       fun = summarize_fun,
                                       quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          valindiv <-
            do.call(rbind, lapply(1:length(valindiv), function(i){
              tmp <- transform(valindiv[[i]],
                               individual = paste0(i),
                               block = paste0("B", leading_zeros(j, n = 2)))
              tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

            }
            ))
          if("coverage_fraction" %in% colnames(valindiv)){
            valindiv$coverage_fraction <- NULL
          }
        }
        if(length(plot_index) == 1){
          colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
        }
        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), gridindiv, valindiv, check.names = FALSE)
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
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

    } else{
      ####### ANY TYPE OF POLYGON ########
      # check if segmentation is performed
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      if(crop_to_shape_ext){
        ext_anal <-
          plot_grid |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
        mind_temp <- terra::crop(mind, terra::ext(ext_anal))
      } else{
        mind_temp <- mind
      }
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(usepickmask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            stop("`segment_index` must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        # compute plot coverage
        dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        dmask[is.na(dmask) == TRUE] <- 1
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        dmask <- EBImage::bwlabel(dmask)
        conts <- EBImage::ocontour(dmask)
        conts <- conts[sapply(conts, nrow) > 2]
        resx <- terra::res(mosaiccr)[1]
        resy <- terra::res(mosaiccr)[1]
        sf_plt <- sf::st_sf(
          geometry = lapply(conts, function(x) {
            tmp <- x
            tmp[, 2] <-  extends[3] + (nrow(mask) - tmp[, 2]) * resy
            tmp[, 1] <- extends[1] + tmp[, 1] * resy
            geometry = sf::st_polygon(list(as.matrix(tmp |> poly_close())))
          }),
          data = data.frame(individual = paste0(1:length(conts))),
          crs = terra::crs(mosaic)
        ) |>
          sf::st_make_valid()
        covered_area <-
          suppressWarnings(
            sapply(1:nrow(plot_grid), function(i){
              plot_grid[i, ] |>
                sf::st_intersection(sf_plt) |>
                sf::st_area() |>
                sum()
            })
          )
        plot_grid <-
          plot_grid |>
          poorman::mutate(covered_area = as.numeric(covered_area),
                          plot_area = as.numeric(sf::st_area(geometry)),
                          coverage = covered_area / plot_area)

        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE, inverse = !invert[j])
      }

      if(segment_individuals[j]){
        if(usepickmask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] < thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] > thresh
          }
        }
        dmask <- EBImage::Image(matrix(matrix(mask), ncol = nrow(mind_temp), nrow = ncol(mind_temp)))
        extends <- terra::ext(mind_temp)
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

        # control noise removing
        if(!is.null(lower_size[j]) & !is.null(topn_lower[j]) | !is.null(upper_size[j]) & !is.null(topn_upper[j])){
          stop("Only one of 'lower_*' or 'topn_*' can be used.")
        }
        ifelse(!is.null(lower_size[j]),
               gridindiv <- gridindiv[gridindiv$area > lower_size[j], ],
               gridindiv <- gridindiv[gridindiv$area > mean(gridindiv$area) * lower_noise[j], ])
        if(!is.null(upper_size[j])){
          gridindiv <- gridindiv[gridindiv$area < upper_size[j], ]
        }
        if(!is.null(topn_lower[j])){
          gridindiv <- gridindiv[order(gridindiv$area),][1:topn_lower[j],]
        }
        if(!is.null(topn_upper[j])){
          gridindiv <- gridindiv[order(gridindiv$area, decreasing = TRUE),][1:topn_upper[j],]
        }

        valindiv <-
          exactextractr::exact_extract(x = mind_temp,
                                       y = gridindiv,
                                       fun = summarize_fun,
                                       quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          valindiv <-
            do.call(rbind, lapply(1:length(valindiv), function(i){
              tmp <- transform(valindiv[[i]],
                               individual = paste0(i),
                               block = paste0("B", leading_zeros(j, n = 2)))
              tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]
            }
            ))
          if("coverage_fraction" %in% colnames(valindiv)){
            valindiv$coverage_fraction <- NULL
          }
        }
        if(length(plot_index) == 1){
          colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
        }

        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), plot_id = leading_zeros(1, n = 4), gridindiv, valindiv, check.names = FALSE)
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
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))
    }

    # bind the results
    if(inherits(vals, "list")){
      vals <-
        do.call(rbind, lapply(1:length(vals), function(i){
          tmp <- transform(vals[[i]],
                           plot_id = paste0(leading_zeros(i, n = 4)),
                           block = paste0("B", leading_zeros(j, 2)))
          tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]
        }
        ))
      if("coverage_fraction" %in% colnames(vals)){
        vals$coverage_fraction <- NULL
      }
      if(length(plot_index) == 1){
        colnames(vals) <- paste0(colnames(vals), ".", plot_index)
      }
      vals <-
        vals |>
        poorman::nest_by(block, plot_id) |>
        poorman::ungroup()
      vals <- cbind(plot_grid, vals)
    } else{
      if(length(plot_index) == 1){
        colnames(vals) <- paste0(colnames(vals), ".", plot_index)
      }
      vals <- transform(vals,
                        plot_id = paste0(leading_zeros(1:nrow(vals), n = 4)),
                        block = paste0("B", leading_zeros(j, 2)))
      vals <- vals[, c(ncol(vals), ncol(vals) - 1, 1:(ncol(vals) - 2))]
    }
    if(!is.null(summarize_fun)){
      results[[j]] <- cbind(plot_grid, vals, check.names = FALSE)
    } else{
      results[[j]] <- vals
    }
    # end
  }


  # bind the results  ## at a level plot
  results <-
    do.call(rbind, lapply(results, function(x){x})) |>
    poorman::relocate(block, plot_id, .before = 1)


  if(any(segment_individuals)){

    result_indiv <- do.call(rbind, lapply(result_indiv, function(x){x}))
    blockid <- unique(result_indiv$block)
    summres <-
      lapply(1:length(blockid), function(i){
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
      as.data.frame(check.names = FALSE) |>
      poorman::select(block, plot_id, plot_area)
    # compute coverage area
    result_plot_summ <-
      do.call(rbind, lapply(summres, function(x){x})) |>
      poorman::left_join(plot_area, by = c("block", "plot_id")) |>
      poorman::mutate(coverage = as.numeric(area_sum / plot_area), .after = area) |>
      poorman::left_join(results |>   poorman::select(block, plot_id, geometry), by = c("block", "plot_id")) |>
      sf::st_as_sf()
    colnames(result_plot_summ)[12:(ncol(result_plot_summ) - 3)] <- colnames(result_indiv)[11:(ncol(result_indiv) - 1)]

    centroid <-
      suppressWarnings(sf::st_centroid(result_indiv) |>
                         sf::st_coordinates() |>
                         as.data.frame() |>
                         setNames(c("x", "y")))
    result_indiv <-
      poorman::bind_cols(result_indiv, centroid) |>
      poorman::relocate(x, y, .after = individual)


    if(map_individuals){
      result_individ_map <- map_individuals(result_indiv, direction = map_direction)

      result_plot_summ <-
        result_plot_summ |>
        poorman::mutate(mean_distance = result_individ_map$means,
                        cv = result_individ_map$cvs,
                        .before = n)
    } else{
      result_individ_map <- NULL
    }


  } else{
    result_plot_summ <- NULL
    result_indiv <- NULL
    result_individ_map <- NULL

  }

  if(!is.null(summarize_fun) & isTRUE(plot)){
    if(verbose){
      cat("\014","\nPreparing to plot...\n")
    }
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- terra::aggregate(mosaiccr, fact = downsample)
    }
    # if(nlyrs < 3){
    #   basemap <-
    #     mapview::mapview(mosaiccr,
    #                      maxpixels = 5e6,
    #                      legend = FALSE,
    #                      map.types = "CartoDB.Positron",
    #                      alpha.regions = 1,
    #                      na.color = "transparent",
    #                      verbose = FALSE)
    # } else{
    #   basemap <-
    #     mapview::viewRGB(
    #       as(mosaiccr, "Raster"),
    #       layer.name = "base",
    #       r = r,
    #       g = g,
    #       b = b,
    #       na.color = "#00000000",
    #       maxpixels = 60000000,
    #       quantiles = quantiles
    #     )
    # }
    if(any(segment_individuals)){
      dfplot <- result_plot_summ
    } else{
      dfplot <- results
    }
    map <-
      basemap +
      suppressWarnings(
        mapview::mapview(dfplot,
                         zcol = attribute,
                         layer.name = attribute,
                         col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                         alpha.regions = 0.75,
                         na.color = "#00000000",
                         maxBytes = 64 * 1024 * 1024,
                         verbose = FALSE)
      )

    if(any(segment_individuals)){
      attribute <- ifelse(!attribute %in% colnames(result_indiv), "area", attribute)
      mapindivid <-
        basemap +
        suppressWarnings(
          mapview::mapview(result_plot_summ,
                           legend = FALSE,
                           alpha.regions = 0.4,
                           zcol = "block",
                           map.types = "OpenStreetMap") +
            mapview::mapview(result_indiv,
                             zcol = attribute,
                             layer.name = attribute,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
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
              result_individ_map = result_individ_map,
              map_plot = map,
              map_indiv = mapindivid,
              shapefile = created_shapes))
}

#' Analyze mosaics iteratively
#'
#' High-resolution mosaics can take a significant amount of time to analyze,
#' especially when `segment_individuals = TRUE` is used in mosaic_analyze().
#' This is because the function needs to create in-memory arrays to segment
#' individual using the watershed algorithm. This process utilizes a for-loop
#' approach, iteratively analyzing each shape within the mosaic one at a time.
#' To speed up processing, the function crops the original mosaic to the extent
#' of the current shape before analyzing it. This reduces the resolution for
#' that specific analysis, sacrificing some detail for faster processing.
#'
#' @inheritParams mosaic_analyze
#' @param ... Further arguments passed on to [mosaic_analyze()]
#'
#' @return A list containing the following objects:
#' * `result_plot`: The results at a plot level.
#' *  `result_plot_summ`: The summary of results at a plot level. When
#'  `segment_individuals = TRUE`, the number of individuals, canopy coverage,
#'  and mean values of some shape statistics such as perimeter, length, width,
#'  and diameter are computed.
#' * `result_individ`: The results at an individual level.
#' * `map_plot`: An object of class `mapview` showing the plot-level results.
#' * `map_individual`: An object of class `mapview` showing the individual-level
#'   results.
#' @export
#'
mosaic_analyze_iter <- function(mosaic,
                                shapefile,
                                r = 3,
                                g = 2,
                                b = 1,
                                plot = TRUE,
                                verbose = TRUE,
                                max_pixels = 3e6,
                                attribute = "area",
                                segment_individuals = FALSE,
                                segment_index = "VARI",
                                plot_index =  "VARI",
                                color_regions = rev(grDevices::terrain.colors(50)),
                                alpha = 0.75,
                                quantiles = c(0, 1),
                                ...){
  bind <- list()
  ind <- mosaic_index(mosaic,
                      index = unique(c(plot_index, segment_index)),
                      r = r,
                      g = g,
                      b = b,
                      re = re,
                      nir = nir,
                      plot = FALSE)
  for (i in 1:nrow(shapefile)) {
    if(verbose){
      cat("\014","\nAnalyzing plot", i, "\n")
    }
    bind[[paste0("P", leading_zeros(i, 4))]] <-
      mosaic_analyze(terra::crop(mosaic, terra::vect(shapefile$geometry[[i]]) |> terra::ext()),
                     r = r, g = g, b = b,
                     indexes = terra::crop(ind, terra::vect(shapefile$geometry[[i]]) |> terra::ext()),
                     shapefile = shapefile[i, ],
                     segment_individuals = segment_individuals,
                     segment_index = segment_index,
                     plot_index = plot_index,
                     build_shapefile = FALSE,
                     plot = FALSE,
                     grid = FALSE,
                     verbose = FALSE,
                     crop_to_shape_ext = FALSE,
                     ...)
  }
  if(is.null(bind[[1]]$result_individ_map)){
    result_individ_map <- NULL
  }
  if(is.null(bind[[1]]$result_indiv)){
    result_indiv <- result_plot_summ <- NULL
  } else{
    result_indiv <- poorman::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_indiv
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      poorman::relocate(plot_id, .after = block) |>
      sf::st_as_sf()

    result_plot_summ <- poorman::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_plot_summ
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      poorman::relocate(plot_id, .after = block) |>
      sf::st_as_sf()
  }

  result_plot <- poorman::bind_rows(
    lapply(bind, function(x){
      tmp <- x$result_plot
      tmp$plot_id <- NULL
      tmp
    }),
    .id = "plot_id"
  ) |>
    poorman::relocate(plot_id, .after = block) |>
    sf::st_as_sf()



  if(isTRUE(plot)){
    if(verbose){
      cat("\014","\nPreparing to plot...\n")
    }
    if(terra::nlyr(mosaic) < 3){
      basemap <-
        suppressWarnings(
          mapview::mapview(mosaic,
                           maxpixels = 5e6,
                           legend = FALSE,
                           map.types = "CartoDB.Positron",
                           alpha.regions = 1,
                           na.color = "transparent",
                           verbose = FALSE)
        )
    } else{
      basemap <-
        suppressWarnings(
          mapview::viewRGB(
            as(mosaic, "Raster"),
            layer.name = "base",
            r = r,
            g = g,
            b = b,
            na.color = "#00000000",
            maxpixels = 5e6,
            quantiles = quantiles
          )
        )
    }
    if(!is.null(result_indiv)){
      dfplot <- result_plot_summ
    } else{
      dfplot <- result_plot
    }
    # plot level
    map <-
      basemap +
      suppressWarnings(
        mapview::mapview(dfplot,
                         zcol = attribute,
                         layer.name = attribute,
                         col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                         alpha.regions = 0.75,
                         na.color = "#00000000",
                         maxBytes = 64 * 1024 * 1024,
                         verbose = FALSE)
      )
    # individual plot
    if(!is.null(result_indiv)){
      mapindivid <-
        basemap +
        suppressWarnings(
          mapview::mapview(result_plot_summ,
                           alpha.regions = 0.4,
                           zcol = attribute,
                           col.regions = custom_palette(c("darkred", "yellow", "darkgreen"), n = 3),
                           map.types = "OpenStreetMap") +
            mapview::mapview(result_indiv,
                             zcol = ifelse(!attribute %in% colnames(result_indiv), "area", attribute),
                             layer.name = ifelse(!attribute %in% colnames(result_indiv), "area", attribute),
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
    } else{
      mapindivid <- NULL
    }
  } else{
    map <- NULL
    mapindivid <- NULL
  }
  if(verbose){
    cat("\014","\nDone", i, "\n")
  }
  return(list(result_plot = result_plot,
              result_plot_summ = result_plot_summ,
              result_indiv = result_indiv,
              result_individ_map = result_individ_map,
              map_plot = map,
              map_indiv = mapindivid))
}



#' Mosaic View
#'
#' @details
#' The function can generate either an interactive map using the 'mapview'
#' package or a static plot using the 'base' package, depending on the `viewer`
#' and `show` parameters. If show = "index" is used, the function first computes
#' an image index that can be either an RGB-based index or a multispectral
#' index, if a multispectral mosaic is provided.
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
#' @param shapefile An optional shapefile of class `sf` to be plotted over the
#'   mosaic. It can be, for example, a plot-level result returned by
#'   [mosaic_analyze()].
#' @param attribute The attribute name(s) or column number(s) in shapefile table
#'   of the column(s) to be rendered.
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
#' @importFrom sf st_crs st_transform st_make_grid st_intersection st_make_valid
#' @importFrom poorman summarise across mutate arrange left_join bind_cols
#'   bind_rows contains ends_with everything between pivot_longer pivot_wider
#'   where select filter relocate rename
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
                        shapefile = NULL,
                        attribute = NULL,
                        viewer = c("mapview", "base"),
                        show = c("rgb", "index"),
                        index = "B",
                        max_pixels = 1000000,
                        downsample = NULL,
                        alpha = 1,
                        quantiles = c(0, 1),
                        color_regions = custom_palette(c("red", "yellow", "forestgreen")),
                        axes = FALSE,
                        ...){
  terra::terraOptions(progress = 0)
  on.exit(terra::terraOptions(progress = 1))
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
  dwspf <- find_aggrfact(mosaic, max_pixels = max_pixels)
  if(dwspf > 0 & is.null(downsample)){
    message(paste0("Using downsample = ", dwspf, " so that the number of rendered pixels approximates the `max_pixels`"))
    mosaic <- terra::aggregate(mosaic, fact = dwspf)
  }
  if(viewopt == "index" & terra::nlyr(mosaic) > 2){
    mosaic <- mosaic_index(mosaic, index = index, plot = FALSE)
  }
  if(viewopt == "rgb"){
    if(terra::nlyr(mosaic) > 2){
      if(is.null(shapefile)){
        map <-
          mapview::viewRGB(as(mosaic, "Raster"),
                           na.color = "#00000000",
                           layer.name = "base",
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
        mapview::viewRGB(as(mosaic, "Raster"),
                         na.color = "#00000000",
                         layer.name = "base",
                         r = r,
                         g = g,
                         b = b,
                         maxpixels = 60000000,
                         quantiles = quantiles) +
          suppressWarnings(
            mapview::mapview(shapefile,
                             zcol = attribute,
                             layer.name = attribute,
                             col.regions = color_regions,
                             alpha.regions = alpha,
                             na.color = "#00000000",
                             maxBytes = 64 * 1024 * 1024,
                             verbose = FALSE))
      }

    } else{
      if(vieweropt == "base"){
        terra::plot(mosaic,
                    axes = axes,
                    colNA = "white",
                    ...)
      } else{
        if(is.null(shapefile)){
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
        } else{
          mapview::mapview(as(mosaic, "Raster"),
                           na.color = "#00000000",
                           layer.name = "",
                           maxpixels = 60000000) +
            suppressWarnings(
              mapview::mapview(shapefile,
                               zcol = attribute,
                               layer.name = attribute,
                               col.regions = color_regions,
                               alpha.regions = alpha,
                               na.color = "#00000000",
                               maxBytes = 64 * 1024 * 1024,
                               verbose = FALSE))

        }
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
#' @param check_16bits Checks if mosaic has maximum value in the 16-bits format
#'   (65535), and replaces it by NA. Defaults to `TRUE`.
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
mosaic_input <- function(mosaic,
                         info = TRUE,
                         check_16bits = TRUE,
                         ...){
  mosaic <- suppressWarnings(terra::rast(mosaic, ...))
  if(terra::crs(mosaic) == ""){
    message("Missing Coordinate Reference System. Setting to EPSG:3857")
    terra::crs(mosaic) <- terra::crs("EPSG:3857")
  }
  if(info){
    print(mosaic)
  }
  if(check_16bits){
    if(max(suppressWarnings(terra::minmax(mosaic)), na.rm = TRUE) == 65535){
      mosaic[mosaic == 65535] <- NA
    }
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


#' A wrapper around terra::resample()
#'
#' Transfers values between SpatRaster objects that do not align (have a
#' different origin and/or resolution). See [terra::resample()] for more details
#'
#' @param mosaic SpatRaster to be resampled
#' @param y SpatRaster with the geometry that x should be resampled to
#' @param ... Further arguments passed on to [terra::resample()].
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' library(pliman)
#' library(terra)
#' r <- rast(nrows=3, ncols=3, xmin=0, xmax=10, ymin=0, ymax=10)
#' values(r) <- 1:ncell(r)
#' s <- rast(nrows=25, ncols=30, xmin=1, xmax=11, ymin=-1, ymax=11)
#' x <- mosaic_resample(r, s, method="bilinear")
#' opar <- par(no.readonly =TRUE)
#' par(mfrow=c(1,2))
#' plot(r)
#' plot(x)
#' par(opar)
mosaic_resample <- function(mosaic, y, ...){
  terra::resample(mosaic, y, ...)
}


#' A wrapper around terra::aggregate()
#'
#' Aggregate a SpatRaster to create a new SpatRaster with a lower resolution
#' (larger cells). See [terra::aggregate()] for more details
#'
#' @param mosaic SpatRaster
#' @param fact Aggregation factor expressed as number of cells in each direction
#'   (horizontally and vertically). See [terra::aggregate()] for more details.
#' @param ... Further arguments passed on to [terra::aggregate()].
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' library(pliman)
#' library(terra)
#' r <- rast()
#' values(r) <- 1:ncell(r)
#' r2 <- mosaic_aggregate(r, fact = 10)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2))
#' mosaic_plot(r)
#' mosaic_plot(r2)
#' par(opar)
mosaic_aggregate <- function(mosaic, fact = 2, ...){
  terra::aggregate(mosaic, fact, ...)
}


#' A wrapper around terra::plot()
#'
#' Plot the values of a SpatRaster
#'
#' @param mosaic SpatRaster
#' @param ... Further arguments passed on to [terra::plot()].
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' library(pliman)
#' r <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_plot(r)
mosaic_plot <- function(mosaic, ...){
  if(!inherits(mosaic, "SpatRaster")){
    stop("'mosaic' must be an object of class 'SpatRaster'")
  }
  terra::plot(mosaic, ...)
}

#' A wrapper around terra::plotRGB()
#'
#' Plot the RGB of a SpatRaster
#'
#' @param mosaic SpatRaster
#' @param ... Further arguments passed on to [terra::plotRGB()].
#'
#' @return A `NULL` object
#' @export
#'
mosaic_plot_rgb <- function(mosaic, ...){
  if(!inherits(mosaic, "SpatRaster")){
    stop("'mosaic' must be an object of class 'SpatRaster'")
  }
  terra::plotRGB(mosaic, ...)
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
                           r = 3,
                           g = 2,
                           b = 1,
                           max_pixels = 3e6){
  shapefile <- shapefile_input(shapefile, info = FALSE)
  if(!is.null(mosaic)){
    downsample <- find_aggrfact(mosaic, max_pixels = max_pixels)
    if(downsample > 0){
      mosaic <- terra::aggregate(mosaic, fact = downsample, progress = FALSE)
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
    edited <- mapedit::editFeatures(shapefile |> sf::st_transform(crs = 4326), map)
  } else{
    edited <- mapedit::editFeatures(shapefile)
  }
  return(edited)
}


#' Crop a mosaic
#'
#' Crop a `SpatRaster` object based on user-defined selection using an
#' interactive map or plot.
#'
#' @details This function uses the `mosaic_view` function to display an
#'   interactive map or plot of the mosaic raster, allowing users to draw a
#'   rectangle to select the cropping area. The selected area is then cropped
#'   from the input mosaic and returned as a new `SpatRaster` object. If
#'   `shapefile` is declared, the mosaic will be cropped to the extent of
#'   `shapefile`.
#' @importFrom terra crs
#' @inheritParams mosaic_view
#' @param shapefile An optional `SpatVector`, that can be created with
#'   [shapefile_input()].
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
                        shapefile = NULL,
                        show = c("rgb", "index"),
                        index = "R",
                        max_pixels = 500000,
                        downsample = NULL,
                        ...){
  if(is.null(shapefile)){
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
  } else{
    cropped <- terra::crop(mosaic, shapefile)
  }
  invisible(cropped)

}


#' Mosaic Index
#'
#' Compute or extract an index layer from a multi-band mosaic raster.
#' @inheritParams mosaic_view
#' @inheritParams image_index
#' @param index A character value (or a vector of characters) specifying the
#'   target mode for conversion to a binary image. Use [pliman_indexes_rgb()]
#'   and [pliman_indexes_me()] to see the available RGB and multispectral
#'   indexes, respectively. Users can also calculate their own index using  `R,
#'   G, B, RE, and NIR` bands (eg., `index = "R+B/G"`) or using the names of the
#'   mosaic's layers (ex., "(band_1 + band_2) / 2").
#' @param plot Plot the computed index? Defaults to `TRUE`.
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
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' names(mosaic)
#' elev2 <- mosaic_index(mosaic, "elevation * 5", plot = FALSE)
#' oldpar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#'
#' mosaic_plot(mosaic)
#' mosaic_plot(elev2)
#'
#' # return the original parameters
#' par(oldpar)
#'

mosaic_index <- function(mosaic,
                         index = "R",
                         r = 3,
                         g = 2,
                         b = 1,
                         re = 4,
                         nir = 5,
                         plot = TRUE){
  if(length(index) == 1){
    if(inherits(mosaic, "Image")){
      ras <- t(terra::rast(mosaic@.Data))
    } else{
      ras <- mosaic
    }
    ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
    checkind <- index %in% ind$Index
    if (!all(checkind)) {
      message(paste("Index '", paste0(index[!checkind], collapse = ", "), "' is not available. Trying to compute your own index.",
                    sep = ""))
    }
    pattern <- "\\b\\w+\\b"
    layersused <- unlist(regmatches(index, gregexpr(pattern, index, perl = TRUE)))
    onlychar <- suppressWarnings(is.na(as.numeric(layersused)))
    layers_used <- layersused[onlychar]
    if(!any(index  %in% ind$Index) & !all(layers_used  %in% c("R", "G", "B", "RE", "NIR"))){
      # Extract individual layers based on the expression
      layers_used <- layers_used[is.na(suppressWarnings(as.numeric(layers_used)))]
      layers <-
        lapply(layers_used, function(x){
          mosaic[[x]]
        })
      names(layers) <- layers_used
      mosaic_gray <- eval(parse(text = index), envir = layers)
    } else{
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
    }
    names(mosaic_gray) <- index
    if(!is.na(terra::crs(mosaic))){
      if(terra::crs(mosaic_gray) != terra::crs(mosaic)){
        suppressWarnings(terra::crs(mosaic_gray) <- terra::crs(mosaic))
      }
    } else{
      suppressWarnings(terra::crs(mosaic_gray) <- "+proj=utm +zone=32 +datum=WGS84 +units=m")
    }
  } else{
    mosaic_gray <- terra::rast(
      Map(c,
          lapply(seq_along(unique(index)), function(i){
            mosaic_index(mosaic,
                         index = unique(index)[[i]],
                         r = r,
                         g = g,
                         b = b,
                         re = re,
                         nir = nir,
                         plot = FALSE)
          })
      )
    )
  }
  if(plot){
    terra::plot(mosaic_gray)
  }
  invisible(mosaic_gray)
}
#' Segment a mosaic
#'
#' Segment a `SpatRaster` using a computed image index. By default, values
#' greater than `threshold` are kept in the mask.
#'
#' @inheritParams mosaic_index
#' @inheritParams mosaic_analyze
#'
#'
#' @return The segmented mosaic (`SpatRaster` object)
#' @export
#'
#' @examples
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' seg <-
#' mosaic_segment(mosaic,
#'                index = "elevation",
#'                threshold = 350)
#' mosaic_plot(seg)
mosaic_segment <- function(mosaic,
                           index = "R",
                           r = 3,
                           g = 2,
                           b = 1,
                           re = 4,
                           nir = 5,
                           threshold = "Otsu",
                           invert = FALSE){
  ind <- mosaic_index(mosaic,
                      index = index,
                      r = r,
                      g = g,
                      b = b,
                      re = re,
                      nir = nir,
                      plot = FALSE)
  thresh <- ifelse(threshold == "Otsu", otsu(na.omit(terra::values(ind)[, index])), threshold)
  if(invert){
    mask <- ind[[index]] < thresh
  } else{
    mask <- ind[[index]] > thresh
  }
  terra::mask(mosaic, mask, maskvalue = TRUE, inverse = TRUE)
}


#' Segments a mosaic interactively
#'
#' The function segments a mosaic using an interative process where the user
#' picks samples from background (eg., soil) and foreground (eg., plants).
#'
#' @inheritParams mosaic_analyze
#' @param basemap An optional `mapview` object.
#' @param return The output of the function. Either 'mosaic' (the segmented
#'   mosaic), or 'mask' (the binary mask).
#'
#' @return An `SpatRaster` object with the segmented `mosaic` (if `return =
#'   'mosaic'`) or a mask (if `return = 'mask'`).
#' @export
#'
#' @examples
#' if(interactive()){
#'  mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#'  seg <- mosaic_segment_pick(mosaic)
#'  mosaic_plot(seg)
#' }
mosaic_segment_pick <- function(mosaic,
                                basemap = NULL,
                                g = 2,
                                r = 3,
                                b = 1,
                                re = 4,
                                nir = 5,
                                max_pixels = 2e6,
                                downsample = NULL,
                                quantiles = c(0, 1),
                                return = c("mosaic", "mask")){
  if(!return[[1]] %in% c("mosaic", "mask")){
    stop("'return' must be one of 'mosaic' or 'mask'.")
  }
  downsample <- ifelse(is.null(downsample), find_aggrfact(mosaic, max_pixels = max_pixels), downsample)
  if(downsample > 0){
    mosaicp <- terra::aggregate(mosaic, downsample, progress = FALSE)
  }
  if(is.null(basemap)){
    basemap <-
      mosaic_view(mosaicp,
                  r = r,
                  g = g,
                  b = b,
                  re = re,
                  nir = nir,
                  max_pixels = max_pixels,
                  downsample = downsample,
                  quantiles = quantiles)
  }
  soil <- mapedit::editMap(basemap,
                           title = "Use the 'Draw Rectangle' tool to pick up background fractions",
                           editor = "leafpm")$finished
  soil <- soil |> sf::st_transform(sf::st_crs(mosaic))
  soil_sample <-
    exactextractr::exact_extract(mosaic, soil, progress = FALSE) |>
    poorman::bind_rows() |>
    poorman::select(-coverage_fraction) |>
    poorman::mutate(class = 0)

  mapview::mapview() |> mapedit::editMap()

  plant <- mapedit::editMap(basemap,
                            title = "Use the 'Draw Rectangle' tool to pick up foreground fractions",
                            editor = "leafpm")$finished
  plant <- plant |> sf::st_transform(sf::st_crs(mosaic))
  plant_sample <-
    exactextractr::exact_extract(mosaic, plant, progress = FALSE) |>
    poorman::bind_rows() |>
    poorman::select(-coverage_fraction) |>
    poorman::mutate(class = 1)
  df_train <- poorman::bind_rows(plant_sample, soil_sample)
  if(ncol(df_train) == 2){
    names(df_train)[[1]] <- names(mosaic)
  }
  mod <- suppressWarnings(
    glm(class ~.,
        data = df_train,
        family = binomial("logit"))
  )
  mask <- terra::predict(mosaic, mod, type = "response")
  mask[mask < 0.5] <- 0
  mask[mask > 0.5] <- 1
  if(return[[1]] == 'mosaic'){
    terra::mask(mosaic, mask, maskvalue = TRUE, inverse = TRUE)
  } else{
    mask
  }
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
  if(nlr == 5){
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))[,, c(r, g, b, re, nir)]
  } else if(nlr == 3){
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))[,, c(r, g, b)]
  } else{
    mosaic <- EBImage::Image(terra::as.array(terra::trans(mosaic)))
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
                         nir = nir,
                         plot = FALSE)
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
                       scale = ifelse(scl < 255, 255, scl),
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
                                     fun = summarize_fun,
                                     quantiles = summarize_quantiles,
                                     progress = FALSE,
                                     force_df = TRUE,
                                     summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))
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
                           layer.name = "base",
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


