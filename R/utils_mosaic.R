validate_and_replicate <- function(argument, created_shapes) {
  if (length(argument) != length(created_shapes)) {
    warning(paste0("`", deparse(substitute(argument)), "` must have length 1 or ", length(created_shapes), " (the number of drawn polygons)."))
  }
  if (length(argument) == 1 & length(created_shapes) != 1) {
    argument <- rep(argument, length(created_shapes))
  }
  return(argument)
}
validate_and_replicate2 <- function(argument, created_shapes) {
  if (!is.null(argument) & (length(argument) != nrow(created_shapes))) {
    warning(paste0("`", deparse(substitute(argument)), "` must have length 1 or ", nrow(created_shapes), " (the number of drawn polygons)."), call. = FALSE)
  }
  if (length(argument) == 1 & nrow(created_shapes) != 1) {
    argument <- rep(argument, nrow(created_shapes))
  }
  return(argument)
}
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
compute_dists <- function(subset_coords, direction = c("horizontal", "vertical")){
  optdirec <- c("horizontal", "vertical")
  optdirec <- pmatch(direction[[1]], optdirec)
  n <- nrow(subset_coords)
  subset_coords <- subset_coords |> dplyr::select(x, y) |> as.data.frame()
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
  subset_distances
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
    magg <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
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
#' @inheritParams image_binary
#' @inheritParams plot_id
#' @param r,g,b,re,nir,swir,tir The red, green, blue, red-edge,  near-infrared,
#'   shortwave Infrared, and thermal infrared bands of the image, respectively.
#'   By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
#'   multispectral image is provided up to seven bands can be used to compute
#'   built-in indexes. There are no limitation of band numbers if the index is
#'   computed using the band name.
#' @param crop_to_shape_ext Crop the mosaic to the extension of shapefile?
#'   Defaults to `TRUE`. This allows for a faster index computation when the
#'   region of the built shapefile is much smaller than the entire mosaic
#'   extension.
#' @param grid Logical, indicating whether to use a grid for segmentation
#'   (default: TRUE).
#' @param nrow Number of rows for the grid (default: 1).
#' @param ncol Number of columns for the grid (default: 1).
#' @param plot_width,plot_height The width and height of the plot shape (in the
#'   mosaic unit). It is mutually exclusiv with `buffer_col` and `buffer_row`.
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
#' @param mask An optional mask (SpatRaster) to mask the mosaic.
#' @param simplify Removes vertices in polygons to form simpler shapes. The
#'   function implementation uses the Douglas–Peucker algorithm using
#'   [sf::st_simplify()] for simplification.
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
                           re = NA,
                           nir = NA,
                           swir = NA,
                           tir = NA,
                           crop_to_shape_ext = TRUE,
                           grid = TRUE,
                           nrow = 1,
                           ncol = 1,
                           plot_width = NULL,
                           plot_height = NULL,
                           layout = "lrtb",
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
                           mask = NULL,
                           simplify = FALSE,
                           map_individuals = FALSE,
                           map_direction = c("horizontal", "vertical"),
                           watershed = TRUE,
                           tolerance = 1,
                           extension = 1,
                           include_if = "centroid",
                           plot_index = "GLI",
                           segment_index = NULL,
                           threshold = "Otsu",
                           opening = FALSE,
                           closing = FALSE,
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
    terra::crs(mosaic) <- terra::crs("EPSG:4326")
    # terra::ext(mosaic) <- c(0, 1, 0, 1)
  }
  nlyrs <- terra::nlyr(mosaic)
  if(verbose){
    message("\014","\nBuilding the mosaic...\n")
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
                    swir = swir,
                    tir = tir,
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
                        grid = grid,
                        nrow = nrow,
                        ncol = ncol,
                        plot_width = plot_width,
                        plot_height = plot_height,
                        layout = layout,
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
      ress <- terra::res(mosaic)
      if(sum(ress) != 2){
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      } else{
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      }
      mosaiccr <- terra::crop(mosaic, poly_ext)
    } else{
      mosaiccr <- mosaic
    }
  } else{
    if(inherits(shapefile, "list")){
      created_shapes <- lapply(shapefile, function(x){
        # x[, "geometry"]
        x
      })
    } else{
      if(inherits(shapefile, "SpatVector")){
        created_shapes <- sf::st_as_sf(shapefile) |> sf_to_polygon()
      }
      created_shapes <- list(shapefile |> sf_to_polygon())
      names(created_shapes) <- paste(1:length(created_shapes))
    }
    if(crop_to_shape_ext){
      ress <- terra::res(mosaic)
      if(sum(ress) != 2){
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          sf::st_transform(crs = sf::st_crs(terra::crs(mosaic))) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      } else{
        poly_ext <-
          do.call(rbind, lapply(created_shapes, function(x){
            x
          })) |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()
      }

      mosaiccr <- terra::crop(mosaic, poly_ext)

    } else{
      mosaiccr <- mosaic
    }
  }
  segment_plot <- validate_and_replicate(segment_plot, created_shapes)
  segment_individuals <- validate_and_replicate(segment_individuals, created_shapes)
  threshold <- validate_and_replicate(threshold, created_shapes)
  watershed <- validate_and_replicate(watershed, created_shapes)
  segment_index <- validate_and_replicate(segment_index, created_shapes)
  invert <- validate_and_replicate(invert, created_shapes)
  includeopt <- validate_and_replicate(includeopt, created_shapes)
  opening <- validate_and_replicate(opening, created_shapes)
  closing <- validate_and_replicate(closing, created_shapes)
  filter <- validate_and_replicate(filter, created_shapes)
  grid <- validate_and_replicate(grid, created_shapes)
  lower_noise <- validate_and_replicate(lower_noise, created_shapes)

  if(!is.null(lower_size)){
    lower_size <- validate_and_replicate(lower_size, created_shapes)
  }
  if(!is.null(upper_size)){
    upper_size <- validate_and_replicate(upper_size, created_shapes)
  }
  if(!is.null(topn_lower)){
    topn_lower <- validate_and_replicate(topn_lower, created_shapes)
  }
  if(!is.null(topn_upper)){
    topn_upper <- validate_and_replicate(topn_upper, created_shapes)
  }
  #



  if(is.null(indexes)){
    if(verbose){
      message("\014","\nComputing the indexes...\n")
    }
    if(nlyrs > 1 | !all(plot_index %in% names(mosaiccr))){
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
                           swir = swir,
                           tir = tir,
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
    if(!all(segment_index %in% names(mind))){
      stop("`segment_index` must be present in `indexes`")
    }
  }

  results <- list()
  result_indiv <- list()
  extends <- terra::ext(mosaiccr)
  usepickmask <- segment_pick & (segment_individuals[[1]] | segment_plot[[1]])
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
                          max_pixels = max_pixels,
                          return = "mask")
    )
  }
  ihaveamask <- !is.null(mask) & (segment_individuals[[1]] | segment_plot[[1]])
  if(ihaveamask){
    mask <- mask
  }
  for(j in seq_along(created_shapes)){
    if(segment_plot[j] & segment_individuals[j]){
      stop("Only `segment_plot` OR `segment_individuals` can be used", call. = FALSE)
    }
    if(verbose){
      message("\014","\nExtracting data from block", j, "\n")
    }
    if(inherits(created_shapes[[j]]$geometry, "sfc_POLYGON") & nrow(sf::st_coordinates(created_shapes[[j]]$geometry[[1]])) == 5 & grid[[j]]){
      plot_grid <- created_shapes[[j]]
      sf::st_geometry(plot_grid) <- "geometry"
      if(crop_to_shape_ext){
        ress <- terra::res(mosaic)
        if(sum(ress) != 2){
          plot_grid <-
            plot_grid |>
            sf::st_transform(crs = sf::st_crs(terra::crs(mosaic)))
        }
        ext_anal <-
          plot_grid |>
          terra::vect() |>
          terra::buffer(buffer_edge) |>
          terra::ext()

        mind_temp <- terra::crop(mind, terra::ext(ext_anal))
        if(!is.null(mask)){
          mask <- terra::crop(mask, terra::ext(ext_anal))
        }
      } else{
        mind_temp <- mind
      }
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            stop("`segment_index` must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] > thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] < thresh
          }
        }
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
        # compute plot coverage

        tmp <- exactextractr::exact_extract(mind_temp,
                                            plot_grid,
                                            coverage_area = TRUE,
                                            force_df = TRUE,
                                            progress = FALSE)
        covered_area <-
          dplyr::bind_rows(
            lapply(seq_along(tmp), function(i){
              data.frame(covered_area = sum(na.omit(tmp[[i]])[, 2]),
                         plot_area = sum(tmp[[i]][, 2])) |>
                dplyr::mutate(coverage = covered_area / plot_area)
            })
          )
        plot_grid <- dplyr::bind_cols(plot_grid, covered_area)
        if(simplify){
          plot_grid <- plot_grid |> sf::st_simplify(preserveTopology = TRUE)
        }
        rm(tmp)

      }

      # check if segmentation is performed (analyze individuals)
      if(segment_individuals[j]){
        if(usepickmask | ihaveamask){
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
        if(is.numeric(opening[j]) & opening[j] > 0){
          dmask <- image_opening(dmask, size = opening[j])
        }
        if(is.numeric(closing[j]) & closing[j] > 0){
          dmask <- image_closing(dmask, size = closing[j])
        }
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
        if(simplify){
          sf_df <- sf_df |> sf::st_simplify(preserveTopology = TRUE)
        }
        centroids <- suppressWarnings(sf::st_centroid(sf_df))
        intersects <-
          switch (includeopt[j],
                  "intersect" = sf::st_intersects(sf_df, plot_grid),
                  "centroid" =  sf::st_within(centroids, plot_grid),
                  "covered" = sf::st_covered_by(sf_df, plot_grid),
                  "overlap" = sf::st_overlaps(sf_df, plot_grid),

          )
        plot_gridtmp <-
          plot_grid |>
          dplyr::mutate(plot_id_seq = paste0("P", leading_zeros(1:nrow(plot_grid), 4)))
        plot_id <- data.frame(plot_id_seq = paste0(intersects))
        valid_rows <- plot_id$plot_id_seq != "integer(0)"
        sf_df <- sf_df[valid_rows, ]
        plot_id <- paste0("P", leading_zeros(as.numeric(plot_id[valid_rows, ]), n = 4))

        gridindiv <-
          do.call(rbind,
                  lapply(1:nrow(sf_df), function(i){
                    compute_measures_mosaic(as.matrix(sf_df$geometry[[i]]))
                  })) |>
          dplyr::mutate(plot_id_seq = plot_id,
                        individual = paste0(1:nrow(sf_df)),
                        geometry = sf_df$geometry) |>
          dplyr::left_join(plot_gridtmp |> sf::st_drop_geometry(), by = dplyr::join_by(plot_id_seq)) |>
          dplyr::select(-plot_id_seq) |>
          dplyr::relocate(block, plot_id, individual, .before = 1) |>
          sf::st_sf()

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
                                       y = sf::st_sf(gridindiv),
                                       fun = summarize_fun,
                                       quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          if(is.null(summarize_fun)){
            valindiv <- dplyr::bind_rows(valindiv, .id = "individual")
            if("coverage_fraction" %in% colnames(valindiv)){
              valindiv$coverage_fraction <- NULL
            }
            if("value" %in% colnames(valindiv)){
              colnames(valindiv)[2] <- plot_index
            }
            valindiv <- valindiv |> dplyr::nest_by(individual)
          } else{
            valindiv <-
              do.call(rbind, lapply(1:length(valindiv), function(i){
                tmp <- transform(valindiv[[i]],
                                 individual = paste0(i))
                tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

              }
              ))

            if(length(plot_index) == 1){
              colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
            } else{
              # colnames(valindiv) <- c("block", "plot_id", plot_index)
              colnames(vals) <- paste0(colnames(vals), ".", plot_index)
            }
          }
        } else{
          if(length(plot_index) == 1){
            colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
          }
        }
        if(!is.null(summarize_fun)){
          valindiv <-
            dplyr::bind_cols(gridindiv, valindiv) |>
            dplyr::mutate(individual = paste0(1:nrow(gridindiv)), .before = area) |>
            sf::st_sf()
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        } else{
          valindiv <- dplyr::bind_cols(dplyr::left_join(gridindiv, valindiv, by = dplyr::join_by(individual))) |> sf::st_sf()
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        }
      } else{
        dmask <- NULL
        result_indiv[[j]] <- NULL
      }

      # extract the values for the individual plots
      # check if a mask is used and no segmentation
      if(!is.null(mask) & (!segment_individuals[[1]] & !segment_plot[[1]])){
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
      }
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
        if(!is.null(mask)){
          mask <- terra::crop(mask, terra::ext(ext_anal))
        }
      } else{
        mind_temp <- mind
      }
      extends <- terra::ext(mind_temp)
      if(segment_plot[j]){
        if(usepickmask | ihaveamask){
          if(crop_to_shape_ext){
            mask <- terra::crop(mask, terra::ext(ext_anal))
          }
        } else{
          if(!segment_index[j] %in% names(mind_temp)){
            stop("`segment_index` must be one of used in `plot_index`.")
          }
          thresh <- ifelse(threshold[j] == "Otsu", otsu(na.omit(terra::values(mind_temp)[, segment_index[j]])), threshold[j])
          if(invert[j]){
            mask <- mind_temp[[segment_index[j]]] > thresh
          } else{
            mask <- mind_temp[[segment_index[j]]] < thresh
          }
        }
        # compute plot coverage
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
        tmp <- exactextractr::exact_extract(mind_temp,
                                            plot_grid,
                                            coverage_area = TRUE,
                                            force_df = TRUE,
                                            progress = FALSE)
        covered_area <-
          dplyr::bind_rows(
            lapply(seq_along(tmp), function(i){
              data.frame(covered_area = sum(na.omit(tmp[[i]])[, 2]))
            })
          )|>
          dplyr::mutate(plot_area = as.numeric(sf::st_area(plot_grid)),
                        coverage = covered_area / plot_area)
        plot_grid <- dplyr::bind_cols(plot_grid, covered_area)
        if(simplify){
          plot_grid <- plot_grid |> sf::st_simplify(preserveTopology = TRUE)
        }
        rm(tmp)

      }

      if(segment_individuals[j]){
        if(usepickmask | ihaveamask){
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
        if(is.numeric(opening[j]) & opening[j] > 0){
          dmask <- image_opening(dmask, opening[j])
        }
        if(is.numeric(closing[j]) & closing[j] > 0){
          dmask <- image_closing(dmask, closing[j])
        }
        if(!isFALSE(filter[j]) & filter[j] > 1){
          dmask <- EBImage::medianFilter(dmask, filter[j])
        }
        if(watershed[j]){
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
        if(simplify){
          sf_df <- sf_df |> sf::st_simplify(preserveTopology = TRUE)
        }
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
                                       # quantiles = summarize_quantiles,
                                       progress = FALSE,
                                       force_df = TRUE,
                                       summarize_df = ifelse(is.function(summarize_fun), TRUE, FALSE))

        if(inherits(valindiv, "list")){
          if(is.null(summarize_fun)){
            valindiv <- dplyr::bind_rows(valindiv, .id = "individual")
            if("coverage_fraction" %in% colnames(valindiv)){
              valindiv$coverage_fraction <- NULL
            }
            if("value" %in% colnames(valindiv)){
              colnames(valindiv)[2] <- plot_index
            }
            valindiv <- valindiv |> dplyr::nest_by(individual) |> dplyr::ungroup()
          } else{
            valindiv <-
              do.call(rbind, lapply(1:length(valindiv), function(i){
                tmp <- transform(valindiv[[i]],
                                 individual = paste0(i),
                                 block = paste0("B", leading_zeros(j, n = 2)))
                tmp[, c(ncol(tmp), ncol(tmp) - 1, 1:(ncol(tmp) - 2))]

              }
              ))

            if(length(plot_index) == 1){
              colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
            } else{
              # colnames(valindiv) <- c("block", "plot_id", plot_index)
              colnames(vals) <- paste0(colnames(vals), ".", plot_index)
            }
          }
        } else{
          if(length(plot_index) == 1){
            colnames(valindiv) <- paste0(colnames(valindiv), ".", plot_index)
          }
        }

        if(!is.null(summarize_fun)){
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), plot_id = "P0001", gridindiv, valindiv, check.names = FALSE)
          result_indiv[[j]] <- valindiv
        } else{
          valindiv <- cbind(block = paste0("B", leading_zeros(j, n = 2)), plot_id = "P0001", dplyr::left_join(gridindiv, valindiv, by = dplyr::join_by(individual)), check.names = FALSE)
          result_indiv[[j]] <- valindiv[order(valindiv$plot_id), ]
        }
      } else{
        result_indiv[[j]] <- NULL
      }

      # extract the values for the individual plots
      # check if a mask is used and no segmentation
      if(!is.null(mask) & (!segment_individuals[[1]] & !segment_plot[[1]])){
        mind_temp <- terra::mask(mind_temp, mask, maskvalues = TRUE)
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
      names(vals) <- paste0(plot_grid$block, "_", plot_grid$plot_id)
      vals <- dplyr::bind_rows(vals, .id = "plot") |> pliman::separate_col(plot, c("block", "plot_id"))
      if("coverage_fraction" %in% colnames(vals)){
        vals$coverage_fraction <- NULL
      }
      if(length(plot_index) == 1){
        if(ncol(vals) == 3){
          colnames(vals)[3] <- plot_index
        }
      }
      vals <-
        vals |>
        dplyr::nest_by(block, plot_id, row, column) |>
        dplyr::ungroup() |>
        dplyr::left_join(plot_grid, by = dplyr::join_by(block, plot_id, row, column))
    } else{
      if(length(plot_index) == 1){
        if(ncol(vals) == 1){
          colnames(vals) <- paste0(colnames(vals), ".", plot_index)
        } else{
          colnames(vals) <- paste0(colnames(vals), ".", plot_index)
        }
      }
      vals <- dplyr::bind_cols(plot_grid, vals)
    }
    results[[j]] <- vals
  }


  # bind the results  ## at a level plot
  results <- dplyr::bind_rows(results) |> sf::st_sf()
  # return(list(results = results, result_indiv = result_indiv))
  if(any(segment_individuals)){
    result_indiv <- do.call(rbind, result_indiv)
    blockid <- unique(result_indiv$block)
    summres <-
      lapply(1:length(blockid), function(i){
        result_indiv |>
          dplyr::filter(block == blockid[i]) |>
          as.data.frame() |>
          dplyr::group_by(plot_id) |>
          dplyr::summarise(area_sum = sum(area, na.rm = TRUE),
                           n = length(area),
                           dplyr::across(dplyr::where(is.numeric), \(x){mean(x, na.rm = TRUE)})
          ) |>
          dplyr::mutate(block = blockid[i], .before = 1) |>
          dplyr::ungroup() |>
          dplyr::relocate(n, .after = plot_id)
      })
    names(summres) <- blockid
    # compute plot area
    plot_area <-
      results |>
      dplyr::mutate(plot_area = sf::st_area(geometry)) |>
      as.data.frame(check.names = FALSE) |>
      dplyr::select(block, plot_id, row, column, plot_area)
    # compute coverage area
    result_plot_summ <-
      do.call(rbind, lapply(summres, function(x){x})) |>
      dplyr::left_join(plot_area, by = dplyr::join_by(block, plot_id, row, column)) |>
      dplyr::mutate(coverage = as.numeric(area_sum / plot_area), .after = area) |>
      dplyr::left_join(results |>   dplyr::select(block, plot_id, row, column, geometry),
                       by = dplyr::join_by(block, plot_id, row, column)) |>
      sf::st_as_sf()

    centroid <-
      suppressWarnings(sf::st_centroid(sf::st_sf(result_indiv)) |>
                         sf::st_coordinates() |>
                         as.data.frame() |>
                         setNames(c("x", "y")))
    result_indiv <-
      dplyr::bind_cols(result_indiv, centroid) |>
      dplyr::relocate(x, y, .after = individual)

    if(map_individuals){
      dists <-
        result_indiv |>
        sf::st_drop_geometry() |>
        dplyr::select(block, plot_id, row, column, x, y) |>
        dplyr::group_by(block, plot_id, row, column)

      splits <- dplyr::group_split(dists)
      names(splits) <- dplyr::group_keys(dists) |> dplyr::mutate(key = paste0(block, "_", plot_id)) |> dplyr::pull()
      dists <- lapply(splits, compute_dists)

      cvs <- sapply(dists, function(x){
        (sd(x) / mean(x)) * 100
      })
      means <- sapply(dists, mean)

      result_plot_summ <-
        result_plot_summ |>
        dplyr::mutate(mean_distance = means,
                      cv = cvs,
                      .before = n)
      result_individ_map <- list(distances = dists,
                                 means = means,
                                 cvs = cvs)
    } else{
      result_individ_map <- NULL
    }


  } else{
    result_plot_summ <- NULL
    result_indiv <- NULL
    result_individ_map <- NULL

  }

  if(isTRUE(plot)){
    if(verbose){
      message("\014","\nPreparing to plot...\n")
    }
    downsample <- find_aggrfact(mosaiccr, max_pixels = max_pixels)
    if(downsample > 0){
      mosaiccr <- mosaic_aggregate(mosaiccr, pct = round(100 / downsample))
    }
    if(any(segment_individuals)){
      dfplot <- result_plot_summ
      if(!attribute %in% colnames(dfplot)){
        attribute <- "area"
      }
    } else{
      dfplot <- results
      if(!attribute %in% colnames(dfplot)){
        attribute <- NULL
      }
    }
    if(is.null(summarize_fun)){
      dfplot <-
        dfplot |>
        sf::st_drop_geometry() |>
        tidyr::unnest(cols = data) |>
        dplyr::group_by(block, plot_id, row, column) |>
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x){mean(x, na.rm = TRUE)}), .groups = "drop") |>
        dplyr::left_join(dfplot |> dplyr::select(block, plot_id, row, column, geometry),
                         by = dplyr::join_by(block, plot_id, row, column)) |>
        sf::st_sf()

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
    message("\014","Done!\n")
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
                                basemap = NULL,
                                r = 3,
                                g = 2,
                                b = 1,
                                re = NA,
                                nir = NA,
                                swir = NA,
                                tir = NA,
                                plot = TRUE,
                                verbose = TRUE,
                                max_pixels = 3e6,
                                attribute = NULL,
                                summarize_fun = "mean",
                                segment_plot = FALSE,
                                segment_individuals = FALSE,
                                segment_index = "VARI",
                                plot_index =  "VARI",
                                color_regions = rev(grDevices::terrain.colors(50)),
                                alpha = 0.75,
                                quantiles = c(0, 1),
                                ...){
  pind <- unique(c(plot_index, segment_index))
  if(is.null(attribute)){
    attribute <- paste(summarize_fun, pind[[1]], sep = ".")
  }
  bind <- list()
  for (i in 1:nrow(shapefile)) {
    if(verbose){
      message("\014","\nAnalyzing plot", i, "\n")
    }
    bind[[paste0("P", leading_zeros(i, 4))]] <-
      mosaic_analyze(terra::crop(mosaic, terra::vect(shapefile$geometry[[i]]) |> terra::ext()),
                     basemap = basemap,
                     r = r, g = g, b = b, re = re, nir = nir, swir = swir, tir = tir,
                     shapefile = shapefile[i, ],
                     segment_individuals = segment_individuals,
                     segment_index = segment_index,
                     segment_plot = segment_plot,
                     plot_index = pind,
                     build_shapefile = FALSE,
                     plot = FALSE,
                     grid = FALSE,
                     verbose = FALSE,
                     crop_to_shape_ext = TRUE,
                     ...)
  }
  if(is.null(bind[[1]]$result_individ_map)){
    result_individ_map <- NULL
  }
  if(is.null(bind[[1]]$result_indiv)){
    result_indiv <- result_plot_summ <- NULL
  } else{
    result_indiv <- dplyr::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_indiv
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      dplyr::relocate(plot_id, .after = block) |>
      sf::st_as_sf()

    result_plot_summ <- dplyr::bind_rows(
      lapply(bind, function(x){
        tmp <- x$result_plot_summ
        tmp$plot_id <- NULL
        tmp
      }),
      .id = "plot_id"
    ) |>
      dplyr::relocate(plot_id, .after = block) |>
      sf::st_as_sf()
  }

  result_plot <- dplyr::bind_rows(
    lapply(bind, function(x){
      tmp <- x$result_plot
      tmp$plot_id <- NULL
      tmp
    }),
    .id = "plot_id"
  ) |>
    dplyr::relocate(plot_id, .after = block) |>
    sf::st_as_sf()



  if(isTRUE(plot)){
    if(verbose){
      message("\014","\nPreparing to plot...\n")
    }
    if(is.null(basemap)){
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
                         zcol = ifelse(!attribute %in% colnames(dfplot), "area", attribute),
                         layer.name = ifelse(!attribute %in% colnames(dfplot), "area", attribute),
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
    message("\014","\nDone", i, "\n")
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
#' @inheritParams mosaic_index
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
#' @param downsample_fun The resampling function. Defaults to nearest. See further details in [mosaic_aggregate()].
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
#' @importFrom dplyr summarise across mutate arrange left_join bind_cols
#'   bind_rows contains ends_with everything between where select filter
#'   relocate rename
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
                        edit = FALSE,
                        title = "",
                        shapefile = NULL,
                        attribute = NULL,
                        viewer = c("mapview", "base"),
                        show = c("rgb", "index"),
                        index = "B",
                        max_pixels = 1000000,
                        downsample = NULL,
                        downsample_fun = "nearest",
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
    mosaic <- mosaic_aggregate(mosaic, pct = round(100 / dwspf), fun = downsample_fun)
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
#'   (65535), and replaces it by NA. Defaults to `FALSE`.
#' @param check_datatype Logical. If \code{TRUE}, checks and suggests the
#'   appropriate data type based on the raster values.
#' @param filename character. The Output filename.
#' @param datatype The datatype. By default, the function will try to guess the
#'   data type that saves more memory usage and file size. See
#'   [terra::writeRaster()] and [terra::datatype()] for more details.
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
#' x <- system.file("ex/logo.tif", package="terra")
#' rast <- mosaic_input(x)
#' mosaic_plot(rast)
#'
#' # create a temporary filename for the example
#' f <- file.path(tempdir(), "test.tif")
#' mosaic_export(rast, f, overwrite=TRUE)
#' list.files(tempdir())
#'
mosaic_input <- function(mosaic,
                         mosaic_pattern = NULL,
                         info = TRUE,
                         check_16bits = FALSE,
                         check_datatype = FALSE,
                         ...){
  if(!is.null(mosaic_pattern)){
    if(mosaic_pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      mosaic_pattern <- "^[0-9].*$"
    }
    path <- getwd()
    imgs <- list.files(pattern = mosaic_pattern, path)
    if(length(grep(mosaic_pattern, imgs)) == 0){
      stop(paste("'", mosaic_pattern, "' mosaic_pattern not found in '",
                 paste0(dir)),
           call. = FALSE)
    }
    list_img <-
      lapply(imgs, function(x){
        mosaic_input(x, info = FALSE)
      })
    names(list_img) <- imgs
    invisible(list_img)
  } else{
    mosaic <- suppressWarnings(terra::rast(mosaic, ...))
    if(terra::crs(mosaic) == ""){
      message("Missing Coordinate Reference System. Setting to EPSG:3857")
      terra::crs(mosaic) <- terra::crs("EPSG:3857")
    }
    if(terra::is.lonlat(mosaic)){
      eps <- mosaic_epsg(mosaic)
      warning(paste0("The current raster is in the lat/lon coordinate system, which may result in processing errors when trying to segment individuals in the `mosaic_analyze()` function. It is highly suggested to reproject the raster using mosaic_project() with ", eps), call. = FALSE)
    }
    if(check_16bits | check_datatype){
      cels <- sample(1:terra::ncell(mosaic), 2000, replace = TRUE)
      a <- na.omit(unlist(terra::extract(mosaic, cels)))
      a <- a[!is.infinite(a)]
      if(inherits(a, "numeric")){
        if(check_datatype){
          if(length(a[a - floor(a) != 0]) == 0){
            minv <- min(a)
            maxv <- max(a)
            if(all(minv >= 0) & all(maxv <= 255)){
              datatype <- "INT1U"
            } else{
              datatype <- "INT2U"
            }
          } else{
            datatype <- "FLT4S"
          }
          dtterra <- terra::datatype(mosaic)[[1]]
          if(datatype != dtterra){
            warning(paste("Based on the mosaic values, the datatype should be", datatype, "but it is ", dtterra,
                          ". Consider exporting it with `mosaic_export()` to assign the suggested datatype, which can save file size and memory usage during index computation. "))
          }
        }
        if(check_16bits){
          if(max(suppressWarnings(terra::minmax(mosaic)), na.rm = TRUE) == 65535){
            mosaic[mosaic == 65535] <- NA
          }
        }
      }
    }
    if(info){
      print(mosaic)
    }
    return(mosaic)
  }
}

#' @export
#' @name mosaic_input
mosaic_export <- function(mosaic,
                          filename,
                          datatype = NULL,
                          overwrite = FALSE,
                          ...){
  cels <- sample(1:terra::ncell(mosaic), 2000)
  a <- na.omit(unlist(terra::extract(mosaic, cels)))
  a <- a[!is.infinite(a)]

  if(is.null(datatype)){
    if(length(a[a - floor(a) != 0]) == 0){
      minv <- min(a)
      maxv <- max(a)
      if(all(minv >= 0) & all(maxv <= 255)){
        datatype <- "INT1U"
      } else{
        datatype <- "INT2U"
      }
    } else{
      datatype <- "FLT4S"
    }
  }
  message(paste0("Exporting the mosaic using datatype = ", datatype))
  terra::writeRaster(mosaic,
                     filename = filename,
                     overwrite = overwrite,
                     datatype = datatype,
                     gdal=c("COMPRESS=DEFLATE", "BIGTIFF=IF_NEEDED"),
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


#' SpatRaster aggregation
#'
#' Aggregate a SpatRaster to create a new SpatRaster with a lower resolution
#' (larger cells), using the GDAL's gdal_translate utility
#' https://gdal.org/programs/gdal_translate.html
#'
#' @param mosaic SpatRaster
#' @param pct The size as a fraction (percentage) of the input image size.
#'   Either a scalar (eg., 50), or a length-two numeric vector. In the last,
#'   different percentage reduction/expansion can be used for columns, and rows,
#'   respectively.
#' @param fun The resampling function. Defaults to `nearest`, which applies the
#'   nearest neighbor (simple sampling) resampler. Other accepted values are:
#'   'average', 'rms', 'bilinear', 'cubic', 'cubicspline', 'lanczos', and
#'   'mode'. See Details for a detailed explanation.
#' @param in_memory Wheter to return an 'in-memory' `SpatRaster`. If `FALSE`,
#'   the aggregated raster will be returned as an 'in-disk' object.
#' @return SpatRaster
#' @export
#'
#' @examples
#' library(pliman)
#' library(terra)
#' r <- rast()
#' values(r) <- 1:ncell(r)
#' r2 <- mosaic_aggregate(r, pct = 10)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2))
#' mosaic_plot(r)
#' mosaic_plot(r2)
#' par(opar)
mosaic_aggregate <- function(mosaic,
                             pct = 50,
                             fun = "nearest",
                             in_memory = TRUE){
  outsize <- compute_outsize(pct)
  td <- tempdir()
  if(terra::inMemory(mosaic)[[1]]){
    in_raster <- file.path(td, "tmp_aggregate.tif")
    terra::writeRaster(mosaic, in_raster, overwrite = TRUE)
    on.exit({
      file.remove(in_raster)
      if(in_memory){
        file.remove(out_raster)
      }
    })
  } else{
    in_raster <- terra::sources(mosaic)[[1]]
    on.exit(
      if(in_memory){
        file.remove(out_raster)
      }
    )
  }
  out_raster <- file.path(td, "tmp_aggregate_small.tif")
  sf::gdal_utils(
    util = "translate",
    source = in_raster,
    destination = out_raster,
    options = strsplit(paste("-r", fun, "-outsize", outsize[1], outsize[2]), split = "\\s")[[1]]
  )
  if(in_memory){
    terra::rast(out_raster) |> terra::wrap() |> terra::unwrap()
  } else{
    terra::rast(out_raster)
  }
}


#' A wrapper around terra::plot()
#'
#' Plot the values of a SpatRaster
#'
#' @param mosaic SpatRaster
#' @param col character vector to specify the colors to use. Defaults to
#'   `custom_palette(c("red", "yellow", "forestgreen"))`.
#' @param ... Further arguments passed on to [terra::plot()].
#' @param smooth logical. If TRUE (default) the cell values are smoothed (only
#'   if a continuous legend is used).
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' library(pliman)
#' r <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_plot(r)
mosaic_plot <- function(mosaic,
                        col = custom_palette(c("red", "yellow", "forestgreen"), n = 200),
                        smooth = TRUE,
                        ...){
  if(!inherits(mosaic, "SpatRaster")){
    stop("'mosaic' must be an object of class 'SpatRaster'")
  }
  terra::plot(mosaic,
              col = col,
              smooth = smooth,
              ...)
}

#' A wrapper around terra::hist()
#'
#' Create a histogram of the values of a `SpatRaster`.
#'
#' @param mosaic SpatRaster
#' @param layer positive integer or character to indicate layer numbers (or
#'   names). If missing, all layers are used
#' @param ... Further arguments passed on to [terra::hist()].
#'
#' @return A `NULL` object
#' @export
#'
#' @examples
#' library(pliman)
#' r <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' mosaic_hist(r)
mosaic_hist <- function(mosaic, layer, ...){
  if(!inherits(mosaic, "SpatRaster")){
    stop("'mosaic' must be an object of class 'SpatRaster'")
  }
  terra::hist(mosaic, layer, ...)
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
#' @inheritParams mosaic_index
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
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
#'   G, B, RE, NIR, SWIR, and TIR` bands (eg., `index = "R+B/G"`) or using the
#'   names of the mosaic's layers (ex., "(band_1 + band_2) / 2").
#' @param r,g,b,re,nir,swir,tir The red, green, blue, red-edge,  near-infrared,
#'   shortwave Infrared, and thermal infrared bands of the image, respectively.
#'   By default, the function assumes a BGR as input (b = 1, g = 2, r = 3). If a
#'   multispectral image is provided up to seven bands can be used to compute
#'   built-in indexes. There are no limitation of band numbers if the index is
#'   computed using the band name.
#' @param plot Plot the computed index? Defaults to `TRUE`.
#' @param in_memory Logical, indicating whether the indexes should be computed
#'   in memory. Defaults to `TRUE`. In most cases, this is 2-3 times faster, but
#'   errors can occur if `mosaic` is a large `SpatRaster`. If `FALSE`, raster
#'   algebra operations are performed on temporary files.
#' @param workers numeric. The number of workers you want to use for parallel
#'   processing when computing multiple indexes.
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
                         re = NA,
                         nir = NA,
                         swir = NA,
                         tir = NA,
                         plot = TRUE,
                         in_memory = TRUE,
                         workers = 1){
  indices <- c(r = r, g = g, b = b, re = re, nir = nir, swir = swir, tir = tir)
  valid_indices <- indices[!is.na(indices)]
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
    dimmo <- prod(dim(mosaic)[1:2])
    pattern <- "\\b\\w+\\b"
    reserved <- c("exp", "abs", "min", "max", "median", "sum", "sqrt", "cos", "sin", "tan", "log", "log10")
    layersused <- setdiff(unlist(regmatches(index, gregexpr(pattern, index, perl = TRUE))), reserved)
    onlychar <- suppressWarnings(is.na(as.numeric(layersused)))
    layers_used <- layersused[onlychar]
    if(!any(index  %in% ind$Index) & !all(layers_used  %in% c("R", "G", "B", "RE", "NIR", "SWIR", "TIR"))){
      # Extract individual layers based on the expression
      layers_used <- layers_used[is.na(suppressWarnings(as.numeric(layers_used)))]
      layers <-
        lapply(layers_used, function(x){
          mosaic[[x]]
        })
      names(layers) <- layers_used
      mosaic_gray <- eval(parse(text = index), envir = layers)
    } else{
      if(in_memory){
        if(index %in% ind$Index){
          formula <- as.character(ind$Equation[as.character(ind$Index)==index])
          mosaic_gray <- terra::lapp(mosaic[[valid_indices]], parse_formula(formula, valid_indices))
        } else{
          mosaic_gray <- terra::lapp(mosaic[[valid_indices]], parse_formula(index, valid_indices))
        }
      } else{
        R <- try(ras[[indices[["r"]]]], TRUE)
        G <- try(ras[[indices[["g"]]]], TRUE)
        B <- try(ras[[indices[["b"]]]], TRUE)
        RE <- try(ras[[indices[["re"]]]], TRUE)
        NIR <- try(ras[[indices[["nir"]]]], TRUE)
        SWIR <- try(ras[[indices[["swir"]]]], TRUE)
        TIR <- try(ras[[indices[["tir"]]]], TRUE)
        if(index %in% ind$Index){
          mosaic_gray <- eval(parse(text = as.character(ind$Equation[as.character(ind$Index)==index])))
        } else{
          mosaic_gray <- eval(parse(text = index))
        }
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
    if(workers > 1){
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential))
      `%dofut%` <- doFuture::`%dofuture%`
      if(terra::inMemory(mosaic)){
        tf <- paste0(tempfile(), ".tif")
        on.exit(file.remove(tf))
        terra::writeRaster(mosaic, filename = tf)
        tempf <- tf
      } else{
        tempf <- terra::sources(mosaic)
      }
      d <-
        foreach::foreach(i = seq_along(index),
                         .options.future = list(
                           seed = TRUE
                         )) %dofut%{
                           tfil <- paste0(tempfile(), ".tif")
                           mosaic_index(terra::rast(tempf),
                                        index = unique(index)[[i]],
                                        r = r,
                                        g = g,
                                        b = b,
                                        re = re,
                                        nir = nir,
                                        swir = swir,
                                        tir = tir,
                                        in_memory = in_memory,
                                        plot = FALSE) |>
                             terra::writeRaster(tfil, overwrite = TRUE)
                           tfil
                         }
      mosaic_gray <- terra::rast(unlist(d)) |> terra::wrap() |> terra::unwrap()
      on.exit(file.remove(unlist(d)))
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
                           swir = swir,
                           tir = tir,
                           in_memory = in_memory,
                           plot = FALSE)
            })
        )
      )
    }
  }
  if(plot){
    mosaic_plot(mosaic_gray)
  }
  invisible(mosaic_gray)
}

#' Mosaic Index with GDAL
#'
#' Compute or extract an index layer from a multi-band mosaic raster using
#' gdal_calc.py (https://gdal.org/programs/gdal_calc.html). This requires a
#' Python and GDAL installation.
#' @inheritParams mosaic_index
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
#' @param python The PATH for python.exe
#' @param gdal The PATH for gdal_calc.py
#' @return An index layer extracted/computed from the mosaic raster.
#' @export
#' @examples
#' if((Sys.which('python.exe') != '' ) & (Sys.which('gdal_calc.py') != '' )){
#' library(pliman)
#' mosaic <- mosaic_input(system.file("ex/elev.tif", package="terra"))
#' names(mosaic) <- "R"
#' elev2 <- mosaic_index2(mosaic, "R * 5", plot = FALSE)
#' oldpar <- par(no.readonly=TRUE)
#' mosaic_plot(mosaic)
#' mosaic_plot(elev2)
#' par(mfrow=c(1,2))
#' }

mosaic_index2 <- function(mosaic,
                          index = "B",
                          r = 3,
                          g = 2,
                          b = 1,
                          re = 4,
                          nir = 5,
                          plot = TRUE,
                          python = Sys.which('python.exe'),
                          gdal = Sys.which('gdal_calc.py')) {
  if(python == ''){
    stop('Error: Python executable (python.exe) not found on the system. Please, install it and add it to the PATH variable')
  }
  if(gdal==''){
    stop("Error: gdal_calc.py not found on the system. Make sure GDAL is installed and available in the system PATH.
         You can install GDAL by installing miniconda (https://docs.conda.io/projects/miniconda/en/latest/)
         and then running 'conda install -c conda-forge gdal' in the conda terminal.
         Additionally, make sure that the Python Scripts directory is added to the system PATH.")
  }
  ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
  checkind <- index %in% ind$Index
  if (!checkind) {
    message(paste("Index '", paste0(index[!checkind], collapse = ", "), "' is not available. Trying to compute your own index.",
                  sep = ""))
  } else{
    index <- as.character(ind$Equation[as.character(ind$Index)==index])
  }
  if(terra::inMemory(mosaic)){
    tf <- tempfile(fileext = ".tif")
    mosaic_export(mosaic, tf)
    on.exit(file.remove(tf))
    infile <- tf
  } else{
    infile <- terra::sources(mosaic)
  }
  mosaicbands <- c("R", "G", "B", "E", "I")
  outfile <- tempfile(fileext = ".tif")
  nbands <- terra::nlyr(mosaic)
  inputs <- paste0('-',
                   mosaicbands[seq_len(nbands)], ' ', infile, ' --',
                   mosaicbands[seq_len(nbands)], '_band ', seq_len(nbands), collapse=' ')
  system2(python,
          args=c(gdal,
                 inputs,
                 sprintf("--outfile=%s", outfile),
                 sprintf('--calc="%s"', index),
                 '--co="COMPRESS=DEFLATE"',
                 '--co="BIGTIFF=IF_NEEDED"',
                 '--overwrite'),
          stdout=FALSE
  )
  res <- terra::rast(outfile)
  if(plot){
    terra::plot(res)
  }
  return(res)
}

#' Segment a mosaic
#'
#' Segment a `SpatRaster` using a computed image index. By default, values
#' greater than `threshold` are kept in the mask.
#'
#' @inheritParams mosaic_index
#' @inheritParams mosaic_analyze
#' @param return The output of the function. Either 'mosaic' (the segmented
#'   mosaic), or 'mask' (the binary mask).
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
                           re = NA,
                           nir = NA,
                           swir = NA,
                           tir = NA,
                           threshold = "Otsu",
                           invert = FALSE,
                           return = c("mosaic", "mask")){
  if(!return[[1]] %in% c("mosaic", "mask")){
    stop("'return' must be one of 'mosaic' or 'mask'.")
  }
  ind <- mosaic_index(mosaic,
                      index = index,
                      r = r,
                      g = g,
                      b = b,
                      re = re,
                      nir = nir,
                      swir = swir,
                      tir = tir,
                      plot = FALSE)
  thresh <- ifelse(threshold == "Otsu", otsu(na.omit(terra::values(ind)[, index])), threshold)
  if(invert){
    mask <- ind[[index]] > thresh
  } else{
    mask <- ind[[index]] < thresh
  }
  if(return[[1]] == 'mosaic'){
    terra::mask(mosaic, mask, maskvalue = TRUE)
  } else{
    mask
  }
}


#' Segments a mosaic interactively
#'
#' The function segments a mosaic using an interative process where the user
#' picks samples from background (eg., soil) and foreground (eg., plants).
#'
#' @inheritParams mosaic_index
#' @inheritParams mosaic_view
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
                                max_pixels = 2e6,
                                downsample = NULL,
                                quantiles = c(0, 1),
                                return = c("mosaic", "mask")){
  if(!return[[1]] %in% c("mosaic", "mask")){
    stop("'return' must be one of 'mosaic' or 'mask'.")
  }
  downsample <- ifelse(is.null(downsample), find_aggrfact(mosaic, max_pixels = max_pixels), downsample)
  if(downsample > 0){
    mosaic <- mosaic_aggregate(mosaic, pct = round(100 / downsample))
  }
  if(is.null(basemap)){
    basemap <-
      mosaic_view(mosaic,
                  r = r,
                  g = g,
                  b = b,
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
    dplyr::bind_rows() |>
    dplyr::select(-coverage_fraction) |>
    dplyr::mutate(class = 0)

  mapview::mapview() |> mapedit::editMap()

  plant <- mapedit::editMap(basemap,
                            title = "Use the 'Draw Rectangle' tool to pick up foreground fractions",
                            editor = "leafpm")$finished
  plant <- plant |> sf::st_transform(sf::st_crs(mosaic))
  plant_sample <-
    exactextractr::exact_extract(mosaic, plant, progress = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::select(-coverage_fraction) |>
    dplyr::mutate(class = 1)
  df_train <- dplyr::bind_rows(plant_sample, soil_sample)
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
    terra::mask(mosaic, mask, maskvalue = TRUE)
  } else{
    mask
  }
}


#' Mosaic to pliman
#'
#' Convert an `SpatRaster` object to a `Image` object with optional scaling.
#' @inheritParams mosaic_view
#' @inheritParams mosaic_index
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
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
#' @param r,g,b The red, green, blue bands.
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
#' @inheritParams mosaic_index
#' @inheritParams mosaic_to_pliman
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
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
#' @inheritParams mosaic_index
#' @inheritParams analyze_objects
#' @param r,g,b,re,nir The red, green, blue, red-edge, and  near-infrared bands
#'   of the image, respectively. By default, the function assumes a BGR as input
#'   (b = 1, g = 2, r = 3). If a multispectral image is provided up to seven
#'   bands can be used to compute built-in indexes. There are no limitation of
#'   band numbers if the index is computed using the band name.
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
      terra::extractAlong(x = mind,
                          y = polygons_ext,
                          ID = FALSE)
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
      if(invert){
        mask <- mind[[1]] < otsu(na.omit(terra::values(mind)[, index[1]]))
      } else{
        mask <- mind[[1]] < otsu(na.omit(terra::values(mind)[, index[1]]))
      }
      mind <- terra::mask(mind, mask, maskvalues = TRUE)
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


#' Convert Sentinel data to GeoTIFF format
#'
#' This function converts Sentinel satellite data files to GeoTIFF format.
#'
#' @param layers (character) Vector of file paths to Sentinel data files. If
#'   NULL, the function searches for files in the specified path with names
#'   containing "B".
#' @param path (character) Directory path where Sentinel data files are located.
#'   Default is the current directory.
#' @param destination (character) File path for the output GeoTIFF file.
#' @param spat_res (numeric) Spatial resolution of the output GeoTIFF file.
#'   Default is 10 meters.
#'
#' @details The function converts Sentinel satellite data files to GeoTIFF
#'   format using GDAL utilities. It builds a virtual raster file (VRT) from the
#'   input files and then translates it to GeoTIFF format. Compression is
#'   applied to the output GeoTIFF file using DEFLATE method.
#'
#'
#'
#' @export
sentinel_to_tif <- function(layers = NULL,
                            path = ".",
                            destination,
                            spat_res = 10){

  if(is.null(layers)){
    files <- list.files(path = path, pattern = "B")
  } else{
    files <- layers
  }
  tf <- tempfile(fileext = ".vrt")
  sf::gdal_utils(
    util = "buildvrt",
    source  = files,
    destination = tf,
    options = strsplit(paste("-separate -tr", spat_res, spat_res), split = "\\s")[[1]]
  )

  sf::gdal_utils(
    util = "translate",
    source  = tf,
    destination = paste0(path, "/", destination),
    options = strsplit(paste("-co COMPRESS=DEFLATE", "-of GTiff"), split = "\\s")[[1]]
  )
}

#' Calculate Canopy Height Model and Volume
#'
#' This function calculates the canopy height model (CHM) and the volume for a
#' given digital surface model (DSM) raster layer. Optionally, a digital terrain
#' model (DTM) can be provided or interpolated using a set of points or a moving
#' window.
#'
#' @param dsm A `SpatRaster` object representing the digital surface model. Must
#'   be a single-layer raster.
#' @param dtm (optional) A `SpatRaster` object representing the digital terrain
#'   model. Must be a single-layer raster. If not provided, it can be
#'   interpolated from points or created using a moving window.
#' @param points (optional) An `sf` object representing sample points for DTM
#'   interpolation. If provided, `dtm` will be interpolated using these points.
#' @param interpolation (optional) A character string specifying the
#'   interpolation method to use when `points` are provided. Options are
#'   "Kriging" (default) or "Tps" (Thin Plate Spline).
#' @param window_size An integer  (meters) specifying the window size (rows and
#'   columns, respectively) for creating a DTM using a moving window. Default is
#'   c(10, 10).
#' @param mask (optional) A `SpatRaster` object used to mask the CHM and volume
#'   results. Default is NULL.
#' @param mask_soil Is `mask` representing a soil mask (eg., removing plants)? Default is TRUE.
#'
#' @return A `SpatRaster` object with three layers: `dtm` (digital terrain
#'   model), `height` (canopy height model), and `volume`.
#'
#' @details
#' The function first checks if the input `dsm` is a valid single-layer
#' `SpatRaster` object. If `dtm` is not provided, The function generates a
#' Digital Terrain Model (DTM) from a Digital Surface Model (DSM) by
#' downsampling and smoothing the input raster data. It iterates over the DSM
#' matrix in windows of specified size, finds the minimum value within each
#' window, and assigns these values to a downsampled matrix. After downsampling,
#' the function applies a mean filter to smooth the matrix, enhancing the visual
#' and analytical quality of the DTM. Afterwards, DTM is resampled with the
#' original DSM.
#'
#' If both `dsm` and `dtm` are provided, the function ensures they have the same
#' extent and number of cells, resampling `dtm` if necessary. The CHM is then
#' calculated as the difference between `dsm` and `dtm`, and the volume is
#' calculated by multiplying the CHM by the pixel size. The results are
#' optionally masked using the provided `mask`.
#'
#' @importFrom fields Krig Tps
#' @export

mosaic_chm <- function(dsm,
                       dtm = NULL,
                       points = NULL,
                       interpolation = c("Tps", "Kriging"),
                       window_size = c(10, 10),
                       mask = NULL,
                       mask_soil = TRUE,
                       verbose = TRUE){
  sampp <- NULL
  ch1 <- !inherits(dsm,"SpatRaster") || !terra::nlyr(dsm) == 1 || terra::is.bool(dsm) || is.list(dsm)
  if(ch1){
    stop("dsm must be single-layer SpatRaster objects")
  }
  # interpolate dtm using sample of points
  if(is.null(dtm) & !is.null(points)){
    # sampling points
    points <- points |> sf::st_transform(sf::st_crs(dsm))
    if(verbose){
      message("\014","\nExtracting values...\n")
    }
    vals <- terra::extract(dsm, terra::vect(points), xy = TRUE)
    xy <- cbind(vals$x,vals$y)
    z <- vals[, 2]
    if(verbose){
      message("\014","\nInterpolating the raster...\n")
    }
    if(interpolation[[1]] == "Kriging"){
      fit <- suppressMessages(suppressWarnings(fields::Krig(xy, z, aRange=20)))
    }
    if(interpolation[[1]] == "Tps"){
      fit <- suppressMessages(suppressWarnings(fields::Tps(xy, z)))
    }
    sampp <- NULL
    # low resolution to interpolate
    aggr <- find_aggrfact(dsm, 4e5)
    if(aggr > 0){
      mosaicintp <- mosaic_aggregate(dsm, round(100 / aggr))
    } else{
      mosaicintp <- dsm
    }
    dtm <- terra::interpolate(terra::rast(mosaicintp), fit)
    gc()
    terra::crs(dtm) <- terra::crs(dsm)
    if(verbose){
      message("\014","\nResampling and masking the interpolated raster...\n")
    }
    dtm <- terra::resample(dtm, dsm)
    gc()
    dtm <- terra::mask(dtm, dsm)
  }
  # create a dtm using a moving window that extract the minimum values
  # from dsm
  if(is.null(dtm) & is.null(points)){
    resolu <- terra::res(dsm)
    extens <- terra::ext(dsm)
    wide <- extens[2] - extens[1]
    heig <- extens[4] - extens[3]
    nr <- ceiling( heig / window_size[[1]])
    nc <- ceiling( wide / window_size[[2]])
    if(verbose){
      message("\014","\nExtracting minimum value for each moving window...\n")
    }
    shp <- shapefile_build(dsm,
                           nrow = nr,
                           ncol = nc,
                           build_shapefile = FALSE,
                           verbose = FALSE)
    vals <- exactextractr::exact_extract(dsm,
                                         shp[[1]],
                                         fun = "min",
                                         progress = FALSE)
    gc()
    cent <- suppressWarnings(sf::st_centroid(shp[[1]]))
    sampp <-
      cent |>
      dplyr::mutate(dtm = vals) |>
      dplyr::filter(!is.na(dtm))

    xy <- sf::st_coordinates(sampp)
    z <- sampp$dtm

    if(verbose){
      message("\014","\nInterpolating the raster...\n")
    }
    if(interpolation[[1]] == "Kriging"){
      fit <- suppressMessages(suppressWarnings(fields::Krig(xy, z, aRange=20)))
    }
    if(interpolation[[1]] == "Tps"){
      fit <- suppressMessages(suppressWarnings(fields::Tps(xy, z)))
    }

    # low resolution to interpolate
    aggr <- find_aggrfact(dsm, 4e5)
    if(aggr > 0){
      mosaicintp <- mosaic_aggregate(dsm, round(100 / aggr))
    } else{
      mosaicintp <- dsm
    }
    dtm <- terra::interpolate(terra::rast(mosaicintp), fit)
    terra::crs(dtm) <- terra::crs(dsm)
    if(verbose){
      message("\014","\nResampling and masking the interpolated raster...\n")
    }
    dtm <- terra::resample(dtm, dsm)
    gc()
    dtm <- terra::mask(dtm, dsm)
    gc()
  }
  # now, create a chm
  if(!is.null(dtm)){
    ch2 <- !inherits(dtm,"SpatRaster") || !terra::nlyr(dtm) == 1 || terra::is.bool(dtm) || is.list(dtm)
    if(ch2){
      stop("dtm must be single-layer SpatRaster objects")
    }
    if((terra::ext(dsm) != terra::ext(dtm)) || (terra::ncell(dtm) != terra::ncell(dsm))){
      dtm <- terra::resample(dtm, dsm)
    }
    if(verbose){
      message("\014","\nBuilding the digital terrain model...\n")
    }
    chm <- dsm - dtm
    gc()
    psize <- prod(terra::res(chm))
    volume <- chm * psize
    chm <- c(chm, volume)
    gc()
    if(!is.null(mask)){
      if((terra::ext(mask) != terra::ext(dsm)) || (terra::ncell(mask) != terra::ncell(dsm))){
        mask <- terra::resample(mask, dsm)
      }
      chm <- terra::mask(chm,  mask, maskvalues = mask_soil)
    }
    chm <- c(dtm, chm)
    names(chm) <- c("dtm", "height", "volume")
  }
  if(verbose){
    message("\014","\nDone!\n")
  }
  return(list(chm = chm, sampling_points = sampp, mask = ifelse(is.null(mask), FALSE, TRUE)))
}

#' Extract Canopy Height and Volume
#'
#' This function extracts canopy height and volume metrics for given plots
#' within a specified shapefile.
#' @param chm A list object containing the Canopy Height Model (CHM) generated
#'   by the [mosaic_chm()] function.
#' @param shapefile An `sf` object representing the plot boundaries for which
#'   the metrics will be extracted.
#'
#' @return A `sf` object with extracted metrics including minimum, 10th
#'   percentile, median (50th percentile), 90th percentile, interquartile range
#'   (IQR), mean, maximum canopy height, coefficient of variation (CV) of canopy
#'   height, canopy height entropy, total volume, covered area, plot area, and
#'   coverage percentage. Centroid coordinates (x, y) of each plot are also
#'   included.
#' @details
#' The function uses the `exactextractr` package to extract canopy height and
#' volume metrics from the CHM. For each plot in the shapefile, the function
#' computes various statistics on the canopy height values (e.g., min, max,
#' percentiles, mean, CV, entropy) and sums the volume values. If a mask was
#' applied in the CHM calculation, the covered area and plot area are also
#' computed.
#' @export
mosaic_chm_extract <- function(chm, shapefile){
  custom_summary <- function(values, coverage_fractions, ...) {
    valids <- na.omit(values)
    entropy <- function(values) {
      freq <- table(round(values, 2))
      prob <- freq / sum(freq)
      entropy <- -sum(prob * log(prob))
      return(entropy)
    }
    quantiles <- quantile(valids, c(0.1, 0.5, 0.9))
    data.frame(
      min = min(valids),
      q10 = quantiles[[1]],
      q50 = quantiles[[2]],
      q90 = quantiles[[3]],
      iqr = IQR(valids),
      mean = sum(valids) / length(valids),
      max = max(valids),
      cv = sd(valids) / mean(valids),
      entropy = entropy(valids)
    )
  }
  height <- exactextractr::exact_extract(chm$chm[[2]],
                                         shapefile,
                                         fun = custom_summary,
                                         force_df = TRUE,
                                         progress = FALSE)
  vol <- exactextractr::exact_extract(chm$chm[[3]],
                                      shapefile,
                                      fun = "sum",
                                      force_df = TRUE,
                                      progress = FALSE)
  names(vol) <- c("volume")

  # include check here if mask is not present
  if(chm$mask){
    area <- exactextractr::exact_extract(chm$chm[[3]],
                                         shapefile,
                                         coverage_area = TRUE,
                                         force_df = TRUE,
                                         progress = FALSE)
    covered_area <-
      purrr::map_dfr(area, function(x){
        data.frame(covered_area = sum(na.omit(x)[, 2]),
                   plot_area = sum(x[, 2]))
      }) |>
      dplyr::mutate(coverage = covered_area / plot_area)
  } else{
    area <- as.numeric(sf::st_area(shapefile))
    covered_area <- data.frame(covered_area = area,
                               plot_area = area,
                               coverage = 1)
  }
  centroids <- suppressWarnings(sf::st_centroid(shapefile)) |> sf::st_coordinates()
  colnames(centroids) <- c("x", "y")
  dftmp <-
    dplyr::bind_cols(height, vol, covered_area, centroids, shapefile) |>
    sf::st_as_sf() |>
    dplyr::relocate(unique_id, block, plot_id, row, column, x, y, .before = 1)
  return(dftmp)
}
#' Determine EPSG Code for a Mosaic
#'
#' This function calculates the EPSG code for a given mosaic based on its
#' geographic extent.
#'
#' @param mosaic A raster object representing the mosaic for which the EPSG code
#'   is to be determined.
#'
#' @return A character string representing the EPSG code corresponding to the
#'   UTM zone and hemisphere of the mosaic's centroid. If the mosaic is not in
#'   the lon/lat coordinate system, a warning is issued.
#'
#' @details The function calculates the centroid of the mosaic's extent,
#'   determines the UTM zone based on the centroid's longitude, and identifies
#'   the hemisphere based on the centroid's latitude. The EPSG code is then
#'   constructed accordingly.
#'
#' @examples
#' \dontrun{
#' library(pliman)
#' library(terra)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#'
#' # Get the EPSG code for the mosaic
#' mosaic_epsg(mosaic)
#' }
#'
#' @export
mosaic_epsg <- function(mosaic) {
  if(terra::is.lonlat(mosaic)){
    extens <- terra::ext(mosaic)
    latitude <- mean(c(extens[3], extens[4]))
    longitude <- mean(c(extens[1], extens[2]))
    utm_zone <- floor((longitude + 180) / 6) + 1
    hemisphere <- ifelse(latitude >= 0, "N", "S")
    epsg_code <- if (hemisphere == "N") {
      32600 + utm_zone
    } else {
      32700 + utm_zone
    }
    return(paste0("EPSG:", epsg_code))
  } else{
    warning("`mosaic` is not in the lon/lat coordinate system.")
  }
}

#' Project a Mosaic to a New Coordinate Reference System (CRS)
#'
#' This function projects a given mosaic to a specified CRS.
#'
#' @param mosaic A raster object representing the mosaic to be projected.
#' @param y The target CRS to which the mosaic should be projected. This can be
#'   specified in various formats accepted by the [terra::project()] function.
#' @param ... Additional arguments passed to the [terra::project()] function.
#'
#' @return A raster object representing the projected mosaic.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(pliman)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#' mosaic
#' # Define target CRS (EPSG code for WGS 84 / UTM zone 33N)
#' target_crs <- "EPSG:32633"
#'
#' # Project the mosaic
#' projected_mosaic <- mosaic_project(mosaic, "EPSG:32633")
#' projected_mosaic
#' }
#'
#' @export
mosaic_project <- function(mosaic, y, ...){
  return(terra::project(mosaic, y, ...))
}

#' Project a Mosaic from Lon/Lat to EPSG-based CRS
#'
#' This function projects a given mosaic from the lon/lat coordinate system to
#' an EPSG-based CRS determined by the mosaic's extent.
#'
#' @param mosaic A raster object representing the mosaic to be projected. The
#'   mosaic must be in the lon/lat coordinate system.
#'
#' @return A raster object representing the projected mosaic. If the mosaic is
#'   not in the lon/lat coordinate system, a warning is issued.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(pliman)
#'
#' # Create a sample mosaic
#' mosaic <- rast(nrow=10, ncol=10, xmin=-120, xmax=-60, ymin=30, ymax=60)
#'
#' # Project the mosaic to the appropriate UTM zone
#' mosaic_lonlat2epsg(mosaic)
#' }
#'
#' @export
mosaic_lonlat2epsg <- function(mosaic){
  if(terra::is.lonlat(mosaic)){
    epsg <- mosaic_epsg(mosaic)
    return(terra::project(mosaic, epsg))
  } else{
    warning("`mosaic` is not in the lon/lat coordinate system.")
  }
}

#' Extract Values from a Raster Mosaic Using a Shapefile
#'
#' This function extracts values from a raster mosaic based on the regions
#' defined in a shapefile using [exactextractr::exact_extract()].
#'
#' @param mosaic A `SpatRaster` object representing the raster mosaic from which
#'   values will be extracted.
#' @param shapefile A shapefile, which can be a `SpatVector` or an `sf` object,
#'   defining the regions of interest for extraction.
#' @param fun A character string specifying the summary function to be used for
#'   extraction. Default is `"median"`.
#' @param ... Additional arguments to be passed to [exactextractr::exact_extract()].
#' @return A data frame containing the extracted values for each region defined in the shapefile.
#' @export
#'
mosaic_extract <- function(mosaic,
                           shapefile,
                           fun = "median",
                           ...){
  if(inherits(shapefile, "SpatVector")){
    shapefile <- sf::st_as_sf(shapefile)
  }
  exactextractr::exact_extract(mosaic,
                               shapefile,
                               fun = fun,
                               force_df = TRUE,
                               ...)
}
