#' Analyzes objects in an image
#'
#' * [analyze_objects()] provides tools for counting and extracting object
#' features (e.g., area, perimeter, radius, pixel intensity) in an image. See
#' more at **Details** section.
#' * [plot.anal_obj()] Produces an histogram for the R, G, and B values when
#' argument `object_index` is used in the function [analyze_objects()].
#'
#' @details
#'A binary image is first generated to segment the foreground and background.
#'The argument `index` is useful to choose a proper index to segment the image
#'(see [image_binary()] for more details). Then, the number of objects in the
#'foreground is counted. By setting up arguments such as `lower_size`,
#'`upper_size` it is possible to set a threshold for lower and upper sizes of
#'the objects, respectively. The argument `object_size` can be used to set up
#'pre-defined values of `tolerance` and `extension` depending on the image
#'resolution. This will influence the watershed-based object segmentation. Users
#'can also tune-up `tolerance` and `extension` explicitly to a better precision
#'of watershed segmentation.
#'
#'If `watershed = FALSE` is used, all pixels for each connected set of
#'foreground pixels in `img` are set to a unique object. This is faster
#'(specially for a large number of objects) but is not able to segment touching
#'objects.
#'
#'If color palettes samples are provided, a general
#'linear model (binomial family) fitted to the RGB values is used to segment
#'fore- and background.
#'
#'By using `pattern` it is possible to process several images with common
#'pattern names that are stored in the current working directory or in the
#'subdirectory informed in `dir_original`'. To speed up the computation time,
#'one can set `parallel = TRUE`.
#' @param img The image to be analyzed.
#' @param foreground A color palette of the foreground (optional).
#' @param background A color palette of the background (optional).
#' @param pattern A pattern of file name used to identify images to be imported.
#'   For example, if `pattern = "im"` all images in the current working
#'   directory that the name matches the pattern (e.g., img1.-, image1.-, im2.-)
#'   will be imported as a list. Providing any number as pattern (e.g., `pattern
#'   = "1"`) will select images that are named as 1.-, 2.-, and so on. An error
#'   will be returned if the pattern matches any file that is not supported
#'   (e.g., img1.pdf).
#' @param parallel If `TRUE` processes the images asynchronously (in parallel)
#'   in separate R sessions running in the background on the same machine. It
#'   may speed up the processing time, especially when `pattern` is used is
#'   informed. When `object_index` is informed, multiple sections will be used
#'   to extract the RGB values for each object in the image. This may
#'   significantly speed up processing time when an image has lots of objects
#'   (say >1000).
#' @param workers A positive numeric scalar or a function specifying the number
#'   of parallel processes that can be active at the same time. By default, the
#'   number of sections is set up to 50% of available cores.
#' @param watershed If `TRUE` (default) performs watershed-based object
#'   detection. This will detect objects even when they are touching one other.
#'   If `FALSE`, all pixels for each connected set of foreground pixels are set
#'   to a unique object. This is faster but is not able to segment touching
#'   objects.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param trim Number of pixels removed from edges in the analysis. The edges of
#'   images are often shaded, which can affect image analysis. The edges of
#'   images can be removed by specifying the number of pixels. Defaults to
#'   `FALSE` (no trimmed edges).
#' @param fill_hull Fill holes in the binary image? Defaults to `FALSE`. This is
#'   useful to fill holes in objects that have portions with a color similar to
#'   the background. IMPORTANT: Objects touching each other can be combined into
#'   one single object, which may underestimate the number of objects in an
#'   image.
#' @param filter Performs median filtering after image processing? defaults to
#'   `FALSE`. See more at [image_filter()].
#' @param invert Inverts the binary image, if desired. This is useful to process
#'   images with black background. Defaults to `FALSE`.
#' @param object_size The size of the object. Used to automatically set up
#'   `tolerance` and `extension` parameters. One of the following. `"small"`
#'   (e.g, wheat grains), `"medium"` (e.g, soybean grains), `"large"`(e.g,
#'   peanut grains), and `"elarge"` (e.g, soybean pods)`.
#' @param index,my_index A character value specifying the target mode for
#'   conversion to binary image when `foreground` and `background` are not
#'   declared. Defaults to `"NB"` (normalized blue). See [image_index()] for
#'   more details.
#' @param object_index Defaults to `FALSE`. If an index is informed, the average
#'   value for each object is returned. It can be the R, G, and B values or any
#'   operation involving them, e.g., `object_index = "R/B"`. In this case, it
#'   will return for each object in the image, the average value of the R/B
#'   ratio. Use [pliman_indexes_eq()] to see the equations of available indexes.
#' @param threshold By default (`threshold = "Otsu"`), a threshold value based
#'   on Otsu's method is used to reduce the grayscale image to a binary image.
#'   If a numeric value is informed, this value will be used as a threshold.
#'   Inform any non-numeric value different than "Otsu" to iteratively chosen
#'   the threshold based on a raster plot showing pixel intensity of the index.
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Higher value smooths out small objects.
#' @param lower_size,upper_size Lower and upper limits for size for the image
#'   analysis. Plant images often contain dirt and dust. To prevent dust from
#'   affecting the image analysis, objects with lesser than 10% of the mean of
#'   all objects are removed. Upper limit is set to `NULL`, i.e., no upper
#'   limit used. One can set a known area or use `lower_limit = 0` to select all
#'   objects (not advised). Objects that matches the size of a given range of
#'   sizes can be selected by setting up the two arguments. For example, if
#'   `lower_size = 120` and `upper_size = 140`, objects with size greater than
#'   or equal 120 and less than or equal 140 will be considered.
#' @param topn_lower,topn_upper Select the top `n` objects based on its area.
#'   `topn_lower` selects the `n` elements with the smallest area whereas
#'   `topn_upper` selects the `n` objects with the largest area.
#' @param lower_eccent,upper_eccent,lower_circ,upper_circ Lower and upper limit
#'   for object eccentricity/circularity for the image analysis. Users may use
#'   these arguments to remove objects such as square papers for scale (low
#'   eccentricity) or cut petioles (high eccentricity) from the images. Defaults
#'   to `NULL` (i.e., no lower and upper limits).
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step. Defaults to
#'   2000.
#' @param show_image Show image after processing?
#' @param show_original Show the count objects in the original image?
#' @param show_chull Show the convex hull around the objects? Defaults to
#'   `FALSE`.
#' @param show_contour Show a contour line around the objects? Defaults
#'   to `TRUE`.
#' @param contour_col,contour_size The color and size for the contour line
#'   around objects. Defaults to `contour_col = "red"` and `contour_size = 1`.
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param show_segmentation Shows the object segmentation colored with random
#'   permutations. Defaults to `FALSE`.
#' @param col_foreground,col_background Foreground and background color after
#'   image processing. Defaults to `NULL`, in which `"black"`, and `"white"` are
#'   used, respectively.
#' @param marker,marker_col,marker_size The type, color and size of the object
#'   marker. Defaults to `NULL`, which plots the object id. Use `marker =
#'   "point"` to show a point in each object or `marker = FALSE` to omit object
#'   marker.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named as `proc_*` where `*` is the image name
#'   given in `img`.
#' @param prefix The prefix to be included in the processed images. Defaults to
#'   `"proc_"`.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory. It can be either a full path, e.g.,
#'   `"C:/Desktop/imgs"`, or a subfolder within the current working directory,
#'   e.g., `"/imgs"`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return `analyze_objects()` returns a list with the following objects:
#'  * `results` A data frame with the following variables for each object in the
#'  image:
#'     - `id`:  object identification.
#'     - `x`,`y`:  x and y coordinates for the center of mass of the object.
#'     - `area`:  area of the object (in pixels).
#'     - `area_ch`:  the area of the convex hull around object (in pixels).
#'     - `perimeter`: perimeter (in pixels).
#'     - `radius_min`, `radius_mean`, and `radius_max`: The minimum, mean, and
#'     maximum radius (in pixels), respectively.
#'     - `radius_sd`: standard deviation of the mean radius (in pixels).
#'     - `radius_ratio`: radius ratio given by `radius_max / radius_min`.
#'     - `diam_min`, `diam_mean`, and `diam_max`: The minimum, mean, and
#'     maximum diameter (in pixels), respectively.
#'     - `major_axis`, `minor_axis`: elliptical fit for major and minor axes (in
#'     pixels).
#'     - `eccentricity`: elliptical eccentricity defined by
#'     sqrt(1-minoraxis^2/majoraxis^2). Circle eccentricity is 0 and straight
#'     line eccentricity is 1.
#'     - `theta`: object angle (in radians).
#'     - `solidity`: object solidity given by `area / area_ch`.
#'     - `circularity`: the object circularity given by \eqn{4*pi *(area /
#'     perimeter^2)}.
#'  * `statistics`: A data frame with the summary statistics for the area of the
#'  objects.
#'  * `count`: If `pattern` is used, shows the number of objects in each image.
#'  * `object_rgb`: If `object_index` is used, returns the R, G, and B values
#'  for each pixel of each object.
#'  * `object_index`: If `object_index` is used, returns the index computed for
#'  each object.
#' @references
#' Gupta, S., Rosenthal, D. M., Stinchcombe, J. R., & Baucom, R. S. (2020). The
#' remarkable morphological diversity of leaf shape in sweet potato (Ipomoea
#' batatas): the influence of genetics, environment, and G×E. New Phytologist,
#' 225(5), 2183–2195. \doi{10.1111/NPH.16286}
#'
#' Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area of a Polygon
#' and the Vector Cross Product. The Mathematics Teacher, 110(8), 631–636.
#' \doi{10.5951/mathteacher.110.8.0631}
#'
#' @export
#' @name analyze_objects
#' @importFrom  utils install.packages
#' @importFrom grDevices col2rgb dev.off jpeg png
#' @importFrom graphics lines par points rect text
#' @importFrom stats aggregate binomial glm kmeans predict sd
#' @importFrom utils menu
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' obj <- analyze_objects(img)
#' obj$statistics
#'
#' # Enumerate the objects in the original image
#' # Return the top-5 grains with the largest area
#'
#' top <-
#'  analyze_objects(img,
#'                  marker = "id",
#'                  topn_upper = 5)
#' top$results
#' }
#'
analyze_objects <- function(img,
                            foreground = NULL,
                            background = NULL,
                            pattern = NULL,
                            parallel = FALSE,
                            workers = NULL,
                            watershed = TRUE,
                            resize = FALSE,
                            trim = FALSE,
                            fill_hull = FALSE,
                            filter = FALSE,
                            invert = FALSE,
                            object_size = "medium",
                            index = "NB",
                            my_index = NULL,
                            object_index = NULL,
                            threshold = "Otsu",
                            tolerance = NULL,
                            extension = NULL,
                            lower_size = NULL,
                            upper_size = NULL,
                            topn_lower = NULL,
                            topn_upper = NULL,
                            lower_eccent = NULL,
                            upper_eccent = NULL,
                            lower_circ = NULL,
                            upper_circ = NULL,
                            randomize = TRUE,
                            nrows = 2000,
                            show_image = TRUE,
                            show_original = TRUE,
                            show_chull = FALSE,
                            show_contour = TRUE,
                            contour_col = "red",
                            contour_size = 1,
                            show_background = TRUE,
                            show_segmentation = FALSE,
                            col_foreground = NULL,
                            col_background = NULL,
                            marker = FALSE,
                            marker_col = NULL,
                            marker_size = NULL,
                            save_image = FALSE,
                            prefix = "proc_",
                            dir_original = NULL,
                            dir_processed = NULL,
                            verbose = TRUE){
  check_ebi()
  if(!object_size %in% c("small", "medium", "large", "elarge")){
    stop("'object_size' must be one of 'small', 'medium', 'large', or 'elarge'")
  }
  if(!missing(img) & !missing(pattern)){
    stop("Only one of `img` or `pattern` arguments can be used.", call. = FALSE)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  help_count <-
    function(img, foreground, background, resize, fill_hull, threshold, filter, tolerance, extension,
             randomize, nrows, show_image, show_original, show_background, marker,
             marker_col, marker_size, save_image, prefix,
             dir_original, dir_processed, verbose){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      if(trim != FALSE){
        if(!is.numeric(trim)){
          stop("Argument `trim` must be numeric.", call. = FALSE)
        }
        img <- image_trim(img, trim)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          stop("Argument `resize` must be numeric.", call. = FALSE)
        }
        img <- image_resize(img, resize)
      }
      if(filter != FALSE){
        if(!is.numeric(filter)){
          stop("Argument `filter` must be numeric.", call. = FALSE)
        }
        img <- image_filter(img, size = filter)
      }
      if(!is.null(foreground) && !is.null(background)){
        if(is.character(foreground)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = foreground)
          check_names_dir(foreground, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          foreground <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        if(is.character(background)){
          all_files <- sapply(list.files(diretorio_original), file_name)
          imag <- list.files(diretorio_original, pattern = background)
          check_names_dir(background, all_files, diretorio_original)
          name <- file_name(imag)
          extens <- file_extension(imag)
          background <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
        }
        original <-
          data.frame(CODE = "img",
                     R = c(img@.Data[,,1]),
                     G = c(img@.Data[,,2]),
                     B = c(img@.Data[,,3]))
        foreground <-
          data.frame(CODE = "foreground",
                     R = c(foreground@.Data[,,1]),
                     G = c(foreground@.Data[,,2]),
                     B = c(foreground@.Data[,,3]))
        background <-
          data.frame(CODE = "background",
                     R = c(background@.Data[,,1]),
                     G = c(background@.Data[,,2]),
                     B = c(background@.Data[,,3]))
        back_fore <-
          transform(rbind(foreground[sample(1:nrow(foreground)),][1:nrows,],
                          background[sample(1:nrow(background)),][1:nrows,]),
                    Y = ifelse(CODE == "background", 0, 1))
        modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = back_fore))
        pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
        foreground_background <- matrix(pred1, ncol = dim(img)[[2]])
        foreground_background <- image_correct(foreground_background, perc = 0.02)
        ID <- c(foreground_background == 1)
        ID2 <- c(foreground_background == 0)
        if(isTRUE(watershed)){
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
          res <- length(foreground_background)
          parms2 <- parms[parms$object_size == object_size,]
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          nmask <- EBImage::watershed(EBImage::distmap(foreground_background),
                                      tolerance = tol,
                                      ext = ext)
        } else{
          nmask <- EBImage::bwlabel(foreground_background)
        }
      } else{
        img2 <- image_binary(img,
                             index = index,
                             my_index = my_index,
                             invert = invert,
                             fill_hull = fill_hull,
                             threshold = threshold,
                             resize = FALSE,
                             show_image = FALSE)[[1]]
        if(isTRUE(watershed)){
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
          res <- length(img2)
          parms2 <- parms[parms$object_size == object_size,]
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          nmask <- EBImage::watershed(EBImage::distmap(img2),
                                      tolerance = tol,
                                      ext = ext)
        } else{
          nmask <- EBImage::bwlabel(img2)
        }
        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
      }
      shape <-
        cbind(data.frame(EBImage::computeFeatures.shape(nmask)),
              data.frame(EBImage::computeFeatures.moment(nmask))
        )
      object_contour <- EBImage::ocontour(nmask)
      ch <- conv_hull(object_contour)
      area_ch <- trunc(as.numeric(unlist(poly_area(ch))))
      shape <- transform(shape,
                         id = 1:nrow(shape),
                         radius_ratio = s.radius.max / s.radius.min,
                         diam_mean = s.radius.mean * 2,
                         diam_min = s.radius.min * 2,
                         diam_max = s.radius.max * 2,
                         area_ch =   area_ch,
                         solidity = s.area / area_ch,
                         circularity = 4*pi*(s.area / s.perimeter^2),
                         minor_axis = m.majoraxis*sqrt(1-m.eccentricity^2))
      shape <- shape[, c("id", "m.cx", "m.cy", "s.area", "area_ch", "s.perimeter", "s.radius.mean",
                         "s.radius.min", "s.radius.max", "s.radius.sd", "radius_ratio", "diam_mean",
                         "diam_min", "diam_max", "m.majoraxis", "minor_axis", "m.eccentricity",
                         "m.theta", "solidity",  "circularity")]
      colnames(shape) <- c("id", "x", "y", "area", "area_ch", "perimeter", "radius_mean",
                           "radius_min", "radius_max", "radius_sd", "radius_ratio", "diam_mean",
                           "diam_min", "diam_max", "major_axis", "minor_axis", "eccentricity",
                           "theta", "solidity", "circularity")
      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        stop("Only one of 'lower_*' or 'topn_*' can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$area > lower_size, ],
             shape <- shape[shape$area > mean(shape$area) * 0.1, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$area, decreasing = TRUE),][1:topn_upper,]
      }
      if(!is.null(lower_eccent)){
        shape <- shape[shape$eccentricity > lower_eccent, ]
      }
      if(!is.null(upper_eccent)){
        shape <- shape[shape$eccentricity < upper_eccent, ]
      }
      if(!is.null(lower_circ)){
        shape <- shape[shape$circularity > lower_circ, ]
      }
      if(!is.null(upper_circ)){
        shape <- shape[shape$circularity < upper_circ, ]
      }
      object_contour <- object_contour[shape$id]
      if(!is.null(object_index)){
        if(!is.character(object_index)){
          stop("`object_index` must be a character.", call. = FALSE)
        }
        ind_formula <- object_index
        data_mask <- nmask@.Data
        get_rgb <- function(img, data_mask, index){
          data.frame(id = index,
                     R = img@.Data[,,1][which(data_mask == index)],
                     G = img@.Data[,,2][which(data_mask == index)],
                     B = img@.Data[,,3][which(data_mask == index)])
        }
        if(isTRUE(parallel)){
          nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
          clust <- makeCluster(nworkers)
          clusterExport(clust,
                        varlist = c("img", "data_mask", "get_rgb"),
                        envir=environment())
          on.exit(stopCluster(clust))
          object_rgb <-
            do.call(rbind,
                    parLapply(clust, 1:max(data_mask),
                              function(i){
                                get_rgb(img, data_mask, i)
                              })
            )
        } else{
          object_rgb <-
            do.call(rbind,
                    lapply(1:max(data_mask), function(i){
                      get_rgb(img, data_mask, i)
                    }))
        }
        object_rgb <- subset(object_rgb, id %in% shape$id)
        # indexes by id
        indexes <-
          by(object_rgb,
             INDICES = object_rgb$id,
             FUN = function(x){
               data.frame(
                 do.call(cbind,
                         lapply(seq_along(ind_formula), function(i){
                           data.frame(transform(x, index = eval(parse(text = ind_formula[i])))[,5])
                         })
                 )
               )
             }
          )
        indexes <-
          do.call(rbind,
                  lapply(indexes, data.frame)
          )
        indexes <- data.frame(cbind(id = object_rgb$id, indexes))
        colnames(indexes) <- c("id", ind_formula)
        indexes <- aggregate(. ~ id, indexes, mean, na.rm = TRUE)
      } else{
        object_rgb <- NULL
        indexes <- NULL
      }
      stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area",
                                   "sd_area", "sum_area"),
                          value = c(length(shape$area),
                                    min(shape$area),
                                    mean(shape$area),
                                    max(shape$area),
                                    sd(shape$area),
                                    sum(shape$area)))
      results <- list(results = shape,
                      statistics = stats,
                      object_rgb = object_rgb,
                      object_index = indexes)
      class(results) <- "anal_obj"
      if(show_image == TRUE | save_image == TRUE){
        backg <- !is.null(col_background)
        col_background <- col2rgb(ifelse(is.null(col_background), "white", col_background))
        col_foreground <- col2rgb(ifelse(is.null(col_foreground), "black", col_foreground))
        if(show_original == TRUE & show_segmentation == FALSE){
          im2 <- img
          if(backg){
            im3 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- col_background[3]
          }
        }
        if(show_original == TRUE & show_segmentation == TRUE){
          im2 <- EBImage::colorLabels(nmask)
          if(backg){
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
          } else{
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- img@.Data[,,1][which(im2@.Data[,,1]==0)]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- img@.Data[,,2][which(im2@.Data[,,2]==0)]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- img@.Data[,,3][which(im2@.Data[,,3]==0)]
          }
        }
        if(show_original == FALSE){
          if(show_segmentation == TRUE){
            im2 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
          } else{
            im2 <- img
            im2@.Data[,,1][ID] <- col_foreground[1]
            im2@.Data[,,2][ID] <- col_foreground[2]
            im2@.Data[,,3][ID] <- col_foreground[3]
            im2@.Data[,,1][ID2] <- col_background[1]
            im2@.Data[,,2][ID2] <- col_background[2]
            im2@.Data[,,3][ID2] <- col_background[3]
          }
        }
        show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
        marker <- ifelse(is.null(marker), "id", marker)
        if(!isFALSE(show_mark) & marker != "point" & !marker %in% colnames(shape)){
          warning("Accepted 'marker' are: {", paste(colnames(shape), collapse = ", "),
                  "}. Drawing the object id.", call. = FALSE)
          marker <- "id"
        }
        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
        if(show_image == TRUE){
          if(marker != "point"){
            plot(im2)
            if(show_mark){
              text(shape[, 2],
                   shape[, 3],
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(isTRUE(show_chull)){
              plot_contour(ch, col = "black")
            }
          } else{
            plot(im2)
            if(show_mark){
              points(shape[, 2],
                     shape[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
            if(isTRUE(show_contour)  & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
          }
        }
        if(save_image == TRUE){
          if(dir.exists(diretorio_processada) == FALSE){
            dir.create(diretorio_processada, recursive = TRUE)
          }
          png(paste0(diretorio_processada, "/",
                     prefix,
                     name_ori, ".",
                     extens_ori),
              width = dim(im2@.Data)[1],
              height = dim(im2@.Data)[2])
          if(marker != "point"){
            plot(im2)
            if(show_mark){
              text(shape[, 2],
                   shape[, 3],
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
          } else{
            plot(im2)
            if(show_mark){
              points(shape[, 2],
                     shape[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
          }
          dev.off()
        }
      }
      invisible(results)
    }
  if(missing(pattern)){
    help_count(img, foreground, background, resize, fill_hull, threshold, filter,
               tolerance , extension, randomize, nrows, show_image, show_original,
               show_background, marker, marker_col, marker_size, save_image, prefix,
               dir_original, dir_processed, verbose)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("Pattern '", pattern, "' not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_count(x,
                               foreground, background, resize, fill_hull, threshold,
                               filter, tolerance , extension, randomize,
                               nrows, show_image, show_original, show_background,
                               marker, marker_col, marker_size, save_image, prefix,
                               dir_original, dir_processed, verbose =  FALSE)
                  })

    } else{
      pb <- progress(max = length(plants), style = 4)
      foo <- function(plants, ...){
        run_progress(pb, ...)
        help_count(img  = plants,
                   foreground, background, resize, fill_hull, threshold, filter,
                   tolerance, extension, randomize, nrows, show_image, show_original,
                   show_background, marker, marker_col, marker_size, save_image,
                   prefix, dir_original, dir_processed, verbose)
      }
      results <-
        lapply(seq_along(names_plant), function(i){
          foo(names_plant[i],
              actual = i,
              text = paste("Processing image", names_plant[i]))
        })
    }
    names(results) <- names_plant
    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["statistics"]],
                          id =  names(results[i]))[,c(3, 1, 2)]
              })
      )
    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["results"]],
                          img =  names(results[i]))
              })
      )
    if("img" %in% colnames(results)){
      results <- results[, c(ncol(results), 1:ncol(results) - 1)]
    }
    summ <- stats[stats$stat == "n", c(1, 3)]
    if(verbose == TRUE){
      names(summ) <- c("Image", "Objects")
      cat("--------------------------------------------\n")
      print(summ, row.names = FALSE)
      cat("--------------------------------------------\n")
      message("Done!")

    }
    invisible(
      structure(
        list(statistics = stats,
             count = summ,
             results = results),
        class = "anal_obj_ls"
      )
    )
  }
}


#' @name analyze_objects
#' @param x An object of class `anal_obj`.
#' @param which Which to plot. Either 'measure' (object measures) or 'index'
#'   (object index). Defaults to `"measure"`.
#' @param measure The measure to plot. Defaults to `"area"`.
#' @param type The type of plot. Either `"hist"` or `"density"`. Partial matches
#'   are recognized.
#' @param facet Create a facet plot for each object when `which = "index"` is
#'   used?. Defaults to `FALSE`.
#' @param ... Further argument passed on to [lattice::histogram()] or
#'   [lattice::densityplot()]
#' @method plot anal_obj
#' @importFrom lattice densityplot levelplot
#' @export
#' @return `plot.anal_obj()` returns a `trellis` object containing the
#'   distribution of the pixels, optionally  for each object when `facet = TRUE`
#'   is used.
#' @examples
#' \donttest{
#' library(pliman)
#'
#' img <- image_pliman("soy_green.jpg")
#' # Segment the foreground (grains) using the normalized blue index (NB, default)
#' # Shows the average value of the blue index in each object
#'
#' rgb <-
#'    analyze_objects(img,
#'                    marker = "id",
#'                    object_index = "B")
#' # density of area
#' plot(rgb)
#'
#' # histogram of perimeter
#' plot(rgb, measure = "perimeter", type = "histogram") # or 'hist'
#'
#' # density of the blue (B) index
#' plot(rgb, which = "index")
#' }
plot.anal_obj <- function(x,
                          which = "measure",
                          measure = "area",
                          type = "density",
                          facet = FALSE,
                          ...){
  if(!which %in% c("measure", "index")){
    stop("'which' must be one of 'measure' or 'index'", call. = FALSE)
  }
  if(which == "measure"){
    nam <- colnames(x$results)
    if(!measure %in% nam){
      stop("Measure '", measure, "' not available in 'x'. Try one of the '",
           paste0(nam, collapse = ", "), call. = FALSE)
    }
    temp <- x$results[[measure]]
    types <- c("density", "histogram")
    matches <- grepl(type, types)
    type <- types[matches]
    if(type == "histogram"){
      lattice::histogram(temp,
                         type = "count",
                         xlab = paste(measure, "(pixels)"),
                         ...)
    } else{
      lattice::densityplot(temp,
                           xlab = paste(measure, "(pixels)"),
                           col = "blue",
                           ...)
    }
  } else{
    rgb <- x$object_rgb
    if(is.null(rgb)){
      stop("RGB values not found. Use `object_index` in the function `analyze_objects()`.", call. = FALSE)
    }
    rgb$id <- rownames(rgb)
    rgb <-
      reshape(rgb,
              direction = "long",
              varying = list(names(rgb)[2:4]),
              v.names = "value",
              idvar = "id",
              timevar = "Spectrum",
              times = c("r", "g", "b"))
    rgb$Spectrum <- factor( rgb$Spectrum, levels = unique( rgb$Spectrum))
    if(isTRUE(facet)){
      densityplot(~value | factor(id),
                  data = rgb,
                  groups = Spectrum,
                  par.settings = list(superpose.line = list(col = c("red", "green","blue"))),
                  xlab = "Pixel value",
                  plot.points = FALSE,
                  ...)
    } else{
      densityplot(~value,
                  data = rgb,
                  groups = Spectrum,
                  par.settings = list(superpose.line = list(col = c("red", "green","blue"))),
                  xlab = "Pixel value",
                  plot.points = FALSE,
                  ...)
    }
  }
}
