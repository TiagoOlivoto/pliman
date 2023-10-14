#' Performs plant disease measurements
#'
#' @description
#'* `measure_disease()` computes the percentage of symptomatic leaf area and
#'(optionally) counts and compute shapes (area, perimeter, radius, etc.) of
#'lesions in a sample or entire leaf using color palettes. See more at
#'**Details**.
#'
#'* `measure_disease_iter()` provides an iterative section for
#'`measure_disease()`, where the user picks up samples in the image to create
#'the needed color palettes.
#'
#' @details
#'
#' In `measure_disease()`, a general linear model (binomial family) fitted to
#'the RGB values is used to segment the lesions from the healthy leaf. If a
#'pallet of background is provided, the function takes care of the details to
#'isolate it before computing the number and area of lesions. By using `pattern`
#'it is possible to process several images with common pattern names that are
#'stored in the current working directory or in the subdirectory informed in
#'`dir_original`.
#'
#'If `img_healthy` and `img_symptoms` are not declared, RGB-based phenotyping of
#'foliar disease severity is performed using the index informed in `index_lb` to
#'first segment leaf from background and `index_dh` to segment diseased from
#'healthy tissues.
#'
#'`measure_disease_iter()` only run in an interactive section. In this function,
#'users will be able to pick up samples of images to iteratively create the
#'needed color palettes. This process calls [`pick_palette()`] internally. If
#'`has_background` is TRUE (default) the color palette for the background is
#'first created. The sample of colors is performed in each left-button mouse
#'click and continues until the user press Esc. Then, a new sampling process is
#'performed to sample the color of healthy tissues and then diseased tissues.
#'The generated palettes are then passed on to measure_disease(). All the
#'arguments of such function can be passed using the ... (three dots).
#'
#' When `show_features = TRUE`, the function computes a total of 36 lesion
#' features (23 shape features and 13 texture features). The Haralick texture
#' features for each object based on a gray-level co-occurrence matrix (Haralick
#' et al. 1979). See more details in [analyze_objects()].
#'
#' @param img The image to be analyzed.
#' @param img_healthy A color palette of healthy tissues.
#' @param img_symptoms A color palette of lesioned tissues.
#' @param img_background A color palette of the background (if exists). These
#'   arguments can be either an `Image` object stored in the global environment
#'   or a character value. If a chacarceter is used (eg., `img_healthy =
#'   "leaf"`), the function will search in the current working directory a valid
#'   image that contains "`leaf"` in the name. Note that if two images matches
#'   this pattern, an error will occour.
#' @param pattern A pattern of file name used to identify images to be
#'   processed. For example, if `pattern = "im"` all images that the name
#'   matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
#'   Providing any number as pattern (e.g., `pattern = "1"`) will select
#'   images that are named as 1.-, 2.-, and so on.
#' @param filter Performs median filtering in the binary image that segments the
#'   leaf from background? By default, a median filter of `size = 10` is applied.
#'   This is useful to reduce the noise and segment the leaf and background more
#'   accurately. See more at [image_filter()]. Set to `FALSE` to cancel median
#'   filtering.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time, especially when `pattern` is used is informed. The
#'   number of sections is set up to 30% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'   numeric value of range 0-100 (proportion of the size of the original
#'   image).
#' @param fill_hull Fill holes in the image? Defaults to `TRUE`. This is useful
#'   to fill holes in leaves, e.g., those caused by insect attack, ensuring the
#'   hole area will be accounted for the leaf, not background.
#' @param index_lb The index used to segment the foreground (e.g., leaf) from
#'   the background. If not declared, the entire image area (pixels) will be
#'   considered in the computation of the severity.
#' @param index_dh The index used to segment diseased from healthy tissues when
#'   `img_healthy` and `img_symptoms` are not declared. Defaults to `"GLI"`. See
#'   [image_index()] for more details.
#' @param has_white_bg Logical indicating whether a white background is present.
#'   If `TRUE`, pixels that have R, G, and B values equals to 1 will be
#'   considered as `NA`. This may be useful to compute an image index for
#'   objects that have, for example, a white background. In such cases, the
#'   background will not be considered for the threshold computation.
#' @param threshold By default (`threshold = NULL`), a threshold value based on
#'   Otsu's method is used to reduce the grayscale image to a binary image. If a
#'   numeric value is informed, this value will be used as a threshold. Inform
#'   any non-numeric value different than "Otsu" to iteratively choose the
#'   threshold based on a raster plot showing pixel intensity of the index. Must
#'   be a vector of length 2 to indicate the threshold for `index_lb` and
#'   `index_dh`, respectively.
#' @param invert Inverts the binary image if desired. This is useful to process
#'   images with black background. Defaults to `FALSE`.
#' @param lower_size Lower limit for size for the image analysis. Leaf images
#'   often contain dirt and dust. To prevent dust from affecting the image
#'   analysis, the lower limit of analyzed size is set to 0.1, i.e., objects
#'   with lesser than 10% of the mean of all objects are removed. One can set a
#'   known area or use `lower_limit = 0` to select all objects (not advised).
#' @param upper_size Upper limit for size for the image analysis. Defaults to
#'   `NULL`, i.e., no upper limit used.
#' @param topn_lower,topn_upper Select the top `n` lesions based on its area.
#'   `topn_lower` selects the `n` lesions with the smallest area whereas
#'   `topn_upper` selects the `n` lesions with the largest area.
#' @param randomize Randomize the lines before training the model? Defaults to
#'   `TRUE`.
#' @param nsample The number of sample pixels to be used in training step.
#'   Defaults to `3000`.
#' @param watershed If `TRUE` (Default) implements the Watershed Algorithm to
#'   segment lesions connected by a fairly few pixels that could be considered
#'   as two distinct lesions. If `FALSE`, lesions that are connected by any
#'   pixel are considered unique lesions. For more details see
#'   [EBImage::watershed()].
#' @param lesion_size The size of the lesion. Used to automatically tune
#'   `tolerance` and `extension` parameters. One of the following. `"small"`
#'   (2-5 mm in diameter, e.g, rust pustules), `"medium"` (0.5-1.0 cm in
#'   diameter, e.g, wheat leaf spot), `"large"` (1-2 cm in diameter, and
#'   `"elarge"` (2-3 cm in diameter, e.g, target spot of soybean).
#' @param tolerance The minimum height of the object in the units of image
#'   intensity between its highest point (seed) and the point where it contacts
#'   another object (checked for every contact pixel). If the height is smaller
#'   than the tolerance, the object will be combined with one of its neighbors,
#'   which is the highest. Defaults to `NULL`, i.e., starting values are set up
#'   according to the argument `lesion_size`.
#' @param extension Radius of the neighborhood in pixels for the detection of
#'   neighboring objects. Defaults to 20. Higher value smooths out small
#'   objects.
#' @param show_features If `TRUE` returnS the lesion features such as number,
#'   area, perimeter, and radius. Defaults to `FALSE`.
#' @param show_segmentation Shows the object segmentation colored with random
#'   permutations. Defaults to `TRUE`.
#' @param plot Show image after processing? Defaults to `TRUE`.
#' @param show_original Show the symptoms in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param show_contour Show a contour line around the lesions? Defaults
#'   to `TRUE`.
#' @param contour_col,contour_size The color and size for the contour line
#'   around objects. Defaults to `contour_col = "white"` and `contour_size = 1`.
#' @param col_leaf Leaf color after image processing. Defaults to `"green"`
#' @param col_lesions Symptoms color after image processing. Defaults to
#'   `"red"`.
#' @param col_background Background color after image processing. Defaults to
#'   `"NULL"`.
#' @param marker,marker_col,marker_size The type, color and size of the object
#'   marker. Defaults to `NULL`, which shows nothing. Use `marker = "point"` to
#'   show a point in each lesion or `marker = "*"` where `"*"` is any variable
#'   name of the `shape` data frame returned by the function.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named as `proc_*` where `*` is the image name
#'   given in `img`.
#' @param prefix The prefix to be included in the processed images. Defaults to
#'   `"proc_"`.
#' @param name The name of the image to save. Use this to overwrite the name of
#'   the image in `img`.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory. It can be either a full path, e.g.,
#'   `"C:/Desktop/imgs"`, or a subfolder within the current working directory,
#'   e.g., `"/imgs"`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param has_background A logical indicating if the image has a background to
#'   be segmented before processing.
#' @param r The radius of neighborhood pixels. Defaults to `2`. A square is
#'   drawn indicating the selected pixels.
#' @param viewer The viewer option. If not provided, the value is retrieved
#'   using [get_pliman_viewer()]. This option controls the type of viewer to use
#'   for interactive plotting. The available options are "base" and "mapview".
#'   If set to "base", the base R graphics system is used for interactive
#'   plotting. If set to "mapview", the mapview package is used. To set this
#'   argument globally for all functions in the package, you can use the
#'   [set_pliman_viewer()] function. For example, you can run
#'   `set_pliman_viewer("mapview")` to set the viewer option to "mapview" for
#'   all functions.
#' @param show The show option for the mapview viewer, either `"rgb"` or
#'   `"index"`.
#' @param index The index to be shown when `show = "rgb"`.
#' @param ... Further parameters passed on to `measure_disease()`.
#' @return
#' * `measure_disease()` returns a list with the following objects:
#'   - `severity` A data frame with the percentage of healthy and symptomatic
#'  areas.
#'   - `shape`,`statistics` If `show_features = TRUE` is used, returns the shape
#'  (area, perimeter, etc.) for each lesion and a summary statistic of the
#'  results.
#' * `measure_disease_iter()` returns a list with the following objects:
#'    - `results` A list with the objects returned by `measure_disease()`.
#'    - `leaf` The color palettes for the healthy leaf.
#'    - `disease` The color palettes for the diseased leaf.
#'    - `background` The color palettes for the background.
#' @name measure_disease
#' @export
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("sev_leaf_nb.jpg")
#' healthy <- image_pliman("sev_healthy.jpg")
#' lesions <- image_pliman("sev_sympt.jpg")
#' image_combine(img, healthy, lesions, ncol = 3)
#'
#' sev <-
#'  measure_disease(img = img,
#'                  img_healthy = healthy,
#'                  img_symptoms = lesions,
#'                  lesion_size = "large",
#'                  plot = TRUE)
#'
#' # an interactive section
#' measure_disease_iter(img)
#' }
#'
measure_disease <- function(img,
                            img_healthy = NULL,
                            img_symptoms = NULL,
                            img_background = NULL,
                            pattern = NULL,
                            filter = 10,
                            parallel = FALSE,
                            workers = NULL,
                            resize = FALSE,
                            fill_hull = TRUE,
                            index_lb = NULL,
                            index_dh = "GLI",
                            has_white_bg = FALSE,
                            threshold = NULL,
                            invert = FALSE,
                            lower_size = NULL,
                            upper_size = NULL,
                            topn_lower = NULL,
                            topn_upper = NULL,
                            randomize = TRUE,
                            nsample = 3000,
                            watershed = FALSE,
                            lesion_size = "medium",
                            tolerance = NULL,
                            extension = NULL,
                            show_features = FALSE,
                            show_segmentation = FALSE,
                            plot = TRUE,
                            show_original = TRUE,
                            show_background = TRUE,
                            show_contour = TRUE,
                            contour_col = "white",
                            contour_size = 1,
                            col_leaf = NULL,
                            col_lesions = NULL,
                            col_background = NULL,
                            marker = FALSE,
                            marker_col = NULL,
                            marker_size = NULL,
                            save_image = FALSE,
                            prefix = "proc_",
                            name = NULL,
                            dir_original = NULL,
                            dir_processed = NULL,
                            verbose = TRUE){
  # check_ebi()
  if(!missing(img) & !missing(pattern)){
    stop("Only one of `img` or `pattern` arguments can be used.", call. = FALSE)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste("./", sep = "")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  if(is.character(img_healthy)){
    all_files <- sapply(list.files(diretorio_original), file_name)
    imag <- list.files(diretorio_original, pattern = img_healthy)
    check_names_dir(img_healthy, all_files, "")
    name_h <- file_name(imag)
    extens <- file_extension(imag)
    img_healthy <- image_import(paste(diretorio_original, "/", name_h, ".", extens, sep = ""))
  }
  if(is.character(img_symptoms)){
    all_files <- sapply(list.files(diretorio_original), file_name)
    imag <- list.files(diretorio_original, pattern = img_symptoms)
    check_names_dir(img_symptoms, all_files, "")
    name_h <- file_name(imag)
    extens <- file_extension(imag)
    img_symptoms <- image_import(paste(diretorio_original, "/", name_h, ".", extens, sep = ""))
  }
  if(is.character(img_background)){
    all_files <- sapply(list.files(diretorio_original), file_name)
    imag <- list.files(diretorio_original, pattern = img_background)
    check_names_dir(img_background, all_files, "")
    name_h <- file_name(imag)
    extens <- file_extension(imag)
    img_background <- image_import(paste(diretorio_original, "/", name_h, ".", extens, sep = ""))
  }
  help_count <-
    function(img, img_healthy, img_symptoms, img_background, resize, fill_hull, invert,
             index_lb, index_dh, has_white_bg, lesion_size, tolerance, extension,
             randomize, nsample, plot, show_original, show_background,
             col_leaf, col_lesions, col_background,
             save_image, dir_original, dir_processed, marker, marker_col, marker_size){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "jpg"
      }
      backg <- !is.null(col_background)
      # color for background
      if (is.null(col_background)){
        col_background <- col2rgb("white") / 255
      } else{
        ifelse(is.character(col_background),
               col_background <- col2rgb(col_background) / 255,
               col_background <- col_background / 255)
      }
      # color for lesions
      if (is.null(col_lesions)){
        col_lesions <- col2rgb("black") / 255
      } else{
        ifelse(is.character(col_lesions),
               col_lesions <- col2rgb(col_lesions) / 255,
               col_lesions <- col_lesions / 255)
      }
      # color for leaf
      if (is.null(col_leaf)){
        col_leaf <- col2rgb("green") / 255
      } else{
        ifelse(is.character(col_leaf),
               col_leaf <- col2rgb(col_leaf) / 255,
               col_leaf <- col_leaf / 255)
      }
      if(!is.null(img_healthy) && !is.null(img_symptoms)){

        original <-
          data.frame(CODE = "img",
                     R = c(img@.Data[,,1]),
                     G = c(img@.Data[,,2]),
                     B = c(img@.Data[,,3]))
        sadio <-
          data.frame(CODE = "img_healthy",
                     R = c(img_healthy@.Data[,,1]),
                     G = c(img_healthy@.Data[,,2]),
                     B = c(img_healthy@.Data[,,3]))
        sintoma <-
          data.frame(CODE = "img_symptoms",
                     R = c(img_symptoms@.Data[,,1]),
                     G = c(img_symptoms@.Data[,,2]),
                     B = c(img_symptoms@.Data[,,3]))
        ncol_img <- dim(img)[[2]]
        ################## no background #############
        if(is.null(img_background)){
          sadio_sintoma <-
            transform(rbind(sadio[sample(1:nrow(sadio)),][1:nsample,],
                            sintoma[sample(1:nrow(sintoma)),][1:nsample,]),
                      Y = ifelse(CODE == "img_healthy", 1, 0))
          sadio_sintoma$CODE <- NULL
          usef_area <- nrow(original)
          model <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = sadio_sintoma))
          # isolate plant
          pred1 <- round(predict(model, newdata = original, type="response"), 0)
          plant_symp <- 1 - matrix(pred1, ncol = ncol_img)
          ID <- c(plant_symp == 0)
          pix_sympt <- length(which(ID == FALSE))
          ID2 <- c(plant_symp == 1)
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE),
                            header = T, sep = ";")
          parms2 <- parms[parms$object_size == lesion_size,]
          res <- length(plant_symp)
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          ifelse(watershed == FALSE,
                 nmask <- EBImage::bwlabel(plant_symp),
                 nmask <- EBImage::watershed(EBImage::distmap(plant_symp),
                                             tolerance = tol,
                                             ext = ext))
          if(plot == TRUE | save_image == TRUE){
            if(show_original == TRUE & show_segmentation == FALSE){
              im2 <- img
              if(isFALSE(show_contour)){
                im2@.Data[,,1][!ID] <- col_lesions[1]
                im2@.Data[,,2][!ID] <- col_lesions[2]
                im2@.Data[,,3][!ID] <- col_lesions[3]
              }
              if(backg){
                im3 <- EBImage::colorLabels(nmask)
                im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- img@.Data[,,1][which(im3@.Data[,,1]==0)]
                im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- img@.Data[,,2][which(im3@.Data[,,2]==0)]
                im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- img@.Data[,,3][which(im3@.Data[,,3]==0)]
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
                im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_leaf[1]
                im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_leaf[2]
                im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_leaf[3]
              } else{
                im2 <- img
                im2@.Data[,,1][!ID] <- col_lesions[1]
                im2@.Data[,,2][!ID] <- col_lesions[2]
                im2@.Data[,,3][!ID] <- col_lesions[3]
                im2@.Data[,,1][ID] <- col_leaf[1]
                im2@.Data[,,2][ID] <- col_leaf[2]
                im2@.Data[,,3][ID] <- col_leaf[3]
              }
            }
          }
        } else{

          fundo <-
            data.frame(CODE = "img_background",
                       R = c(img_background@.Data[,,1]),
                       G = c(img_background@.Data[,,2]),
                       B = c(img_background@.Data[,,3]))
          # separate image from background
          fundo_resto <-
            transform(rbind(sadio[sample(1:nrow(sadio)),][1:nsample,],
                            sintoma[sample(1:nrow(sintoma)),][1:nsample,],
                            fundo[sample(1:nrow(fundo)),][1:nsample,]),
                      Y = ifelse(CODE == "img_background", 0, 1))
          modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"),
                                          data = fundo_resto))
          pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
          ifelse(fill_hull == TRUE,
                 plant_background <- EBImage::fillHull(matrix(pred1, ncol = ncol_img)),
                 plant_background <- matrix(pred1, ncol = ncol_img))
          if(is.numeric(filter) & filter > 1){
            plant_background <- EBImage::medianFilter(plant_background, size = filter)
          }
          plant_background[plant_background == 1] <- 2
          sadio_sintoma <-
            transform(rbind(sadio[sample(1:nrow(sadio)),][1:nsample,],
                            sintoma[sample(1:nrow(sintoma)),][1:nsample,]),
                      Y = ifelse(CODE == "img_healthy", 1, 0))
          sadio_sintoma$CODE <- NULL
          modelo2 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"),
                                          data = sadio_sintoma))
          # isolate plant
          ID <- c(plant_background == 2)
          usef_area <- nrow(original[ID,])
          pred2 <- round(predict(modelo2, newdata = original[ID,], type="response"), 0)
          pred3 <- round(predict(modelo2, newdata = original, type="response"), 0)
          pix_sympt <- length(which(pred2 == 0))
          pred3[!ID] <- 1
          leaf_sympts <- 1 - matrix(pred3, ncol = ncol_img)
          plant_background[leaf_sympts == 1] <- 3
          parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE),
                            header = T, sep = ";")
          parms2 <- parms[parms$object_size == lesion_size,]
          res <- length(leaf_sympts)
          rowid <-
            which(sapply(as.character(parms2$resolution), function(x) {
              eval(parse(text=x))}))
          ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
          tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
          if(isTRUE(fill_hull)){
            leaf_sympts <- EBImage::fillHull(leaf_sympts)
          }
          ifelse(watershed == FALSE,
                 nmask <- EBImage::bwlabel(leaf_sympts),
                 nmask <- EBImage::watershed(EBImage::distmap(leaf_sympts),
                                             tolerance = tol,
                                             ext = ext))
          if(plot == TRUE | save_image == TRUE){
            if(show_original == TRUE & show_segmentation == TRUE){
              im2 <- EBImage::colorLabels(nmask)
              if(backg){
                im2@.Data[,,1][!ID] <- col_background[1]
                im2@.Data[,,2][!ID] <- col_background[2]
                im2@.Data[,,3][!ID] <- col_background[3]
                im2@.Data[,,1][ID][which(pred2 != 0)] <- img@.Data[,,1][ID][which(pred2 != 0)]
                im2@.Data[,,2][ID][which(pred2 != 0)] <- img@.Data[,,2][ID][which(pred2 != 0)]
                im2@.Data[,,3][ID][which(pred2 != 0)] <- img@.Data[,,3][ID][which(pred2 != 0)]
              } else{
                im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- img@.Data[,,1][which(im2@.Data[,,1]==0)]
                im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- img@.Data[,,2][which(im2@.Data[,,2]==0)]
                im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- img@.Data[,,3][which(im2@.Data[,,3]==0)]
              }
            }
            if(show_original == TRUE & show_segmentation == FALSE){
              im2 <- img
              if(isFALSE(show_contour)){
                im2@.Data[,,1][ID][which(pred2 == 0)] <- col_lesions[1]
                im2@.Data[,,2][ID][which(pred2 == 0)] <- col_lesions[2]
                im2@.Data[,,3][ID][which(pred2 == 0)] <- col_lesions[3]
              }
              if(backg){
                im2@.Data[,,1][!ID] <- col_background[1]
                im2@.Data[,,2][!ID] <- col_background[2]
                im2@.Data[,,3][!ID] <- col_background[3]
              }
            }
            if(show_original == FALSE){
              if(show_segmentation == TRUE){
                im2 <- EBImage::colorLabels(nmask)
                im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
                im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
                im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
                im2@.Data[,,1][ID][which(pred2 != 0)] <- col_leaf[1]
                im2@.Data[,,2][ID][which(pred2 != 0)] <- col_leaf[2]
                im2@.Data[,,3][ID][which(pred2 != 0)] <- col_leaf[3]
              } else{
                im2 <- img
                im2@.Data[,,1][ID][which(pred2 == 0)] <- col_lesions[1]
                im2@.Data[,,2][ID][which(pred2 == 0)] <- col_lesions[2]
                im2@.Data[,,3][ID][which(pred2 == 0)] <- col_lesions[3]
                im2@.Data[,,1][ID][which(pred2 != 0)] <- col_leaf[1]
                im2@.Data[,,2][ID][which(pred2 != 0)] <- col_leaf[2]
                im2@.Data[,,3][ID][which(pred2 != 0)] <- col_leaf[3]
                im2@.Data[,,1][!ID] <- col_background[1]
                im2@.Data[,,2][!ID] <- col_background[2]
                im2@.Data[,,3][!ID] <- col_background[3]
              }
            }
          }
        }
      } else{

        ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        if(!is.null(index_lb)){
          # segment leaf from background
          if(is.null(threshold)){
            threshold1 <- "Otsu"
          } else{
            threshold1 <- ifelse(length(threshold) == 1, threshold, threshold[1])
          }
          my_thresh <- ifelse(is.na(suppressWarnings(as.numeric(threshold1))),
                              as.character(threshold1),
                              as.numeric(threshold1))
          if(!isFALSE(invert)){
            invert1 <- ifelse(length(invert) == 1, invert, invert[1])
          } else{
            invert1 <- FALSE
          }
          seg <- help_segment(img,
                              index = index_lb,
                              threshold = my_thresh,
                              invert = invert1,
                              filter = filter)

          img <- seg
        }
        # segment disease from leaf
        if(is.null(threshold)){
          threshold2 <- "Otsu"
        } else{
          threshold2 <- ifelse(length(threshold) == 1, threshold, threshold[2])
        }
        my_thresh2 <- ifelse(is.na(suppressWarnings(as.numeric(threshold2))),
                             as.character(threshold2),
                             as.numeric(threshold2))
        if(!isFALSE(invert)){
          invert2 <- ifelse(length(invert) == 1, invert, invert[2])
        } else{
          invert2 <- FALSE
        }
        img2 <- help_binary(img,
                            index = index_dh,
                            threshold = my_thresh2,
                            invert = invert2,
                            has_white_bg = has_white_bg,
                            resize = FALSE)
        img2@.Data[is.na(img2@.Data)] <- FALSE
        # which(is.na(img2@.Data))
        res <- length(img2)
        if(!is.null(index_lb)){
          usef_area <- res - length(which(img[,,1]==1))
          img2@.Data[which(img[,,1]==1)] <- FALSE
        } else{
          usef_area <- res
        }
        pix_sympt <- length(which(img2 == TRUE))
        parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        parms2 <- parms[parms$object_size == lesion_size,]
        rowid <-
          which(sapply(as.character(parms2$resolution), function(x) {
            eval(parse(text=x))}))
        ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
        tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
        if(isTRUE(fill_hull)){
          img2 <- EBImage::fillHull(img2)
        }
        ifelse(watershed == FALSE,
               nmask <- EBImage::bwlabel(img2),
               nmask <- EBImage::watershed(EBImage::distmap(img2),
                                           tolerance = tol,
                                           ext = ext)
        )
        ID <- which(img2 == 1)
        ID2 <- which(img2 == 0)
        if(plot == TRUE | save_image == TRUE){
          if(show_original == TRUE & show_segmentation == FALSE){
            im2 <- img
            im2@.Data[,,1][which(img[,,1]==1)] <- col_background[1]
            im2@.Data[,,2][which(img[,,2]==1)] <- col_background[2]
            im2@.Data[,,3][which(img[,,3]==1)] <- col_background[3]
            if(isFALSE(show_contour)){
              im2@.Data[,,1][ID] <- col_lesions[1]
              im2@.Data[,,2][ID] <- col_lesions[2]
              im2@.Data[,,3][ID] <- col_lesions[3]
            }
          }
          if(show_original == TRUE & show_segmentation == TRUE){
            im2 <- EBImage::colorLabels(nmask)
            if(!is.null(index_lb)){
              im2@.Data[,,1][which(img[,,1]==1)] <- col_background[1]
              im2@.Data[,,2][which(img[,,2]==1)] <- col_background[2]
              im2@.Data[,,3][which(img[,,3]==1)] <- col_background[3]
            }
            im2@.Data[,,1][ID2] <- img@.Data[,,1][ID2]
            im2@.Data[,,2][ID2] <- img@.Data[,,2][ID2]
            im2@.Data[,,3][ID2] <- img@.Data[,,3][ID2]
          }
          if(show_original == FALSE){
            if(show_segmentation == TRUE){
              im2 <- EBImage::colorLabels(nmask)
              im2@.Data[,,1][ID2] <- col_leaf[1]
              im2@.Data[,,2][ID2] <- col_leaf[2]
              im2@.Data[,,3][ID2] <- col_leaf[3]
            } else{
              im2 <- img
              im2@.Data[,,1][ID2] <- col_leaf[1]
              im2@.Data[,,2][ID2] <- col_leaf[2]
              im2@.Data[,,3][ID2] <- col_leaf[3]
              im2@.Data[,,1][ID] <- col_lesions[1]
              im2@.Data[,,2][ID] <- col_lesions[2]
              im2@.Data[,,3][ID] <- col_lesions[3]
            }
            if(!is.null(index_lb)){
              im2@.Data[,,1][which(img[,,1]==1)] <- col_background[1]
              im2@.Data[,,2][which(img[,,2]==1)] <- col_background[2]
              im2@.Data[,,3][which(img[,,3]==1)] <- col_background[3]
            }
          }
        }
      }
      symptomatic <- pix_sympt /  usef_area * 100
      healthy <- 100 - symptomatic
      severity <- data.frame(healthy = healthy,
                             symptomatic = symptomatic)
      has_any_sev <- length(which(nmask != 0)) > 0
      if(has_any_sev){
        shape <- compute_measures_disease(mask = nmask)
        has_lesion <- nrow(shape$shape) > 0
      } else{
        has_lesion <- FALSE
      }

      if(isTRUE(show_features) & has_lesion & has_any_sev){
        object_contour <- shape$cont
        ch <- shape$ch
        shape <- shape$shape
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
        stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area",
                                     "sd_area", "sum_area"),
                            value = c(length(shape$area),
                                      min(shape$area),
                                      mean(shape$area),
                                      max(shape$area),
                                      sd(shape$area),
                                      sum(shape$area)))
      } else{
        shape <- NULL
        stats <- NULL
      }
      if(!isFALSE(marker) & isTRUE(show_features)){
        show_mark <- TRUE
        marker <- ifelse(is.null(marker), "id", marker)
        if(!isFALSE(show_mark) & marker != "point" & !marker %in% colnames(shape)){
          warning("Accepted 'marker' are: {", paste(colnames(shape), collapse = ", "),
                  "}. Drawing the object id.", call. = FALSE)
          marker <- "id"
        }
        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.9, marker_size)
      } else{
        show_mark <- FALSE
      }
      if(isTRUE(show_features) & isTRUE(show_contour) & has_lesion){
        ocont <- object_contour[shape$id]
      }
      if(isTRUE(show_contour) & show_original == TRUE){
        ocont <- EBImage::ocontour(nmask)
      }
      if(plot == TRUE){
        if(marker != "point"){
          plot(im2)
          if(show_features & show_mark & has_lesion){
            text(shape[,2],
                 shape[,3],
                 round(shape[, marker], 2),
                 col = marker_col,
                 cex = marker_size)
          }
          if(isTRUE(show_contour) & show_original == TRUE){
            plot_contour(ocont, col = contour_col, lwd = contour_size)
          }
        } else{
          plot(im2)
          if(show_features & show_mark & has_lesion){
            points(shape[,2],
                   shape[,3],
                   col = marker_col,
                   pch = 16,
                   cex = marker_size)
          }
          if(isTRUE(show_contour) & show_original == TRUE){
            plot_contour(ocont, col = contour_col, lwd = contour_size)
          }
        }
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada)
        }
        name_img <- ifelse(is.null(name), name_ori, name)
        jpeg(paste0(diretorio_processada, "/",
                    prefix,
                    name_img, ".",
                    "jpg"),
             width = dim(im2@.Data)[1],
             height = dim(im2@.Data)[2])
        if(marker != "point"){
          plot(im2)
          if(show_features & show_mark & has_lesion){
            text(shape[,2],
                 shape[,3],
                 round(shape[, marker], 2),
                 col = marker_col,
                 cex = marker_size)
          }
          if(isTRUE(show_contour) & show_original == TRUE){
            plot_contour(ocont, col = contour_col, lwd = contour_size)
          }
        } else{
          plot(im2)
          if(show_features & show_mark & has_lesion){
            points(shape[,2],
                   shape[,3],
                   col = marker_col,
                   pch = 16,
                   cex = marker_size)
          }
          if(isTRUE(show_contour) & show_original == TRUE){
            plot_contour(ocont, col = contour_col, lwd = contour_size)
          }
        }
        dev.off()
      }
      results <- list(severity = severity,
                      shape = shape,
                      statistics = stats)
      class(results) <- "plm_disease"
      invisible(results)
    }

  if(missing(pattern)){
    help_count(img, img_healthy, img_symptoms, img_background, resize, fill_hull, invert,
               index_lb, index_dh, has_white_bg, lesion_size, tolerance, extension, randomize,
               nsample, plot, show_original, show_background, col_leaf,
               col_lesions, col_background,  save_image, dir_original, dir_processed,
               marker, marker_col, marker_size)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("'", pattern, "' pattern not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      init_time <- Sys.time()
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.3), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_count"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_count(x,
                               img_healthy, img_symptoms, img_background, resize, fill_hull, invert,
                               index_lb, index_dh, has_white_bg, lesion_size, tolerance, extension, randomize,
                               nsample, plot, show_original, show_background, col_leaf,
                               col_lesions, col_background,  save_image, dir_original, dir_processed,
                               marker, marker_col, marker_size)
                  })

    } else{
      init_time <- Sys.time()
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        if(verbose == TRUE){
          run_progress(pb, actual = i,
                       text = paste("Processing image", names_plant[i]))
        }
        results[[i]] <-
          help_count(img  = names_plant[i],
                     img_healthy, img_symptoms, img_background, resize, fill_hull, invert,
                     index_lb, index_dh, has_white_bg, lesion_size, tolerance, extension, randomize,
                     nsample, plot, show_original, show_background, col_leaf,
                     col_lesions, col_background,  save_image, dir_original,
                     dir_processed, marker, marker_col, marker_size)
      }
    }
    names(results) <- names_plant
    if(isTRUE(show_features)){
      stats <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["statistics"]],
                            img =  names(results[i]))[,c(3, 1, 2)]
                })
        )
      shape <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["shape"]],
                            img =  names(results[i]))
                })
        )
      if("img" %in% colnames(shape)){
        shape <- shape[, c(ncol(shape), 1:ncol(shape) - 1)]
      }
    } else{
      shape <- NULL
      stats <- NULL
    }
    severity <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["severity"]],
                          img =  names(results[i]))[, c(3, 1:2)]
              })
      )
    message("Done!")
    message("Elapsed time: ", sec_to_hms(as.numeric(difftime(Sys.time(),  init_time, units = "secs"))))
    invisible(
      structure(
        list(severity = severity,
             shape = shape,
             stats = stats,
             parms = list(
               pattern = pattern,
               img_healthy = img_healthy,
               img_symptoms = img_symptoms,
               img_background = img_background,
               dir_original = diretorio_original,
               dir_processed = diretorio_processada,
               save_image = save_image)),
        class = "plm_disease"
      )
    )
  }
}

#' @name measure_disease
#' @export
measure_disease_iter <- function(img,
                                 has_background = TRUE,
                                 r = 2,
                                 viewer = get_pliman_viewer(),
                                 show = "rgb",
                                 index = "NGRDI",
                                 ...){
  viewopt <- c("base", "mapview")
  viewopt <- viewopt[pmatch(viewer[[1]], viewopt)]
  if(viewopt == "base"){
    plot(img)
  }
  if(isTRUE(has_background)){
    # Call the functions independently
    # Call pick_background function
    if(viewopt == "base"){
      message("Use the first mouse button to pick up BACKGROUND colors. Press Est to exit")
    }
    back <- pick_palette(img,
                         r = r,
                         verbose = FALSE,
                         palette = FALSE,
                         plot = FALSE,
                         viewer = viewopt,
                         show = show,
                         index = index,
                         title = "Use the first mouse button to pick up BACKGROUND colors. Click 'Done' to finish",
                         col = "blue")
    if(viewopt != "base"){
      image_view(img[1:10,1:10,])
    }
  } else{
    back <- NULL
  }
  # Call pick_leaf function
  if(viewopt == "base"){
    message("Use the first mouse button to pick up LEAF colors. Press Est to exit")
  }
  leaf <- pick_palette(img,
                       r = r,
                       verbose = FALSE,
                       palette = FALSE,
                       plot = FALSE,
                       viewer = viewopt,
                       show = show,
                       index = index,
                       title = "Use the first mouse button to pick up LEAF colors. Click 'Done' to finish",
                       col = "black")
  if(viewopt != "base"){
    image_view(img[1:10,1:10,])
  }

  # Call pick_disease function
  if(viewopt == "base"){
    message("Use the first mouse button to pick up DISEASE colors. Press Est to exit")
  }
  disease <- pick_palette(img,
                          r = r,
                          verbose = FALSE,
                          palette = FALSE,
                          plot = FALSE,
                          viewer = viewopt,
                          show = show,
                          index = index,
                          title = "Use the first mouse button to pick up DISEASE colors. Click 'Done' to finish",
                          col = "red")

  temp <-
    measure_disease(img = img,
                    img_healthy = leaf,
                    img_symptoms = disease,
                    img_background = back,
                    ...)

  invisible(list(results = temp,
                 leaf = leaf,
                 disease = disease,
                 background = back))
}
