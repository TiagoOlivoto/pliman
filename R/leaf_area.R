#'Calculates the leaf area
#'
#'* `leaf_area ()` Calculates the leaf area using an image with leaves and a
#'template with a known area. A general linear model (binomial family) fitted to
#'the RGB values is used to first separate the leaves and template from the
#'background and then the leaves from the template. The leaf area is then
#'calculated for each leaf based on the pixel area. By using `img_pattern` it is
#'possible to process several images with common pattern names that are stored
#'in the current working directory or in the subdirectory informed in
#'`dir_originals`.
#' @name leaf_area
#' @param img The image to be analyzed.
#' @param img_pattern A pattern of file name used to identify images to be
#'   processed. For example, if `img_pattern = "im"` all images that the name
#'   matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
#'   Providing any number as pattern (e.g., `img_pattern = "1"`) will select
#'   images that are named as 1.-, 2.-, and so on.
#' @param img_leaf A color palette of the leaves.
#' @param img_background A color palette of background area.
#' @param img_template A color palette of the template areas.
#' @param area_template The known area of the template. The leaf area will be
#'   given in the same unit as `area_template`.
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time, especially when `img_pattern` is used is informed. The
#'   number of sections is set up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param lower_size Lower limit for size for the image analysis. Leaf images
#'   often contain dirt and dust. To prevent dust from affecting the image
#'   analysis, the lower limit of analyzed size is set to 0.1, i.e., objects
#'   with lesser than 10% of the mean of all objects are removed. One can set a
#'   known area or use `lower_limit = 0` to select all objects (not advised).
#' @param upper_size Upper limit for size for the image analysis. Defaults to
#'   `NULL`, i.e., no upper limit used.
#' @param randomize Randomize the lines before training the model?
#' @param nrows The number of lines to be used in training step.
#' @param show_image Show image after processing?
#' @param show_original Show the symptoms in the original image?
#' @param show_background Show the background? Defaults to `TRUE`. A white
#'   background is shown by default when `show_original = FALSE`.
#' @param col_background Background color after image processing.
#' @param col_leaf Leaf color after image processing when `show_original =
#'   FALSE`. Defaults to `"green"`.
#' @param text_col,text_size,text_digits The color, size and significant digits
#'   used in the text. The shows the pattern `o|a`, where `o` and `a` are the
#'   object id and its area, respectively.
#' @param save_image Save the image after processing? The image is saved in the
#'   current working directory named as `proc_*` where `*` is the image name
#'   given in `img`.
#' @param dir_original,dir_processed The directory containing the original and
#'   processed images. Defaults to `NULL`. In this case, the function will
#'   search for the image `img` in the current working directory. After
#'   processing, when `save_image = TRUE`, the processed image will be also
#'   saved in such a directory.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @return A data frame with the results for each image.
#' @export
#' @md
#' @importFrom graphics text par points
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_import(image_pliman("la_pattern.JPG"))
#' leaf <- image_import(image_pliman("la_leaf.jpg"))
#' tmpl <- image_import(image_pliman("la_temp.jpg"))
#' background <- image_import(image_pliman("la_back.jpg"))
#'
#' # Combine the images
#' image_combine(img, leaf, tmpl, background)
#'
#' # Computes the leaf area
#' area <-
#' leaf_area(img = img,
#'           img_leaf = leaf,
#'           img_template = tmpl,
#'           img_background = background,
#'           area_template = 4,
#'           text_col = "white")
#' get_measures(area)
#' }
#'
leaf_area <- function(img,
                      img_leaf,
                      img_background,
                      img_template,
                      area_template,
                      parallel = FALSE,
                      workers = NULL,
                      img_pattern = NULL,
                      lower_size = NULL,
                      upper_size = NULL,
                      randomize = TRUE,
                      nrows = 10000,
                      show_image = TRUE,
                      show_original = TRUE,
                      show_background = TRUE,
                      col_background = NULL,
                      col_leaf = "green",
                      text_col = "black",
                      text_size = 1,
                      text_digits = 2,
                      save_image = FALSE,
                      dir_original = NULL,
                      dir_processed = NULL,
                      verbose = TRUE){
  # Some parts adapted from
  # https://github.com/AlcineiAzevedo/Segmentacao-conchonilha2
  # Thanks to Alcinei Azevedo for his tips
  if(is.null(dir_original)){
    diretorio_original <- paste("./", sep = "")
  } else{
    diretorio_original <- paste("./", dir_original, sep = "")
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste("./", sep = "")
  } else{
    diretorio_processada <- paste("./", dir_processed, sep = "")
  }
  help_area <-
    function(img, img_leaf, img_background, img_template, area_template,
             lower_size, upper_size , randomize, nrows, show_image,
             show_original, show_background, col_background, text_col,
             text_size, text_digits, save_image, dir_original, dir_processed){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img)
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(diretorio_original, "/", name_ori, ".", extens_ori, sep = ""))
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "jpg"
      }
      if(is.character(img_leaf)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img_leaf, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img_leaf)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_leaf <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      if(is.character(img_template)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img_template, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img_template)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_template <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      if(is.character(img_background)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img_background, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = img_background)
        name <- file_name(imag)
        extens <- file_extension(imag)
        img_background <- image_import(paste(diretorio_original, "/", name, ".", extens, sep = ""))
      }
      original <- image_to_mat(img)
      leaf <- image_to_mat(img_leaf)
      template <- image_to_mat(img_template)
      background <- image_to_mat(img_background)
      # separate image from background
      background_resto <-
        transform(rbind(leaf$df_in[sample(1:nrow(leaf$df_in)),][1:nrows,],
                        template$df_in[sample(1:nrow(template$df_in)),][1:nrows,],
                        background$df_in[sample(1:nrow(background$df_in)),][1:nrows,]),
                  Y = ifelse(CODE == "img_background", 1, 0))
      background_resto$CODE <- NULL
      modelo1 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = background_resto))
      pred1 <- round(predict(modelo1, newdata = original$df_in, type="response"), 0)
      plant_background <- matrix(pred1, ncol = ncol(original$R))
      plant_background <- image_correct(plant_background, perc = 0.009)
      plant_background[plant_background == 1] <- 2
      # image_show(plant_background!=2)
      # separate leaf from template
      leaf_template <-
        transform(rbind(leaf$df_in[sample(1:nrow(leaf$df_in)),][1:nrows,],
                        template$df_in[sample(1:nrow(template$df_in)),][1:nrows,]),
                  Y = ifelse(CODE == "img_leaf", 0, 1))
      background_resto$CODE <- NULL
      modelo2 <- suppressWarnings(glm(Y ~ R + G + B, family = binomial("logit"), data = leaf_template))
      # isolate plant
      ID <- c(plant_background == 0)
      pred2 <- round(predict(modelo2, newdata = original$df_in[ID,], type="response"), 0)
      pred3 <- round(predict(modelo2, newdata = original$df_in, type="response"), 0)
      leaf_template <- matrix(pred3, ncol = ncol(original$R))
      leaf_template <- image_correct(leaf_template, perc = 0.009)
      plant_background[leaf_template == 1] <- 3
      mpred1 <- bwlabel(leaf_template == 1)
      shape_template <-
        cbind(data.frame(computeFeatures.shape(mpred1)),
              data.frame(computeFeatures.moment(mpred1))[,1:2]
        )
      shape_template$area <- area_template
      shape_template <- shape_template[shape_template$s.area >= mean(shape_template$s.area), ]
      npix_ref <- shape_template[1, 1]
      mpred2 <- bwlabel(plant_background == 0)
      shape_leaf <-
        cbind(data.frame(computeFeatures.shape(mpred2)),
              data.frame(computeFeatures.moment(mpred2))[, c(1, 2)]
        )
      shape_leaf$area <- shape_leaf$s.area * area_template / npix_ref
      if(!is.null(lower_size)){
        shape_leaf <- shape_leaf[shape_leaf$area > lower_size, ]
      } else{
        shape_leaf <- shape_leaf[shape_leaf$area > mean(shape_leaf$area) * 0.1, ]
      }
      if(!is.null(upper_size)){
        shape_leaf <- shape_leaf[shape_leaf$area < upper_size, ]
        shape_template <- shape_template[shape_template$area < upper_size, ]
      }
      shape <- rbind(shape_leaf, shape_template)
      shape$id <- 1:nrow(shape)
      shape <-
        transform(shape[, c(10, 7, 8, 1, 9, 2:6)],
                  label = paste(id, "|", round(area, text_digits), sep = ""))
      if(show_original == TRUE){
        im2 <- img
        if(!is.null(col_background)){
          col_background <- col2rgb(col_background)
          im2@.Data[,,1][!ID] <- col_background[1]
          im2@.Data[,,2][!ID] <- col_background[2]
          im2@.Data[,,3][!ID] <- col_background[3]
        }
      } else{
        if(is.null(col_background)){
          col_background <- col2rgb("white")
        } else{
          col_background <- col2rgb(col_background)
        }
        col_leaf <- col2rgb(col_leaf)
        im2 <- img
        im2@.Data[,,1][ID] <- col_leaf[1]
        im2@.Data[,,2][ID] <- col_leaf[2]
        im2@.Data[,,3][ID] <- col_leaf[3]
        im2@.Data[,,1][!ID] <- col_background[1]
        im2@.Data[,,2][!ID] <- col_background[2]
        im2@.Data[,,3][!ID] <- col_background[3]
      }
      if(show_image == TRUE){
        image_show(im2)
        text(shape[,2],
             shape[,3],
             shape$label,
             col = text_col,
             cex = text_size)
      }
      if(save_image == TRUE){
        if(dir.exists(diretorio_processada) == FALSE){
          dir.create(diretorio_processada)
        }
        png(paste(diretorio_processada,
                  "/proc_", name_ori, ".",
                  extens_ori,
                  collapse = "", sep = ""),
            width = dim(im2@.Data)[1],
            height = dim(im2@.Data)[2])
        image_show(im2)
        text(shape[,2],
             shape[,3],
             shape$label,
             col = text_col,
             cex = text_size)
        dev.off()
      }
      shape <- shape[,c(1:3, 5:7, 9:10, 8, 11)]
      colnames(shape) <- c("id", "x", "y", "area", "perimeter", "radius_mean",
                           "radius_min", "radius_max", "radius_sd", "label")
      id_obj <- shape[which(shape$area == area_template),1]
      class(shape) <- c("data.frame", "plm_la", id_obj)
      return(shape)
    }
  if(missing(img_pattern)){
    help_area(img, img_leaf, img_background, img_template, area_template,
              lower_size, upper_size , randomize, nrows, show_image,
              show_original, show_background, col_background, text_col,
              text_size, text_digits, save_image, dir_original, dir_processed)
  } else{
    if(img_pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      img_pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = img_pattern, diretorio_original)
    extensions <-  as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(img_pattern, names_plant)) == 0){
      stop(paste("'", img_pattern, "' pattern not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust,
                    varlist = c("names_plant", "help_area", "file_name",
                                "check_names_dir", "file_extension", "image_import",
                                "image_binary", "watershed", "distmap", "computeFeatures.moment",
                                "computeFeatures.shape", "colorLabels", "image_show",
                                "image_to_mat", "image_correct", "bwlabel"),
                    envir=environment())
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      results <-
        parLapply(clust, names_plant,
                  function(x){
                    help_area(x,
                              img_leaf, img_background, img_template, area_template,
                              lower_size, upper_size , randomize, nrows, show_image,
                              show_original, show_background, col_background, text_col,
                              text_size, text_digits, save_image, dir_original, dir_processed)
                  })

    } else{
      results <- list()
      pb <- progress(max = length(plants), style = 4)
      for (i in 1:length(plants)) {
        run_progress(pb, actual = i, text = paste("Processing image", names_plant[i]))
        results[[i]] <-
          help_area(img  = names_plant[i],
                    img_leaf, img_background, img_template, area_template,
                    lower_size, upper_size , randomize, nrows, show_image,
                    show_original, show_background, col_background, text_col,
                    text_size, text_digits, save_image, dir_original, dir_processed)
      }
    }
    return(do.call(rbind, results))
  }
}
