#' Image segmentation with pixels proportion
#'
#' Provides (iterative) image segmentation and get the proportion of pixels.
#' This can be used as an alternative way to compute disease severity without
#' using color palettes.
#' @param image An object (or a list of objects) of class `Image`.
#' @param nseg The number of iterative segmentation steps to be performed.
#' @param index The index to be used in the image segmentation. See
#'   [image_index()] for more details.
#' @param fill_hull Fill holes in the objects? Defaults to `FALSE`.
#' @param filter Performs median filtering after image processing? defaults to
#'   `FALSE`. See more at [image_filter()].
#' @param parallel Processes the images asynchronously (in parallel) in separate
#'   R sessions running in the background on the same machine. It may speed up
#'   the processing time when `image` is a list. The number of sections is set
#'   up to 90% of available cores.
#' @param workers A positive numeric scalar or a function specifying the maximum
#'   number of parallel processes that can be active at the same time.
#' @param show_image Show the image results? Defaults to `TRUE`.
#' @param verbose If `TRUE` (default) a summary is shown in the console.
#' @param nrow,ncol Arguments passed on to [image_combine()]. The number of rows
#'   or columns in the plot grid. Defaults to `NULL`, i.e., a square grid is
#'   produced.
#' @return A list with the following objects
#' * `results` A data frame with the number of pixels and proportion of pixels
#' in relation to the previous segmentation.
#' * `images` A list of segmented images.
#' @export
#' @importFrom utils menu
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' img <- image_import(image_pliman("sev_leaf.jpg"))
#' image_show(img)
#'
#' prop_segmented(img,
#'                nseg = 2,
#'                index = c("G", "GLI"),
#'                ncol = 3)
#' }
prop_segmented <- function(image,
                           nseg = 1,
                           index = NULL,
                           fill_hull = FALSE,
                           filter = FALSE,
                           parallel = FALSE,
                           workers = NULL,
                           show_image = TRUE,
                           ncol = NULL,
                           nrow = NULL,
                           verbose = TRUE){
  if(is.list(image)){
    if(!all(sapply(image, class) == "Image")){
      stop("All images must be of class 'Image'")
    }
    if(parallel == TRUE){
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.9), workers)
      clust <- makeCluster(nworkers)
      clusterExport(clust, c("image", "image_segment", "image_combine"))
      on.exit(stopCluster(clust))
      if(verbose == TRUE){
        message("Image processing using multiple sessions (",nworkers, "). Please wait.")
      }
      a <- parLapply(clust, image, prop_segmented, nseg, index, fill_hull, show_image, ncol, nrow,  verbose)
    } else{
      a <- lapply(image, prop_segmented, nseg, index, fill_hull, show_image, ncol, nrow,  verbose)
    }
    results <-
      do.call(rbind, lapply(a, function(x){
        x$results
      }))
    images <-
      lapply(a, function(x){
        x$images
      })
    return(list(results = results,
                images = images))
  } else{
    if(filter == TRUE){
      image <- image_filter(image)
    }
    if(nseg == 1){
      if(is.null(index)){
        image_segment(image,
                      index = "all",
                      fill_hull = fill_hull,
                      show_image = show_image)
        avali_index <- pliman_indexes()
        index <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB",
                 "GR", "BI", "BIM", "SCI", "GLI", "HI", "NGRDI", "NDGBI", "NDRBI", "I",
                 "S", "VARI", "HUE", "HUE2", "BGI", "L")
      } else{
        index <- index
      }
      segmented <-
        image_segment(image,
                      index = index,
                      fill_hull = fill_hull,
                      show_image = FALSE)
      total <- length(image)
      segm <- length(which(segmented[[1]][["image"]] != 1))
      prop <- segm / total * 100
      results <- data.frame(total = total,
                            segmented = segm,
                            prop = prop)
      imgs <- list(image, segmented[[1]][["image"]])
      if(verbose){
        print(results)
      }
      if(show_image == TRUE){
        image_combine(imgs, nrow = nrow, ncol = ncol)
      }
      invisible(list(results = results,
                     images = imgs))
    } else{
      if(is.null(index)){
        image_segment(image,
                      index = "all",
                      fill_hull = fill_hull,
                      show_image = show_image)
        avali_index <- pliman_indexes()
        indx <-
          switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                 "R", "G", "B", "NR", "NG", "NB", "GB", "RB",
                 "GR", "BI", "BIM", "SCI", "GLI", "HI", "NGRDI", "NDGBI", "NDRBI", "I",
                 "S", "VARI", "HUE", "HUE2", "BGI", "L")
      } else{
        if(length(index) != nseg){
          stop("Length of 'index' must be equal 'nseg'.", call. = FALSE)
        }
        indx <- index[1]
      }
      segmented <- list()
      total <- length(image)
      first <-
        image_segment(image,
                      index = indx,
                      fill_hull = fill_hull,
                      show_image = FALSE)
      segmented[[1]] <- first
      for (i in 2:(nseg)) {
        if(is.null(index)){
          image_segment(first,
                        index = "all",
                        fill_hull = fill_hull,
                        show_image = T)
          avali_index <- pliman_indexes()
          indx <-
            switch(menu(avali_index, title = "Choose the index to segment the image, or type 0 to exit"),
                   "R", "G", "B", "NR", "NG", "NB", "GB", "RB",
                   "GR", "BI", "BIM", "SCI", "GLI", "HI", "NGRDI", "NDGBI", "NDRBI", "I",
                   "S", "VARI", "HUE", "HUE2", "BGI", "L")
          if(is.null(indx)){
            break
          }
        } else{
          indx <- index[i]
        }
        second <-
          image_segment(first,
                        index = indx,
                        fill_hull = fill_hull,
                        show_image = FALSE)
        segmented[[i]] <- second
        first <- second
      }
      pixels <-
        rbind(total,
              do.call(rbind,
                      lapply(segmented, function(x){
                        length(which(x[[1]][["image"]] != 1))
                      })
              )
        )
      rownames(pixels) <- NULL
      colnames(pixels) <- "pixels"
      prop <- NULL
      for(i in 2:nrow(pixels)){
        prop[1] <- 100
        prop[i] <- pixels[i] / pixels[i-1] * 100
      }
      pixels <- data.frame(pixels)
      pixels$prop <- prop
      imgs <- lapply(segmented, function(x){
        x[[1]][["image"]]
      })
      imgs <- c(list(image), imgs)
      names <- paste("seg", 1:length(segmented), sep = "")
      names(imgs) <- c("original", names)
      pixels <- transform(pixels, image = c("original",names))
      pixels <- pixels[,c(3, 1, 2)]
      if(verbose){
        print(pixels)
      }
      if(show_image == TRUE){
        image_combine(imgs, ncol = ncol, nrow = nrow)
      }
      invisible(list(results = pixels,
                     images = imgs))
    }
  }
}
