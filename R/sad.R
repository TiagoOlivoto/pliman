#' Produces Santandard Area Diagrams
#'
#' Given an object computed with [measure_disease()] or [measure_disease_byl()]
#' a Standard Area Diagram (SAD) with `n` images are returned with the
#' respective severity values.
#'
#' @details
#' The leaves with the smallest and highest severity will always be in the SAD.
#' If `n = 1`, the leaf with the smallest severity will be returned. The others
#' are sampled sequentially to achieve the `n` images after severity has been
#' ordered in an ascending order. For example, if there are 30 leaves and n is
#' set to 3, the leaves sampled will be the 1st, 15th, and 30th with the
#' smallest severity values.
#'
#' The SAD can be only computed if an image pattern name is used in argument
#' `pattern` of [measure_disease()]. If the images are saved, the `n` images
#' will be retrevied from `dir_processed` directory. Otherwise, the severity
#' will be computed again to generate the images.
#'
#' @param object An object computed with [measure_disease()] or
#'   [measure_disease_byl()].
#' @param n The number of leaves in the Standard Area Diagram.
#' @param nrow,ncol The number of rows and columns in the plot. See
#'   [image_combine())]
#' @param show_original Show original images? Defaults to `FALSE`, i.e., a mask
#'   is returned.
#' @param show_contour Show original images? Defaults to `FALSE`, i.e., a mask
#'   is returned.
#' @param ... Other arguments passed on to [measure_disease()].
#' @references Del Ponte EM, Pethybridge SJ, Bock CH, et al (2017) Standard area
#'   diagrams for aiding severity estimation: Scientometrics, pathosystems, and
#'   methodological trends in the last 25 years. Phytopathology 107:1161–1174.
#'   \doi{10.1094/PHYTO-02-17-0069-FI}
#' @return
#' A data frame with the severity values for the `n` sampled leaves. A plot with
#' the standard area diagram can be saved by wrapping [sad()] with [png()].
#' @export
#'
#' @examples
#' \dontrun{
#' library(pliman)
#' sev <-
#' measure_disease(pattern = "sev_leaf",
#'                 img_healthy = "sev_healthy",
#'                 img_symptoms = "sev_sympt",
#'                 img_background  = "sev_back",
#'                 plot = FALSE,
#'                 save_image = TRUE,
#'                 show_original = FALSE,
#'                 dir_original = image_pliman(),
#'                 dir_processed = tempdir())
#' sad(sev, n = 2)
#' }
sad <- function(object,
                n,
                show_original = FALSE,
                show_contour = FALSE,
                nrow = NULL,
                ncol = NULL,
                ...){
  if(inherits(object, "plm_disease")){
    patt <- object$parms$pattern
    h <- object$parms$img_healthy
    s <- object$parms$img_symptoms
    b <- object$parms$img_background
    dir <- object$parms$dir_original
    dir_proc <- object$parms$dir_processed
    save_image <- object$parms$save_image
    nsamples <- n
    measures <-
      transform(object$severity,
                rank = rank(symptomatic))
    n <- nrow(measures)
    seq <- trunc(seq(1, n, length.out = nsamples))
    seq[c(1, length(seq))] <- c(1, n)
    leaves <- measures[which(measures$rank %in% seq),]
    leaves <- leaves[order(leaves$rank),]
    leaves_name <- paste0("proc_", leaves$img, ".jpg")
    if(isFALSE(save_image)){
      if(is.null(patt) | is.null(h) | is.null(s)){
        stop("'pattern', 'img_healthy', and 'img_symptoms' are mandatory arguments.")
      }
      td <- tempdir()
      temp <-
        measure_disease(pattern = patt,
                        img_healthy = h,
                        img_symptoms = s,
                        img_background = b,
                        dir_original = dir,
                        dir_processed = td,
                        show_original = show_original,
                        show_contour = FALSE,
                        save_image = TRUE,
                        parallel = TRUE,
                        verbose = FALSE,
                        ...)
      sads <- image_import(leaves_name, path = td)
      image_combine(sads,
                    labels = paste0(round(leaves$symptomatic, 1), "%"),
                    ncol = ncol,
                    nrow = nrow)
    } else{
      sads <- image_import(leaves_name, path = dir_proc)
      image_combine(sads,
                    labels = paste0(round(leaves$symptomatic, 1), "%"),
                    ncol = ncol,
                    nrow = nrow)
    }
    invisible(leaves)
  }
  if(inherits(object, "plm_disease_byl")){
    patt <- object$parms$pattern
    h <- object$parms$img_healthy
    s <- object$parms$img_symptoms
    b <- object$parms$img_background
    dir <- object$parms$dir_original
    dir_proc <- object$parms$dir_processed
    save_image <- object$parms$save_image
    nsamples <- n
    measures <-
      transform(object$severity,
                rank = rank(symptomatic),
                img = paste0(img, "_", leaf, ".jpg"))
    n <- nrow(measures)
    seq <- trunc(seq(1, n, length.out = nsamples))
    seq[c(1, length(seq))] <- c(1, n)
    leaves <- measures[which(measures$rank %in% seq),]
    leaves <- leaves[order(leaves$rank),]
    leaves_name <- leaves$img
    if(isFALSE(save_image)){
      stop("Standard Area Diagram can only be generated using `save_image = TRUE` set in `measure_disease_byl()`.")
    } else{
      sads <- image_import(leaves_name, path = dir_proc)
      image_combine(sads,
                    labels = paste0(round(leaves$symptomatic, 1), "%"),
                    ncol = ncol,
                    nrow = nrow)
    }
    invisible(leaves)
  }

}
