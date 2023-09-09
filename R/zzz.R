#' @title Sample images
#' @description Sample images installed with the \pkg{pliman} package
#' @format `*.jpg` format
#' * `flax_leaves.jpg` Flax leaves in a white background
#' * `flax_grains.jpg` Flax grains with background light.
#' *  `la_back.jpg` A cyan palette representing the background of images
#' la_pattern, la_leaves, and soybean_touch.
#' * `la_leaf.jpg` A sample of the leaves in `la_leaves`
#' * `la_leaves.jpg` Tree leaves with a sample of known area.
#' * `mult_leaves.jpg` Three soybean leaflets with soybean rust symptoms.
#' * `objects_300dpi.jpg` An image with 300 dpi resolution.
#' * `potato_leaves.jpg` Three potato leaves, which were gathered from Gupta et
#' al. (2020).
#' * `sev_leaf.jpg` A soybean leaf with a blue background.
#' * `sev_leaf_nb.jpg` A soybean leaf without background.
#' * `sev_back.jpg` A blue palette representing the background of `sev_leaf`.
#' * `sev_healthy.jpg` Healthy area of `sev_leaf`.
#' * `sev_sympt.jpg` The symptomatic area `sev_leaf`.
#' * `shadow.jpg` A shaded leaf, useful to test adaptive thresholding
#' * `soy_green.jpg` Soybean grains with a white background.
#' * `soybean_grain.jpg` A sample palette of the grains in `soy_green`.
#' * `soybean_touch.jpg` Soybean grains with a cyan background touching one each
#' other.
#' * `field_mosaic.jpg` An UVA image from a soybean field.
#' @md
#' @source Personal data, Gupta et al. (2020).
#' @references Gupta, S., Rosenthal, D. M., Stinchcombe, J. R., & Baucom, R. S.
#'   (2020). The remarkable morphological diversity of leaf shape in sweet
#'   potato (Ipomoea batatas): the influence of genetics, environment, and G×E.
#'   New Phytologist, 225(5), 2183–2195. \doi{10.1111/NPH.16286}
#' @name pliman_images
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords images
NULL

.onAttach <- function(libname, pkgname) {
  vers <-  "2.0.1"
  packageStartupMessage("|==========================================================|")
  packageStartupMessage("| Tools for Plant Image Analysis (pliman ", vers,  ")            |")
  packageStartupMessage("| Author: Tiago Olivoto                                    |")
  packageStartupMessage("| Type `citation('pliman')` to know how to cite pliman     |")
  packageStartupMessage("| Visit 'http://bit.ly/pkg_pliman' for a complete tutorial |")
  packageStartupMessage("|==========================================================|")
  check_ebi()
  check_mapview()
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Contorno", "display", "CODE", "dir_original" ,"dir_processada",
      "Spectrum", "value", "area", "id", ".", "object", "s.radius.max",
      "s.radius.min", "y", "s.area", "s.perimeter", "symptomatic", "m.eccentricity",
      "m.majoraxis", "s.radius.mean", "n_greater", "n_less", "setNames", "s.radius.sd",
      "perimeter", "radius_max", "radius_mean", "radius_min", "radius_sd", "X1",
      "X2", "%dopar%", "i", "img", "plotn", "x"))
}


#' @title Contour outlines from five leaves
#' @description A list of contour outlines from five leaves. It may be used as
#'   example in some functions such as [efourier()]
#' @format A list with five objects
#' * `leaf_1`
#' * `leaf_2`
#' * `leaf_3`
#' * `leaf_4`
#' * `leaf_5`
#'
#' Each object is a `data.frame` with the coordinates for the outline perimeter
#' @md
#' @name contours
#' @source Personal data. The images were obtained in the Flavia data set
#'   downlodable at <https://flavia.sourceforge.net/>
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' Global option for controlling the viewer in pliman package
#'
#' Users can set the value of this option using `options("pliman_viewer", value)`.
#' The default value is "base". Use "mapview" to allow image to be
#' plotted/edited using the R packages mapview and mapedit
#'
#' @name pliman_viewer
#' @aliases pliman_viewer
#' @rdname pliman_viewer
NULL


#' Get the value of the pliman_viewer option
#'
#' Retrieves the current value of the pliman_viewer option used in the package.
#'
#' @return The current value of the pliman_viewer option.
#' @export
get_pliman_viewer <- function() {
  option_value <- options("pliman_viewer")[[1]]
  if (!is.null(option_value)) {
    option_value
  } else {
    "base"  # Default value
  }
}

#' Set the value of the pliman_viewer option
#'
#' Sets the value of the pliman_viewer option used in the package.
#'
#' @param value The value to be set for the pliman_viewer option.
#' @export
set_pliman_viewer <- function(value) {
  if(!value %in% c("base", "mapview")){
    stop("`value` must be either 'base' or 'mapview'.", call. = FALSE)
  }
  options(pliman_viewer = value)
}

### Utilities for defining viewer options
.onLoad <- function(libname, pkgname) {
  options(pliman_viewer = "base")
}

