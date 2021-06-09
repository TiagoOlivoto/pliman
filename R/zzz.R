#' @title Sample images
#' @description Sample images installed with the \pkg{pliman} package
#' @format `*.jpg` format
#' *  `la_back.jpg` A cyan palette representing the background of images
#' la_pattern, la_leaves, and soybean_touch.
#' * `la_leaf.jpg` A sample of the leaves in `la_leaves`
#' * `la_leaves.JPG` Tree leaves with a sample of known area.
#' * `la_pattern.JPG` Tree leaves with a yellow sample of known area.
#' * `la_temp.jpg` The yellow sample of `la_pattern.JPG`.
#' * `objects_300dpi.jpg` An image with 300 dpi resolution.
#' * `sev_leaf.jpg` A soybean leaf with a blue background.
#' * `sev_leaf_nb.jpg` A soybean leaf without background.
#' * `sev_back.jpg` A blue palette representing the background of `sev_leaf`.
#' * `sev_healthy.jpg` Healthy area of `sev_leaf`.
#' * `sev_sympt.jpg` The symptomatic area `sev_leaf`.
#' * `soy_green.jpg` Soybean grains with a white background.
#' * `soybean_grain.jpg` A sample palette of the grains in `soy_green`.
#' * `soybean_touch.jpg` Soybean grains with a cyan background touching one each other.
#' @md
#' @source Personal data
#' @name pliman_images
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords images
NULL

.onAttach <- function(libname, pkgname) {
  vers <-  "0.3.0"
  packageStartupMessage("|========================================================|")
  packageStartupMessage("| Tools for Plant Image Analysis (pliman ", vers,  ")          |")
  packageStartupMessage("| Author: Tiago Olivoto                                  |")
  packageStartupMessage("| Type 'vignette('pliman_start')' for a short tutorial   |")
  packageStartupMessage("| Visit 'https://bit.ly/3eL0dF3' for a complete tutorial |")
  packageStartupMessage("|========================================================|")
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Contorno", "display", "CODE", "dir_original" ,"dir_processada",
      "Spectrum", "value", "area", "id", ".", "object", "s.radius.max",
      "s.radius.min"))
  }
