.onAttach <- function(libname, pkgname) {
  vers <-  "v0.1.0"
  packageStartupMessage("|========================================================|")
  packageStartupMessage("| Tools for Plant Image Analysis (pliman) ", vers, "         |")
  packageStartupMessage("| Author: Tiago Olivoto                                  |")
  packageStartupMessage("| Type 'vignette('pliman_start')' for a short tutorial   |")
  packageStartupMessage("| Visit 'https://bit.ly/2RP2kyE' for a complete tutorial |")
  packageStartupMessage("|========================================================|")
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Contorno", "display", "CODE", "dir_original" ,"dir_processada",
      "Spectrum", "value", "area", "id", ".", 'object'))
  }
