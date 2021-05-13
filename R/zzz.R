.onAttach <- function(libname, pkgname) {
  vers <-  "0.2.0"
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
      "Spectrum", "value", "area", "id", ".", 'object'))
  }
