.onAttach <- function(libname, pkgname) {
  vers <-  "0.0.0.9000"
  packageStartupMessage("|===========================================|")
  packageStartupMessage("| Tools for Plant Image Analysis (pliman)   |")
  packageStartupMessage("| Author: Tiago Olivoto                     |")
  packageStartupMessage("|===========================================|")
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Contorno", "display", "CODE", "dir_original" ,"dir_processada",
      "Spectrum", "value", "area", "id", "."))
  }
