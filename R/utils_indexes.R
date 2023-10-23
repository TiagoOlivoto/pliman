#' Utilities for image indexes
#'
#' - `pliman_indexes()`: Get all the available indexes in pliman.
#' - `pliman_indexes_rgb()`: Get all the RGB-based available indexes in pliman.
#' - `pliman_indexes_me()`: Get all the multispectral available indexes in pliman.
#' - `pliman_indexes_eq()`: Get the equations of the available indexes.

#' @name utils_indexes
#' @export
pliman_indexes <- function(){
  read.csv(file = system.file("indexes.csv",
                              package = "pliman",
                              mustWork = TRUE),
           header = T, sep = ";")$Index
}
#' @name utils_indexes
#' @export
pliman_indexes_eq <- function(){
  read.csv(file = system.file("indexes.csv",
                              package = "pliman",
                              mustWork = TRUE),
           header = T, sep = ";")
}

#' @name utils_indexes
#' @export
pliman_indexes_rgb <- function(){
  subset(read.csv(file = system.file("indexes.csv",
                                     package = "pliman",
                                     mustWork = TRUE),
                  header = T, sep = ";"),
         Band == "RGB")$Index
}
#' @name utils_indexes
#' @export
pliman_indexes_me <- function(){
  subset(read.csv(file = system.file("indexes.csv",
                                     package = "pliman",
                                     mustWork = TRUE),
                  header = T, sep = ";"),
         Band == "MULTI")$Index
}
