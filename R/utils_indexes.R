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

#' List Computable Indexes Based on Available Bands
#'
#' This function reads index equations from a CSV file included in the `pliman`
#' package, determines which bands are used in each index equation, and checks
#' which indexes can be computed based on the provided available bands.
#'
#' @param available A character vector of available bands (e.g., `c("R", "G")`).
#' @return A data frame of indexes that can be computed with the available bands.
#' @examples
#' library(pliman)
#' available_bands <- c("R", "G")
#' computable_indexes <- pliman_indexes_ican_compute(available_bands)
#' print(computable_indexes)
#' @export
pliman_indexes_ican_compute <- function(available){
  # Read the indexes CSV file from the package
  ind <- read.csv(file = system.file("indexes.csv", package = "pliman", mustWork = TRUE),
                  header = TRUE, sep = ";")

  # Regular expression pattern to extract words
  pattern <- "\\b\\w+\\b"

  # Reserved words that are not band names
  reserved <- c("exp", "abs", "min", "max", "median", "sum", "sqrt", "cos", "sin", "tan", "log", "log10", "atan")

  test <- list()

  # Loop through each index equation
  for (i in 1:nrow(ind)) {
    index <- ind$Equation[[i]]

    # Extract all words from the equation
    layersused <- setdiff(unlist(regmatches(index, gregexpr(pattern, index, perl = TRUE))), reserved)

    # Filter out words that are not purely alphabetic (i.e., potential band names)
    onlychar <- suppressWarnings(is.na(as.numeric(layersused)))
    test[[i]] <- layersused[onlychar]
  }

  # Set names of the test list to the index names
  names(test) <- ind$Index

  # Function to check if all elements of a vector are in the available bands
  check_bands <- function(bands, available) {
    all(bands %in% available)
  }

  # Apply the check_bands function to each element in the list
  result <- sapply(test, check_bands, available)

  # Return the data frame of computable indexes
  ind[ind$Index %in% ind$Index[result], ]
}

