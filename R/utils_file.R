#' Utilities for file manipulation
#'
#'* `file_extension()` Get the extension of a file.
#'* `file_name()` Get the name of a file.
#'* `file_dir()` Get or directory of a file
#'* `manipulate_files()` Manipulate files in a directory with options to rename
#'(insert prefix or suffix) and save the new files to the same or other provided
#'directory.
#'* `pliman_indexes()` Get the indexes available in pliman.
#'* `pliman_indexes_eq()` Get the equation of the indexes available
#'in pliman.
#' @name utils_file
#' @param file The file name.
#' @param pattern A file name pattern.
#' @param dir The working directory containing the files to be manipulated.
#'   Defaults to the current working directory.
#' @param prefix,suffix A prefix or suffix to be added in the new file names.
#'   Defaults to `NULL` (no prefix or suffix).
#' @param name The name of the new files. Defaults to `NULL` (original names).
#'   `name` can be either a single value or a character vector of the same
#'   length as the number of files manipulated. If one value is informed, a
#'   sequential vector of names will be created as "`name`_1", "`name`_2", and
#'   so on.
#' @param extension The new extension of the file. If not declared (default),
#'   the original extensions will be used.
#' @param sep An optional separator. Defaults to `""`.
#' @param save_to The directory to save the new files. Defaults to the current
#'   working directory. If the file name of a file is not changed, nothing will
#'   occur. If `save_to` refers to a subfolder in the current working directory,
#'   the files will be saved to the given folder. In case of the folder doesn't
#'   exist, it will be created. By default, the files will not be overwritten.
#'   Set `overwrite = TRUE` to overwrite the files.
#' @param overwrite Overwrite the files? Defaults to `FALSE`.
#' @param remove_original Remove original files after manipulation? defaults to
#'   `FALSE`. If `TRUE` the files in `pattern` will be removed.
#' @param verbose If `FALSE`, the code is run silently.
#' @return
#' * `file_extension()`,  `file_name()`, and `file_dir()` return a character
#' string.
#' * `manipulate_files()` No return value. If `verbose == TRUE`, a message is
#' printed indicating which operation succeeded (or not) for each of the files
#' attempted.
#' @export
#' @examples
#' \donttest{
#' library(pliman)
#' # get file name, directory and extension
#' file <- "E:/my_folder/my_subfolder/image1.png"
#' file_dir(file)
#' file_name(file)
#' file_extension(file)
#'
#' # manipulate files
#' dir <- tempdir()
#' list.files(dir)
#' file.create(paste0(dir, "/test.txt"))
#' list.files(dir)
#' manipulate_files("test",
#'                  dir = paste0(dir, "\\"),
#'                 prefix = "chang_",
#'                 save_to = paste0(dir, "\\"),
#'                 overwrite = TRUE)
#' list.files(dir)
#' }
file_extension <- function(file){
  sapply(seq_along(file), function(x){
    strsplit(basename(file[x]), split="\\.")[[1]][-1]
  })
}
#' @export
#' @name utils_file
file_name <- function(file){
  sapply(seq_along(file), function(x){
    strsplit(basename(file[x]), split="\\.")[[1]][1]
  })
}
#' @export
#' @name utils_file
file_dir <- function(file){
  ex <-  ifelse(grepl(".", file, fixed = TRUE),
                sapply(seq_along(file),
                       function(x){
                         paste0(gsub('.$', "", gsub(basename(file[x]), "", file[x])))
                       }),
                file)
  fd <-
    sapply(seq_along(ex),
           function(x){
             ifelse(nchar(x) == 0,
                    paste0("."),
                    ifelse(grepl("^[/]", file[x]), ex[x], paste0("./", ex[x])))
           })
  fd <-
    sapply(seq_along(fd),
           function(x){
             ifelse(grepl(":", fd[x], fixed = TRUE),
                    substring(fd[x], 3, nchar(fd[x])),
                    fd[x])
           })
  return(fd)
}
#' @export
#' @name utils_file
manipulate_files <- function(pattern,
                             dir = NULL,
                             prefix = NULL,
                             name = NULL,
                             suffix = NULL,
                             extension = NULL,
                             sep = "",
                             save_to = NULL,
                             overwrite = FALSE,
                             remove_original = FALSE,
                             verbose = TRUE){
  check_dir <- is.null(dir)
  if(check_dir){
    dir <- paste0("./")
  } else{
    dir <- ifelse(grepl(":", dir, fixed = TRUE),
                  file_dir(dir),
                  paste0("./", dir))
  }
  if(is.null(save_to)){
    save_to <- paste0(ifelse(is.null(save_to), paste0(dir),  paste0(dir, "/")))
    save_to <- ifelse(check_dir, save_to, paste0(save_to, "/"))
  } else{
    save_to <- ifelse(grepl(":", save_to, fixed = TRUE),
                      file_dir(save_to),
                      paste0("./", save_to, "/"))
  }
  if(dir.exists(save_to) == FALSE){
    dir.create(save_to, recursive = TRUE)
  }
  if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
    pattern <- "^[0-9].*$"
  }
  old_files <- list.files(dir, pattern = "00")
  old_files <- paste0(ifelse(nchar(dir) !=2,
                             paste0(dir, "/"),
                             paste(dir)), old_files)
  names <- sapply(old_files, file_name)

  ifelse(!missing(extension),
         extens <- rep(extension, length(old_files)),
         extens <- sapply(old_files, file_extension))
  prefix <- ifelse(is.null(prefix), "", prefix)
  if(is.null(name)){
    name <- names
  } else{
    if(length(name) == 1){
      name <-
        unlist(lapply(seq_along(names),
                      function(i){
                        paste0(name, i, collapse = "_")
                      }))
    } else{
      name <- name
      if(length(name) != length(names)){
        stop("The length of name must be equal to the number of files (", length(names), ").")
      }
    }
  }
  suffix <- ifelse(is.null(suffix), "", suffix)
  new_files <- paste0(save_to, prefix, sep, name, sep, suffix, ".", extens)
  a <- file.copy(from = old_files, to = new_files, overwrite = overwrite)
  if(remove_original == TRUE){
    invisible(file.remove(old_files))
  }
  if(verbose == TRUE){
    if(remove_original == TRUE){
      message(length(old_files), " files successfully deleted from '", dir, "'")
    }
    if(all(a) == TRUE){
      message(length(a), " files successfully copied to '", save_to, "'")
    }
    if(any(a) == FALSE){
      warning("Failed to copy ", length(which(a == FALSE)), " files.", call. = FALSE)
    }
  }
}

#' @name utils_file
#' @export
pliman_indexes <- function(){
  read.csv(file=system.file("indexes.csv",
                            package = "pliman",
                            mustWork = TRUE),
           header = T, sep = ";")$Index
}
#' @name utils_file
#' @export
pliman_indexes_eq <- function(){
  read.csv(file=system.file("indexes.csv",
                            package = "pliman",
                            mustWork = TRUE),
           header = T, sep = ";")$Equation
}
