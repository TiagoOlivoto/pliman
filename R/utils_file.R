#' Utilities for file manipulation
#'
#'* `file_extension()` Get the extension of a file.
#'* `file_name()` Get the name of a file.
#'* `file_dir()` Get or directory of a file
#'* `manipulate_files()` Manipulate files in a directory with options to rename
#'(insert prefix or suffix) and save the new files to the same or other provided
#'directory.
#' @name utils_file
#' @param file The file name.
#' @param pattern A file name pattern.
#' @param dir The working directory containig the files to be manipulated.
#'   Defaults to the current working directory.
#' @param prefix,suffix A prefix or suffix to be added in the new file names.
#'   Defaults to `""`.
#' @param sep An optional separator. Defaults to `""`.
#' @param save_to The directory to save the new files. Defaults to the current
#'   working directory. If the file name of a file is not changed, nothing will
#'   occur. If `save_to` refers to a subfolder in the current working directory,
#'   the files will be saved to the given folder. In case of the folder doesn't
#'   exist, it will be created. By default, the files will not be overwritten.
#'   Set `overwrite = TRUE` to overwrite the files.
#' @param overwrite Overwrite the files? Defaults to `FALSE`.
#' @export
#' @examples
#' \donttest{
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
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
}
#' @export
#' @name utils_file
file_name <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[1])
}
#' @export
#' @name utils_file
file_dir <- function(file){
  ex <-  ifelse(grepl(".", file, fixed = TRUE),
                paste0(gsub('.$', "", gsub(basename(file), "", file))),
                file)
  fd <- ifelse(nchar(ex) == 0,
               paste0("."),
               paste0("./", ex))
  fd <- ifelse(grepl(":", fd, fixed = TRUE),
               substring(fd, 3, nchar(fd)),
               fd)
  return(fd)
}
#' @export
#' @name utils_file
manipulate_files <- function(pattern,
                             dir = NULL,
                             prefix = "",
                             suffix = "",
                             sep = "",
                             save_to = NULL,
                             overwrite = FALSE){
  if(is.null(dir)){
    dir <- paste0("./")
  } else{
    dir <- ifelse(grepl(":", dir, fixed = TRUE),
                  file_dir(dir),
                  paste0("./", dir))
  }
  if(is.null(save_to)){
    save_to <- paste0("./")
  } else{
    save_to <- ifelse(grepl(":", save_to, fixed = TRUE),
                      file_dir(save_to),
                      paste0("./", save_to, "/"))
  }
  if(dir.exists(save_to) == FALSE){
    dir.create(save_to)
  }
  if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
    pattern <- "^[0-9].*$"
  }
  old_files <- list.files(dir, pattern = pattern)
  old_files <- paste0(dir, old_files)
  names <- sapply(old_files, file_name)
  extens <- sapply(old_files, file_extension)
  new_files <- paste0(save_to, prefix, sep, names, sep, suffix, ".", extens)
  file.copy(from = old_files, to = new_files, overwrite = overwrite)
}
