#' Forward-pipe operator
#'
#' Pipe an object forward into a function or call expression.
#'
#' @param lhs The result you are piping.
#' @param rhs Where you are piping the result to.
#'
#' @examples
#' library(pliman)
#'
#' # Basic use:
#'  iris %>% head()
#'
#' # use to apply several functions to an image
#' img <- image_pliman("la_leaves.jpg")
#'
#' img %>%
#'  image_resize(50) %>%        # resize to 50% of the original size
#'  object_isolate(id = 1) %>%  # isolate object 1
#'  image_filter() %>%          # apply a median filter
#'  plot()                      # plot
#'
#' @name pipe
#'
#' @author
#' Nathan Eastwood \email{nathan.eastwood@icloud.com} and Antoine Fabri
#' \email{antoine.fabri@@gmail.com}. The code was obtained from poorman package
#' at \url{https://github.com/nathaneastwood/poorman/blob/master/R/pipe.R}
#'
#' @export
`%>%` <- function(lhs, rhs) {
  rhs_call <- insert_dot(substitute(rhs))
  eval(rhs_call, envir = list(`.` = lhs), enclos = parent.frame())
}

#' @author Antoine Fabri
#' @noRd
insert_dot <- function(expr) {
  if (is.symbol(expr) || expr[[1]] == quote(`(`)) {
    # if a symbol or an expression inside parentheses, make it a call with dot
    # arg
    expr <- as.call(c(expr, quote(`.`)))
  } else if (length(expr) == 1) {
    # if a call without an arg, give it a dot arg
    expr <- as.call(c(expr[[1]], quote(`.`)))
  } else if (
    expr[[1]] != quote(`{`) &&
    !any(vapply(expr[-1], identical, quote(`.`), FUN.VALUE = logical(1))) &&
    !any(vapply(expr[-1], identical, quote(`!!!.`), FUN.VALUE = logical(1)))
  ) {
    # if a call with args but no dot in arg, insert one first
    expr <- as.call(c(expr[[1]], quote(`.`), as.list(expr[-1])))
  }
  expr
}
