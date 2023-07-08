#' @title Utilities for handling with rows and columns
#' @name utils_rows_cols
#' @description
#'
#' * `columns_to_rownames()`: Move a column of `.data` to its row
#' names.
#' * `rownames_to_column()`: Move the row names of `.data` to a new
#' column.
#' * `remove_rownames()`: Remove the row names of `.data`.
#'
#' * `round_cols()` Rounds the values of all numeric variables to the specified
#' number of decimal places (default 2).
#'
#' @param .data A data frame
#' @param var Name of column to use for rownames.
#' @param digits The number of significant figures. Defaults to `2.`
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(pliman)
#' iris2 <- iris |> rownames_to_column()
#' head(iris2)
#' iris2$rowname <- paste0("r", iris2$rowname)
#' iris2 |> column_to_rownames("rowname") |> head()
#' }
#'
column_to_rownames <- function(.data, var = "rowname"){
  df <-
    as.data.frame(.data)  |>
    remove_rownames()
  if(!var %in% colnames(df)){
    stop("Variable '", var, "' not in data.", call. = FALSE)
  }
  rownames(df) <- df[[var]]
  df[[var]] <- NULL
  df
}
#' @name utils_rows_cols
#' @export
rownames_to_column <-  function(.data, var = "rowname"){
  col_names <- colnames(.data)
  if (var %in% col_names) {
    stop("Column `", var, "` already exists in `.data`.")
  }
  .data[, var] <- rownames(.data)
  rownames(.data) <- NULL
  .data[, c(var, setdiff(col_names, var))]
}
#' @name utils_rows_cols
#' @export
remove_rownames <- function(.data){
  rownames(.data) <- NULL
  .data
}

#' @name utils_rows_cols
#' @export
round_cols <- function(.data, digits = 2){
  num_col <- which(sapply(.data,  is.numeric))
  .data[num_col] <- apply(.data[num_col], 2, round, digits  = digits)
  return(.data)
}



#' Utilities for Principal Component Axis analysis
#' @description
#' * `pca()` Computes a Principal Component Analysis. It wrappers
#' [stats::prcomp()], but returns more results such as data, scores,
#' contributions and quality of measurements for individuals and variables.
#' * `get_biplot()`: Produces a biplot for an object computed with `pca()`.
#' * `plot.pca()`: Produces several types of plots, depending on the `type` and `which`
#' arguments.
#'    - `type = "var"` Produces a barplot with the contribution (`which =
#'    "contrib"`), qualitity of adjustment `which = "cos2"`, and a scatter plot
#'    with coordinates (`which = "coord"`) for the variables.
#'    - `type = "ind"` Produces a barplot with the contribution (`which =
#'    "contrib"`), qualitity of adjustment `which = "cos2"`, and a scatter plot
#'    with coordinates (`which = "coord"`) for the individuals.
#'    - `type = "biplot"` Produces a biplot.
#'
#' @param x
#' * For `pca()`, a numeric or complex matrix (or data frame) which provides the
#' data for the principal components analysis.
#' * For `plot.pca()` and `get_biplot()`, an object computed with `pca()`.
#' @param scale A logical value indicating whether the variables should be
#'   scaled to have unit variance before the analysis takes place. Defaults to
#'   `TRUE`.
#' @param axes The principal component axes to plot. Defaults to `axes = c(1, 2)`,
#'   i.e., the first and second interaction principal component axis.
#' @param show Which to show in the biplot. Defaults to `"both"` (both variables
#'   and individuals). One can also use `"var"`, or `"ind"`.
#' @param show_ind_id Shows the labels for individuals? Defaults to `TRUE`.
#' @param show_unit_circle Shows the unit variance circle? Defaults to `TRUE`.
#' @param expand An expansion factor to apply when plotting the second set of
#'   points relative to the first. This can be used to tweak the scaling of the
#'   two sets to a physically comparable scale. Setting to `TRUE` will
#'   automatically compute the expansion factor. Alternatively, a numeric value
#'   can be informed.
#' @param type One of `"var"` (to plot variables), `"ind"` (to plot
#'   individuals), or `"biplot"` to create a biplot.
#' @param which Which measure to plot. Either `which = "contribution"`
#'   (default),  `which = "cos2"` (quality of representation), or `which =
#'   "coord"` (coordinates)
#' @param axis The axist to plot the contribution/cos2. Defaults to 1.
#' @param ... Further arguments passed on to [get_biplot()] when `type =
#'   "biplot"`. Otherwise, When `which = "coord"`, further arguments passed on
#'   to [get_biplot()]. When `which = "contrib"`, or `which = "cos2"` further
#'   arguments passed on to [graphics::barplot()].
#' @return
#'  * `pca()` returns a list including:
#'    - `data`: The raw data used to compute the PCA.
#'    - `variances`: Variances (eigenvalues), and proportion of explained
#'    variance for each component.
#'    - `center,scale`: the centering and scaling used.
#'    - `ind,var` A list with the following objects for individuals/variables, respectively.
#'    - `coord`: coordinates for the individuals/variables (loadings * the
#' component standard deviations)
#'    - `cos2`: cos2 for the individuals/variables (coord^2)
#'    -  `contrib`: The contribution  (in percentage) of a variable to a given
#' principal component: (cos2 * 100) / (total cos2 of the component)
#'
#' * `plot.pca()` returns a list with the coordinates used.
#' * `get_biplot()` returns a `NULL` object
#' @export
#' @name utils_pca
#' @importFrom  stats lm cov
#' @importFrom  graphics arrows
#' @importFrom  graphics  axis box mtext plot.new plot.window polygon title
#'
#' @examples
#' library(pliman)
#' pc <- pca(mtcars[1:10 ,1:6])
#' plot(pc)
#' plot(pc, type = "ind")
#' plot(pc, type = "var", which = "coord")
#' plot(pc, type = "ind", which = "coord")
#' plot(pc, type = "biplot")


pca <- function(x, scale = TRUE){
  if(scale == TRUE){
    df <- scale(x)
    center <- attributes(df)$`scaled:center`
    scale <- attributes(df)$`scaled:scale`
  } else{
    df <- x
    center <- apply(df, 2, mean)
    scale <- NULL
  }
  # eigenvalues and variances
  covdf <- cov(df)
  eigenvalue <- eigen(covdf)$values
  std <- sqrt(eigenvalue)
  prop <- round(eigenvalue / sum(eigenvalue), 3)
  accum <- round(cumsum(prop), digits = 3)
  variances <- data.frame(eigenvalue, std, prop, accum)
  variances$PCA <- paste0("PC", 1:nrow(variances))
  variances <- variances[, c(5, 1:4)]

  s <- svd(df)
  U <- s$u
  D <- diag(s$d)
  V <- s$v
  loadings <- V |> as.data.frame()
  rownames(loadings) <- colnames(df)
  colnames(loadings) <- paste0("PC", 1:ncol(loadings))
  scores <- U %*% D |> as.data.frame()


  # variables
  var_cord <- sweep(V, 2, std, FUN = "*") |> as.data.frame()
  colnames(var_cord) <- paste0("PC", 1:ncol(var_cord))
  rownames(var_cord) <- colnames(df)

  var_coord_func <- function(loadings, comp.sdev){
    loadings*comp.sdev
  }
  var.coord <- t(apply(loadings, 1, var_coord_func, std))
  var.cos2 <- var.coord^2 |> as.data.frame()
  comp.cos2 <- apply(var.cos2, 2, sum)
  contribv <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2,1, contribv, comp.cos2)) |> as.data.frame()

  var <- list(
    contrib = var.contrib,
    cos2 = var.cos2,
    coord = var_cord,
    loadings = loadings
  )

  ind_cord <- data.frame(U %*% D)
  rownames(ind_cord) <- rownames(df)
  colnames(ind_cord) <- paste0("PC", 1:ncol(ind_cord))

  getdistance <- function(ind_row, center, scale){
    return(sum(((ind_row - center)/scale)^2))
  }
  d2 <- apply(df, 1, getdistance, center, scale)
  cos2 <- function(ind, d2){return(ind^2/d2)}
  ind.cos2 <- apply(scores, 2, cos2, d2) |> as.data.frame()
  colnames(ind.cos2) <- paste0("PC", 1:ncol(ind.cos2))
  contribi <- function(ind, comp.sdev, n.ind){
    100*(1/n.ind)*ind^2/comp.sdev^2
  }
  ind.contrib <- t(apply(scores, 1, contribi, std, nrow(scores))) |> as.data.frame()
  rownames(ind.contrib) <- rownames(df)
  colnames(ind.contrib) <- paste0("PC", 1:ncol(ind.contrib))
  ind <- list(
    contrib = ind.contrib,
    cos2 = ind.cos2,
    coord = ind_cord
  )

  return(structure(
    list(
      data = df,
      variances = variances,
      scale = scale,
      center = center,
      ind = ind,
      var = var

    ),
    class = "pca"
  ))
}

#' @export
#' @name utils_pca
get_biplot <- function(x,
                       axes = c(1, 2),
                       show = c("both"),
                       show_ind_id = TRUE,
                       show_unit_circle = TRUE,
                       expand = NULL){
  if(!show %in% c("both", "var", "ind")){
    stop("`show` must be one of 'var', 'ind', or 'both'", call. = FALSE)
  }
  if(!inherits(x, "pca")){
    stop("`x` must be an object computed with 'pca()'", call. = FALSE)
  }
  var <- x$var$coord
  ind <- x$ind$coord

  if(!is.null(expand)){
    if(isTRUE(expand)){
      expand <- min((max(ind[, axes[1]]) - min(ind[, axes[1]])/(max(var[, axes[1]]) - min(var[, axes[1]]))),
                    (max(ind[, axes[2]]) - min(ind[, axes[2]])/(max(var[, axes[2]]) - min(var[, axes[2]]))))

    }
    ind <- ind / expand
  }

  # adapted from https://bit.ly/3Lo34B0
  op <- par(pty = "s",
            cex.main = 1.2,
            cex.lab = 1,
            font.main = 2,
            font.lab = 2,
            family = "sans",
            col.main = "gray10",
            col.lab = "gray10",
            fg = "gray10",
            las = 1)
  on.exit(par(op))

  plot.new()

  if(show %in% c("both", "ind")){
    plot(x = ind[,c(axes[1], axes[2])],
         xlab = NA,
         ylab = NA,
         pch = 16,
         cex = 1,
         col = "red")
    if(isTRUE(show_ind_id)){
      text(x = ind[,c(axes[1], axes[2])],
           labels = row.names(ind),
           pos = 3,
           cex = 0.6)
    }
    if(show == "ind"){
      abline(v = 0, h = 0, lty = 2, col = "grey25")
    }
  }
  summ <- x$variances$prop
  title(xlab = paste("PC1 (", round(summ[axes[1]]*100, digits = 1), "%)", sep = ""),
        ylab = paste("PC1 (", round(summ[axes[2]]*100, digits = 1), "%)", sep = ""),
        line = 2.5,
        adj = 0.5)

  if(show %in% c("both", "var")){
    if (show == "var"){
      plot(x = var[,c(axes[1], axes[2])],
           xlim = c(-1, 1),
           ylim = c(-1, 1),
           xlab = NA,
           ylab = NA,
           asp = 1,
           type = "n")
      arrows(x0 = 0,
             x1 = var[[axes[1]]],
             y0 = 0,
             y1 = var[[axes[2]]],
             col = "blue",
             length = 0.08,
             lwd = 1,
             angle = 30)
      title(xlab = paste("PC1 (", round(summ[axes[1]]*100, digits = 1), "%)", sep = ""),
            ylab = paste("PC1 (", round(summ[axes[2]]*100, digits = 1), "%)", sep = ""),
            line = 2.5,
            adj = 0.5)
    } else{
      op2 <- par(new = TRUE, las = 1)
      on.exit(par(op2))
      plot.window(xlim = c(-1, 1),
                  ylim = c(-1, 1),
                  asp = 1)
      axis(side = 3,
           at = c(-1, 0.5, 0, -0.5, 1),
           labels = TRUE,
           col = "blue")
      axis(side = 4,
           at = c(-1, 0.5, 0, -0.5, 1),
           labels = TRUE,
           col = "blue")
      mtext((text = "PC1 rotations"),
            side = 3,
            cex = 1,
            font = 2,
            family = "sans",
            col = "blue",
            line = 2)
      mtext((text = "PC2 rotations"),
            side = 4,
            cex = 1,
            font = 2,
            family = "sans",
            col = "blue",
            line = 2,
            las = 3)
      box()
    }
    abline(v = 0, h = 0, lty = 2, col = "grey25")
    arrows(x0 = 0,
           x1 = var[[axes[1]]],
           y0 = 0,
           y1 = var[[axes[2]]],
           col = "blue",
           length = 0.08,
           lwd = 1,
           angle = 30)
    x_labs <- var[[axes[1]]]
    y_labs <- var[[axes[2]]]
    ylab_r <- data.frame(y_labs = c(1, -1), adj = c(-2, 2))
    mod <- lm(adj ~ y_labs, data = ylab_r)
    ylab_pos <- predict(mod, newdata = data.frame(y_labs))
    xlab_r <- data.frame(x_labs = c(-1, 1), adj = c(1.2, 0))
    mod <- lm(adj ~ x_labs, data = xlab_r)
    xlab_pos <- predict(mod, newdata = data.frame(x_labs))
    for (i in 1:length(x_labs)){
      text(x = x_labs[i],
           y = y_labs[i],
           labels = row.names(var)[i],
           adj = c(xlab_pos[i], ylab_pos[i]),
           cex = 0.7,
           col = "blue")
    }
    if(isTRUE(show_unit_circle)){
      ucircle = cbind(cos((0:360)/180*pi), sin((0:360)/180*pi))
      polygon(ucircle,
              lty = "solid",
              border = "blue",
              lwd = 1)
    }
  }
  invisible(list(ind, var))
}



#' @export
#' @name utils_pca
plot.pca <- function(x,
                     type = "var",
                     which = "contrib",
                     axis = 1,
                     ...){
  if(!which %in% c("contrib", "cos2", "coord")){
    stop("`show` must be one of 'contrib', 'cos2', or 'coord'", call. = FALSE)
  }
  if(!type %in% c("var", "ind", "biplot")){
    stop("`show` must be one of 'var', 'ind', or 'biplot'", call. = FALSE)
  }
  rotate_x <- function(data, column_to_plot, labels_vec, rot_angle, ...) {
    data <- data[order(data[,column_to_plot], decreasing = TRUE), ]
    plt <- barplot(data[[column_to_plot]],
                   col = 'steelblue',
                   xaxt = "n",
                   ...)
    text(plt,
         par("usr")[3],
         labels = labels_vec,
         srt = rot_angle,
         adj = c(1.1,1.1),
         xpd = TRUE,
         cex = 0.7)
  }


  if (type  %in% c("ind", "var") & which == "contrib"){
    x <- x[[type]][["contrib"]] |> as.data.frame()
    rotate_x(x, axis, row.names(x), 45,
             xlab = "Traits",
             ylab = "Contribution (%)")
    abline(h = mean(x[[axis]]), lty  = 2)
  }
  if (type  %in% c("ind", "var") & which == "cos2"){
    x <- x[[type]][["cos2"]] |> as.data.frame()
    rotate_x(x, axis, row.names(x), 45,
             xlab = "Traits",
             ylab = "Cos2 - Quality of representation")
  }

  if (type == "ind" & which == "coord"){
    get_biplot(x, show = "ind", ...)
  }
  if (type == "var" & which == "coord"){
    get_biplot(x, show = "var", ...)
  }

  if (type == "biplot"){
    get_biplot(x, ...)

  }
}




# Progress bar
# used in metan R package
# https://github.com/TiagoOlivoto/metan/blob/master/R/utils_progress.R
sec_to_hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":"
  )
}
progress <- function(min = 0,
                     max = 100,
                     leftd = "|",
                     rightd = "|",
                     char = "=",
                     style = 2,
                     width = getOption("width"),
                     time = Sys.time()){
  # Adapted from https://stackoverflow.com/a/26920123/15245107
  return(list(min = min,
              max = max,
              leftd = leftd,
              rightd = rightd,
              char = char,
              style = style,
              width = width,
              time = time))
}
run_progress <- function(pb,
                         actual,
                         text = "",
                         digits = 0,
                         sleep = 0){
  Sys.sleep(sleep)
  elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), pb$time, units = "secs")))
  temp <- switch(
    pb$style,
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
         text = paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""))),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
         text =  paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""), '% s%%')),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), elapsed)),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), '% s%%', elapsed))
  )
  step <- round(actual / pb$max * (pb$width - temp$extra))
  cat(sprintf(temp$text,
              strrep(pb$char, step),
              strrep(' ', pb$width - step - temp$extra),
              round(actual / pb$max * 100, digits = digits)), "\r")
  if(actual == pb$max){
    cat("\n")
  }
}

check_names_dir <- function(name, names_dir, dir){
  if(!name %in% names_dir){
    stop(paste("'", name, "' not found in '",
               paste(getwd(), sub(".", "", dir), sep = ""), "'", sep = ""),
         call. = FALSE)
  }
}

check_ebi <- function(){
  if(!requireNamespace("EBImage", quietly = TRUE)) {
    if(interactive() == TRUE){
      inst <-
        switch(menu(c("Yes", "No"), title = "Package {EBImage} required but not installed.\nDo you want to install it now?"),
               "yes", "no")
      if(inst == "yes"){
        if(!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", quiet = TRUE)
        }
        BiocManager::install("EBImage",
                             update = FALSE,
                             ask = FALSE,
                             quiet = TRUE)
      } else{
        message("To use {pliman}, first install {EBImage} following the directions at 'https://bioconductor.org/packages/EBImage'")
      }
    }
  }
}
# correct coordinates in analyze_objects_shp()
correct_coords <- function(coords, nrowimg, ncolimg, nrow, ncol){
  get_row_number <- function(vector, rows, cols) {
    row_numbers <- ((vector - 1) %/% cols) + 1
    return(row_numbers)
  }
  npixperrow <- ncolimg / nrow
  npixpercol <- nrowimg / ncol

  coords <-
    coords |>
    transform(plotn = as.numeric(gsub("[^0-9]", "", img)),
              row = get_row_number(as.numeric(gsub("[^0-9]", "", img)), nrow, ncol)) |>
    transform(col = ifelse(row == 1, plotn, plotn -  ncol * (row - 1))) |>
    transform(y = ifelse(row == 1, y, y + npixperrow * (row - 1)),
              x = ifelse(col == 1, x, x + npixpercol * (col - 1)))

  return(coords[, 1:4])
}




check_mapview <- function() {
  packages <- c("raster", "mapview", "mapedit", "leaflet", "leafem")
  rast <- !requireNamespace("raster", quietly = TRUE)
  mapv <- !requireNamespace("mapview", quietly = TRUE)
  mape <- !requireNamespace("mapedit", quietly = TRUE)
  leafl <- !requireNamespace("leaflet", quietly = TRUE)
  leafe <- !requireNamespace("leafem", quietly = TRUE)
  missing_packages <- packages[c(rast, mapv, mape, leafl, leafe)]
  if (length(missing_packages) > 0) {
    if (interactive()) {
      inst <- switch(menu(c("Yes", "No"),
                          title = paste("Packages", paste(missing_packages, collapse = ", "),
                                        "are required to use the `viewer = 'mapview' option`.\nDo you want to install them now?")),
                     "yes", "no")
      if (inst == "yes") {
        install.packages(missing_packages, quiet = TRUE)
      } else {
        message("To use viewer = 'mapview', first install the required packages:", paste(missing_packages, collapse = ", "))
      }
    } else {
      message("To use viewer = 'mapview', first install the required packages:", paste(missing_packages, collapse = ", "))
    }
  }
}


# get RGB values from a mask computed with EBImage::watershed()
get_rgb <- function(img, data_mask, index){
  data.frame(object = index,
             R = img@.Data[,,1][which(data_mask == index)],
             G = img@.Data[,,2][which(data_mask == index)],
             B = img@.Data[,,3][which(data_mask == index)])
}

# check for infinite values
check_inf <- function(data){
  indx <- apply(data, 2, function(x){
    any(is.na(x) | is.infinite(x))
  })
  if(any(indx) == TRUE){
    warning("Columns ", paste(colnames(data[indx]), collapse = ", "), " with infinite/NA values removed.", call. = FALSE)
  }
  data[,colnames(data[!indx])]
}


clear_td <- function(){
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
}

#' Turns a single character column into multiple columns.
#'
#' Given either a regular expression or a vector of character positions,
#' `separate_col()` turns a single character column into multiple columns.
#'
#' @param .data A data frame
#' @param col Column name
#' @param into Names of new variables to create as character vector
#' @param sep The separator between columns. By default, a regular expression
#'   that matches any sequence of non-alphanumeric values.
#'
#' @return A mutated `.data`
#' @export
#'
#' @examples
#' library(pliman)
#' df <- data.frame(x = paste0("TRAT_", 1:5),
#'                  y = 1:5)
#' df
#' separate_col(df, x, into = c("TRAT", "REP"))
separate_col <- function(.data, col, into, sep = "[^[:alnum:]]+"){
  var <- deparse(substitute(col))
  df <- strsplit(.data[[var]], split = sep)
  df <-
    do.call(rbind,
            lapply(df, function(x)x)) %>%
    as.data.frame()
  .data[[var]] <- NULL
  names(df) <- into
  return(cbind(df, .data))
}



#' Random built-in color names
#'
#' Randomly chooses single or multiple built-in color names which R knows about.
#' See more at [grDevices::colors()]
#'
#' @param n The number of color names. Defaults to 1.
#' @param distinct Logical indicating if the colors returned should all be
#'   distinct. Defaults to `FALSE`.
#'
#' @return A character vector of color names
#' @importFrom grDevices colors
#' @export
#'
#' @examples
#' library(pliman)
#' random_color(n = 3)
random_color <- function(n = 1, distinct = FALSE){
  replace <- ifelse(n > length(colors()), TRUE, FALSE)
  return(sample(colors(distinct = distinct), n, replace = replace))
}


#' ggplot2-like colors generation
#'
#' Generate ggplot2
#'
#' @param n The number of colors. This works well for up to about eight colours,
#'   but after that it becomes hard to tell the different colours apart.
#'
#' @importFrom grDevices hcl
#' @export
#' @examples
#'
#' library(pliman)
#' ggplot_color(n = 3)
ggplot_color <- function(n = 1){
  # adapted from https://stackoverflow.com/a/8197703
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Set and get the Working Directory quicky
#'
#' * [get_wd_here()] gets the working directory to the path of the current script.
#' * [set_wd_here()] sets the working directory to the path of the current script.
#' * [open_wd_here()] Open the File Explorer at the directory path of the current script.
#' * [open_wd()] Open the File Explorer at the current working directory.
#'
#' @param path Path components below the project root. Defaults to `NULL`. This means that
#'   the directory will be set to the path of the file. If the path doesn't exist, the
#'   user will be asked if he wants to create such a folder.
#' @return
#' * [get_wd_here()] returns a full-path directory name.
#' * [get_wd_here()] returns a message showing the current working directory.
#' * [open_wd_here()] Opens the File Explorer of the path returned by `get_wd_here()`.
#' @export
#' @name utils_wd
#' @examples
#'
#' \dontrun{
#' get_wd_here()
#' set_wd_here()
#' open_wd_here()
#' }
set_wd_here <- function(path = NULL){
  if(!requireNamespace("rstudioapi", quietly = TRUE)) {
    if(interactive() == TRUE){
      inst <-
        switch(menu(c("Yes", "No"), title = "Package {rstudioapi} required but not installed.\nDo you want to install it now?"),
               "yes", "no")
      if(inst == "yes"){
        install.packages("rstudioapi", quiet = TRUE)
      } else{
        message("To use `set_wd_here()`, first install {rstudioapi}.")
      }
    }
  } else{
    dir_path <- dirname(rstudioapi::documentPath())
    if(!is.null(path)){
      dir_path <- paste0(dir_path, "/", path)
    }
    d <- try(setwd(dir_path), TRUE)
    if(inherits(d, "try-error")){
      cat(paste0("Cannot change working directory to '", dir_path, "'."))
      done <- readline(prompt = "Do you want to create this folder now? (y/n) ")
      if(done == "y"){
        dir.create(dir_path)
        message("Directory '", dir_path, "' created.")
        setwd(dir_path)
        message("Working directory set to '", dir_path, "'")
      }
    } else{
      message("Working directory set to '", dir_path, "'")
    }
  }
}

#' @export
#' @name utils_wd
get_wd_here <- function(path = NULL){
  if(!requireNamespace("rstudioapi", quietly = TRUE)) {
    if(interactive() == TRUE){
      inst <-
        switch(menu(c("Yes", "No"), title = "Package {rstudioapi} required but not installed.\nDo you want to install it now?"),
               "yes", "no")
      if(inst == "yes"){
        install.packages("rstudioapi", quiet = TRUE)
      } else{
        message("To use `get_wd_here()`, first install {rstudioapi}.")
      }
    }
  } else{
    dir_path <- dirname(rstudioapi::documentPath())
    if(!is.null(path)){
      dir_path <- paste0(dir_path, "/", path)
    }
    dir_path
  }
}
#' @export
#' @name utils_wd
open_wd_here <- function(path = get_wd_here()){
  if(!requireNamespace("utils", quietly = TRUE)) {
    if(interactive() == TRUE){
      inst <-
        switch(menu(c("Yes", "No"), title = "Package {utils} required but not installed.\nDo you want to install it now?"),
               "yes", "no")
      if(inst == "yes"){
        install.packages("utils", quiet = TRUE)
      } else{
        message("To use `open_wd_here()`, first install {utils}.")
      }
    }
  } else{
    utils::browseURL(url = path)
  }
}

#' @export
#' @name utils_wd
open_wd <- function(path = getwd()){
  if(!requireNamespace("utils", quietly = TRUE)) {
    if(interactive() == TRUE){
      inst <-
        switch(menu(c("Yes", "No"), title = "Package {utils} required but not installed.\nDo you want to install it now?"),
               "yes", "no")
      if(inst == "yes"){
        install.packages("utils", quiet = TRUE)
      } else{
        message("To use `open_wd()`, first install {utils}.")
      }
    }
  } else{
    utils::browseURL(url = path)
  }
}

