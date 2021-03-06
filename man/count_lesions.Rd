% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_lesions.R
\name{count_lesions}
\alias{count_lesions}
\title{Counts the number of lesions}
\usage{
count_lesions(
  img,
  img_healthy = NULL,
  img_lesion = NULL,
  img_background = NULL,
  img_pattern = NULL,
  parallel = FALSE,
  workers = NULL,
  resize = FALSE,
  fill_hull = TRUE,
  invert = FALSE,
  index = "NB",
  my_index = NULL,
  lower_size = NULL,
  upper_size = NULL,
  randomize = TRUE,
  nrows = 10000,
  lesion_size = "medium",
  segment = TRUE,
  tolerance = NULL,
  extension = NULL,
  show_segmentation = TRUE,
  show_image = FALSE,
  show_original = TRUE,
  show_background = TRUE,
  col_leaf = NULL,
  col_lesions = NULL,
  col_background = NULL,
  marker = NULL,
  marker_col = NULL,
  marker_size = NULL,
  save_image = FALSE,
  prefix = "proc_",
  dir_original = NULL,
  dir_processed = NULL,
  verbose = TRUE
)
}
\arguments{
\item{img}{The image to be analyzed.}

\item{img_healthy}{A color palette of healthy areas.}

\item{img_lesion}{A color palette of lesioned areas.}

\item{img_background}{An optional color palette of the image background.}

\item{img_pattern}{A pattern of file name used to identify images to be
processed. For example, if \code{img_pattern = "im"} all images that the name
matches the pattern (e.g., img1.-, image1.-, im2.-) will be analyzed.
Providing any number as pattern (e.g., \code{img_pattern = "1"}) will select
images that are named as 1.-, 2.-, and so on.}

\item{parallel}{Processes the images asynchronously (in parallel) in separate
R sessions running in the background on the same machine. It may speed up
the processing time, especially when \code{img_pattern} is used is informed. The
number of sections is set up to 80\% of available cores.}

\item{workers}{A positive numeric scalar or a function specifying the maximum
number of parallel processes that can be active at the same time.}

\item{resize}{Resize the image before processing? Defaults to \code{FALSE}. Use a
numeric value of range 0-100 (proportion of the size of the original
image).}

\item{fill_hull}{Fill holes in the image? Defaults to \code{TRUE}. This is useful
to fill holes in leaves, e.g., those caused by insect attack, ensuring the
hole area will be accounted for the leaf, not background.}

\item{invert}{Inverts the binary image, if desired. This is useful to process
images with black background. Defaults to \code{FALSE}.}

\item{index, my_index}{A character value specifying the target mode for
conversion to binary image when \code{img_healthy} and \code{img_lesion} are not
declared. Defaults to \code{"NB"} (normalized blue). See \code{\link[=image_index]{image_index()}} for
more details.}

\item{lower_size}{Lower limit for size for the image analysis. Leaf images
often contain dirt and dust. To prevent dust from affecting the image
analysis, the lower limit of analyzed size is set to 0.1, i.e., objects
with lesser than 10\% of the mean of all objects are removed. One can set a
known area or use \code{lower_limit = 0} to select all objects (not advised).}

\item{upper_size}{Upper limit for size for the image analysis. Defaults to
\code{NULL}, i.e., no upper limit used.}

\item{randomize}{Randomize the lines before training the model?}

\item{nrows}{The number of lines to be used in training step.}

\item{lesion_size}{The size of the lesion. Used to automatically set up
\code{tolerance} and \code{extension} parameters. One of the following. \code{"small"} (2-5
mm in diameter, e.g, rust pustules), \code{"medium"} (0.5-1.0 cm in diameter,
e.g, wheat leaf spot), \code{"large"} (1-2 cm in diameter, and  \code{"elarge"} (2-3
cm in diameter, e.g, target spot of soybean).}

\item{segment}{If \code{TRUE} (Default) implements the Watershed Algorithm to
segment lesions connected by a fairly few pixels that could be considered
as two distinct lesions. If \code{FALSE}, lesions that are connected by any
pixel are considered unique lesions. For more details see
\code{\link[EBImage:watershed]{EBImage::watershed()}}.}

\item{tolerance}{The minimum height of the object in the units of image
intensity between its highest point (seed) and the point where it contacts
another object (checked for every contact pixel). If the height is smaller
than the tolerance, the object will be combined with one of its neighbors,
which is the highest. Defaults to \code{NULL}, i.e., starting values are set up according to the argument \code{lesion_size}.}

\item{extension}{Radius of the neighborhood in pixels for the detection of
neighboring objects. Defaults to 20. Higher value smooths out small
objects.}

\item{show_segmentation}{Shows the object segmentation colored with random
permutations. Defaults to \code{TRUE}.}

\item{show_image}{Show image after processing?}

\item{show_original}{Show the symptoms in the original image?}

\item{show_background}{Show the background? Defaults to \code{TRUE}. A white
background is shown by default when \code{show_original = FALSE}.}

\item{col_leaf}{Leaf color after image processing. Defaults to \code{"green"}}

\item{col_lesions}{Symptoms color after image processing. Defaults to \code{"red"}.}

\item{col_background}{Background color after image processing. Defaults to
\code{"NULL"}.}

\item{marker, marker_col, marker_size}{The type, color and size of the object
marker. Defaults to \code{NULL}, which shows a red point when \code{show_segmentation = FALSE}. To force a marker to be used with segmented objects, set up to
\code{marker = "point"} (to show a point) or \code{marker = "text"} to enumerate the
objects.}

\item{save_image}{Save the image after processing? The image is saved in the
current working directory named as \verb{proc_*} where \code{*} is the image name
given in \code{img}.}

\item{prefix}{The prefix to be included in the processed images. Defaults to
\code{"proc_"}.}

\item{dir_original, dir_processed}{The directory containing the original and
processed images. Defaults to \code{NULL}. In this case, the function will
search for the image \code{img} in the current working directory. After
processing, when \code{save_image = TRUE}, the processed image will be also
saved in such a directory.}

\item{verbose}{If \code{TRUE} (default) a summary is shown in the console.}
}
\value{
A list with the following objects:
\itemize{
\item \code{results} A data frame with the results (area, perimeter, radius) for
object.
\item \code{statistics} A data frame with the summary statistics for the image.
\item \code{count} (If \code{img_pattern} is used), summarizing the count number for each
image.
}
}
\description{
Counts the number of lesions in a sample or entire leaf based on provided
color palettes samples. A general linear model (binomial family) fitted to the
RGB values is used to segment the lesions from the healthy leaf. If a pallet
of background is provided, the function takes care of the details to isolate
it before computing the number and area of lesions. By using \code{img_pattern} it
is possible to process several images with common pattern names that are
stored in the current working directory or in the subdirectory informed in
\code{dir_original}.
}
\examples{
\donttest{
library(pliman)
img <- image_import(image_pliman("sev_leaf_nb.jpg"))
healthy <- image_import(image_pliman("sev_healthy.jpg"))
lesions <- image_import(image_pliman("sev_sympt.jpg"))
image_combine(img, healthy, lesions, ncol = 3)
a <-
count_lesions(img = img,
              img_healthy = healthy,
              img_lesion = lesions,
              lesion_size = "elarge", # extra large lesions
              show_image = TRUE,
              show_segmentation = FALSE,
              marker = "text")
}

}
\author{
Tiago Olivoto \email{tiagoolivoto@gmail.com}
}
