% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_measures.R
\name{utils_measures}
\alias{utils_measures}
\alias{get_measures}
\alias{plot_measures}
\title{Utilities for object measures}
\usage{
get_measures(
  object,
  id = NULL,
  measure = NULL,
  dpi = NULL,
  verbose = TRUE,
  digits = 3
)

plot_measures(
  object,
  measure = "id",
  digits = 2,
  size = 0.9,
  col = "white",
  ...
)
}
\arguments{
\item{object}{An object computed with \code{\link[=count_objects]{count_objects()}} or \code{\link[=leaf_area]{leaf_area()}}.}

\item{id}{An object in the image to indicate a known value.}

\item{measure}{For \code{plot_measures()}, a character string; for
\code{get_measures()}, a two-sided formula, e.g., \code{measure = area ~ 100}
indicating the known value of object \code{id}. The right-hand side is the known
value and the left-hand side can be one of the following.
\itemize{
\item \code{area} The known area of the object.
\item \code{perimeter} The known perimeter of the object.
\item \code{radius_mean} The known radius of the object.
\item \code{radius_min} The known minimum radius of the object. If the object is a
square, then the \code{radius_min} of such object will be \code{L/2} where \code{L} is the
length of the square side.
\item \code{radius_max} The known maximum radius of the object. If the object is a
square, then the \code{radius_max} of such object according to the Pythagorean
theorem will be \verb{L x sqrt(2) / 2} where \code{L} is the length of the square side.
}}

\item{dpi}{A known resolution of the image in DPI (dots per inch).}

\item{verbose}{If \code{FALSE}, runs the code silently.}

\item{digits}{The number of significant figures. Defaults to \code{2.}}

\item{size}{The size of the text. Defaults to \code{0.9}.}

\item{col}{The color of the text. Defaults to \code{"white"}.}

\item{...}{Further arguments passed on to \code{\link[graphics:text]{graphics::text()}}.}
}
\value{
\itemize{
\item \code{get_measures()} returns a data frame with the object \code{id} and the
measures. If \code{measure} is informed, the pixel values will be corrected by the
value of the known object, given in the unit of the right-hand side of
\code{measure}. If \code{dpi} is informed, then all the measures will be adjusted to
the known \code{dpi}.
\item \code{plot_measures()} returns a \code{NULL} object, drawing the text according to
the x and y coordinates of the objects in \code{object}.
}
}
\description{
\itemize{
\item \code{get_measures()} computes object measures (area, perimeter, radius) by using
either a known resolution (dpi) or an object with known measurements.
\item \code{plot_measures()} draws the object measures given in an object to the
current plot. The object identification (\code{"id"}) is drawn by default.
}
}
\examples{
\donttest{
library(pliman)
img <- image_import(image_pliman("objects_300dpi.jpg"))
image_show(img)
# Image with four objects with a known resolution of 300 dpi
# Higher square: 10 x 10 cm
# Lower square: 5 x 5 cm
# Rectangle: 4 x 2 cm
# Circle: 3 cm in diameter

# Count the objects using the blue band to segment the image
results <-
   count_objects(img,
                 index = "B")
plot_measures(results, measure = "id")

# Get object measures by declaring the known resolution in dots per inch
(measures <- get_measures(results, dpi = 300))

# Calculated diagonal of the object 1
# 10 * sqrt(2) = 14.14

# Observed diagonal of the object 1
measures[1, "radius_max"] * 2


# Get object measures by declaring the known area of object 1
get_measures(results,
             id = 1,
             area ~ 100)
}
}
