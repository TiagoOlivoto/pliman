% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_objects.R
\name{utils_objects}
\alias{utils_objects}
\alias{object_coord}
\alias{object_isolate}
\alias{object_id}
\title{Utilities for working with image objects}
\usage{
object_coord(
  image,
  id = NULL,
  index = "NB",
  invert = FALSE,
  fill_hull = FALSE,
  edge = 2,
  extension = NULL,
  tolerance = NULL,
  object_size = "large",
  show_image = TRUE
)

object_isolate(image, id = NULL, ...)

object_id(image, ...)
}
\arguments{
\item{image}{An image of class \code{Image}.}

\item{id}{\itemize{
\item For \code{object_coord()}, a vector (or scalar) of object \code{id} to compute the
bounding rectangle. Object ids can be obtained with \code{\link[=object_id]{object_id()}}. Set \code{id = "all"} to compute the coordinates for all objects in the image. If \code{id = NULL} (default) a bounding rectangle is drawn including all the objects.
\item For \code{object_isolate()}, a scalar that identifies the object to be extracted.
}}

\item{index}{The index to produce a binary image used to compute bounding
rectangle coordinates. See \code{\link[=image_binary]{image_binary()}} for more details.}

\item{invert}{Inverts the binary image, if desired. Defaults to \code{FALSE}.}

\item{fill_hull}{Fill holes in the objects? Defaults to \code{FALSE}.}

\item{edge}{The number of pixels in the edge of the bounding rectangle.
Defaults to \code{2}.}

\item{extension, tolerance, object_size}{Controls the watershed segmentation of
objects in the image. See \code{\link[=count_objects]{count_objects()}} for more details.}

\item{show_image}{Shows the image with bounding rectangles? Defaults to
\code{TRUE}.}

\item{...}{\itemize{
\item For \code{object_isolate()}, further arguments passed on to \code{\link[=object_coord]{object_coord()}}.
\item For \code{object_id()}, further arguments passed on to \code{\link[=count_objects]{count_objects()}}.
}}
}
\value{
\itemize{
\item \code{object_id()} An image of class \code{"Image"} containing the object's
identification.
\item \code{object_coord()} A list with the coordinates for the bounding rectangles.
If \code{id = "all"} or a numeric vector, a list with a vector of coordinates is
returned.
\item \code{object_isolate()} An image of class \code{"Image"} containing the isolated
object.
}
}
\description{
\itemize{
\item \code{object_id()} get the object identification in an image.
\item \code{object_coord()} get the object coordinates and (optionally) draw a
bounding rectangle around multiple objects in an image.
\item \code{object_isolate()} isolates an object from an image.
}
}
\examples{
\donttest{
library(pliman)
img <- image_import(image_pliman("la_leaves.JPG"))
# Get the object's (leaves) identification
object_id(img)

# Get the coordinates and draw a bounding rectangle around leaves 1 and 3
object_coord(img, id = c(1, 3))

# Isolate leaf 3
isolated <- object_isolate(img, id = 3)
image_show(isolated)
}
}
