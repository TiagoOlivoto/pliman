#'Analyzes objects in an image
#'
#'@description
#' * [analyze_objects()] provides tools for counting and extracting object
#'features (e.g., area, perimeter, radius, pixel intensity) in an image. See
#'more at the **Details** section.
#' * [analyze_objects_iter()] provides an iterative section to measure object
#'features using an object with a known area.
#' * [plot.anal_obj()] produces a histogram for the R, G, and B values when
#'argument `object_index` is used in the function [analyze_objects()].
#'
#'@details A binary image is first generated to segment the foreground and
#'  background. The argument `index` is useful to choose a proper index to
#'  segment the image (see [image_binary()] for more details). It is also
#'  possible to provide color palettes for background and foreground (arguments
#'  `background` and `foreground`, respectively). When this is used, a general
#'  linear model (binomial family) fitted to the RGB values to segment fore- and
#'  background.
#'
#'  Then, the number of objects in the  foreground is counted. By setting up
#'  arguments such as `lower_size` and `upper_size`, it is possible to set a
#'  threshold for lower and upper sizes of the objects, respectively. The
#'  argument `object_size` can be used to set up pre-defined values of
#'  `tolerance` and `extension` depending on the image resolution. This will
#'  influence the watershed-based object segmentation. Users can also tune up
#'  `tolerance` and `extension` explicitly for a better precision of watershed
#'  segmentation.
#'
#'  If `watershed = FALSE` is used, all pixels for each connected set of
#'  foreground pixels in `img` are set to a unique object. This is faster,
#'  especially for a large number of objects, but it is not able to segment
#'  touching objects.
#'
#'  There are some ways to correct the measures based on a reference object. If
#'  a reference object with a known area (`reference_area`) is used in the image
#'  and `reference = TRUE` is used, the measures of the objects will be
#'  corrected, considering the unit of measure informed in `reference_area`.
#'  There are two main ways to work with reference objects.
#'    * The first, is to provide a reference object that has a contrasting color with
#'  both the background and object of interest. In this case, the arguments
#'  `back_fore_index` and `fore_ref_index` can be used to define an index to
#'  first segment the reference object and objects to be measured from the
#'  background, then the reference object from objects to be measured.
#'
#'
#'    * The second one is to use a reference object that has a similar color to the
#'  objects to be measured, but has a contrasting size. For example, if we are
#'  counting small brown grains, we can use a brown reference template that has
#'  an area larger (says 3 times the area of the grains) and then uses
#'  `reference_larger = TRUE`. With this, the larger object in the image will be
#'  used as the reference object. This is particularly useful when images are
#'  captured with background light, such as the example  2. Some types: (i) It
#'  is suggested that the reference object is not too much larger than the
#'  objects of interest (mainly when the `watershed = TRUE`). In some cases, the
#'  reference object can be broken into several pieces due to the watershed
#'  algorithm. (ii) Since the reference object will increase the mean area of
#'  the object, the argument `lower_noise` can be increased. By default
#'  (`lower_noise = 0.1`) objects with lesser than 10% of the mean area of all
#'  objects are removed. Since the mean area will be increased, increasing
#'  `lower_noise` will remove dust and noises more reliably. The argument
#'  `reference_smaller` can be used in the same way
#'
#'  By using `pattern`, it is possible to process several images with common
#'  pattern names that are stored in the current working directory or in the
#'  subdirectory informed in `dir_original`. To speed up the computation time,
#'  one can set `parallel = TRUE`.
#'
#'  [analyze_objects_iter()] can be used to process several images using an
#'  object with a known area as a template. In this case, all the images in the
#'  current working directory that match the `pattern` will be processed. For
#'  each image, the function will compute the features for the objects and show
#'  the identification (id) of each object. The user only needs to inform which
#'  is the id of the known object. Then, given the `known_area`, all the
#'  measures will be adjusted. In the end, a data.frame with the adjusted
#'  measures will be returned. This is useful when the images are taken at
#'  different heights. In such cases, the image resolution cannot be conserved.
#'  Consequently, the measures cannot be adjusted using the argument `dpi` from
#'  [get_measures()], since each image will have a different resolution. NOTE:
#'  This will only work in an interactive section.
#'
#' * Additional measures: By default, some measures are not computed, mainly due to
#'  computational efficiency when the user only needs simple measures such as
#'  area, length, and width.
#'
#'   - If `haralick = TRUE`, The function computes 13 Haralick texture features for
#'  each object based on a gray-level co-occurrence matrix (Haralick et al.
#'  1979). Haralick features depend on the configuration of the parameters
#'  `har_nbins` and `har_scales`. `har_nbins` controls the number of bins used
#'  to compute the Haralick matrix. A smaller `har_nbins` can give more accurate
#'  estimates of the correlation because the number of events per bin is higher.
#'  While a higher value will give more sensitivity. `har_scales` controls the
#'  number of scales used to compute the Haralick features. Since Haralick
#'  features compute the correlation of intensities of neighboring pixels it is
#'  possible to identify textures with different scales, e.g., a texture that is
#'  repeated every two pixels or 10 pixels. By default, the Haralick features
#'  are computed with the R band. To chance this default, use the argument
#'  `har_band`. For example, `har_band = 2` will compute the features with the
#'  green band.
#'
#'   - If `efourier = TRUE` is used, an Elliptical Fourier Analysis (Kuhl and
#'  Giardina, 1982) is computed for each object contour using [efourier()].
#'
#'    - If `veins = TRUE` (experimental), vein features are computed. This will call
#'  [object_edge()] and applies the Sobel-Feldman Operator to detect edges. The
#'  result is the proportion of edges in relation to the entire area of the
#'  object(s) in the image. Note that THIS WILL BE AN OPERATION ON AN IMAGE
#'  LEVEL, NOT an OBJECT LEVEL! So, If vein features need to be computed for
#'  leaves, it is strongly suggested to use one leaf per image.
#'
#'     - If `ab_angles = TRUE` the apex and base angles of each object are
#'  computed with [poly_apex_base_angle()]. By default, the function computes
#'  the angle from the first pixel of the apex of the object to the two pixels
#'  that slice the object at the 25th percentile of the object height (apex
#'  angle). The base angle is computed in the same way but from the first base
#'  pixel.
#'
#'@inheritParams image_binary
#'
#'@param img The image to be analyzed.
#'@param foreground,background A color palette for the foregrond and background,
#'  respectively (optional). If a chacarceter is used (eg., `foreground =
#'  "fore"`), the function will search in the current working directory a valid
#'  image named "fore".
#' @param pick_palettes  Logical argument indicating wheater the user needs to
#'   pick up the color palettes for foreground and background for the image. If
#'   `TRUE` [pick_palette()] will be called internally so that the user can sample
#'   color points representing foreground and background.
#'@param reference Logical to indicate if a reference object is present in the
#'  image. This is useful to adjust measures when images are not obtained with
#'  standard resolution (e.g., field images). See more in the details section.
#'@param reference_area The known area of the reference objects. The measures of
#'  all the objects in the image will be corrected using the same unit of the
#'  area informed here.
#'@param back_fore_index A character value to indicate the index to segment the
#'  foreground (objects and reference) from the background. Defaults to
#'  `"R/(G/B)"`. This index is optimized to segment white backgrounds from green
#'  leaves and a blue reference object.
#'@param fore_ref_index A character value to indicate the index to segment
#'  objects and the reference object. It can be either an available index in
#'  `pliman` (see [pliman_indexes()] or an own index computed with the R, G, and
#'  B bands. Defaults to `"B-R"`. This index is optimized to segment green
#'  leaves from a blue reference object after a white background has been
#'  removed.
#'@param reference_larger,reference_smaller Logical argument indicating when the
#'  larger/smaller object in the image must be used as the reference object.
#'  This only is valid when `reference` is set to `TRUE` and `reference_area`
#'  indicates the area of the reference object. IMPORTANT. When
#'  `reference_smaller` is used, objects with an area smaller than 1% of the
#'  mean of all the objects are ignored. This is used to remove possible noise
#'  in the image such as dust. So, be sure the reference object has an area that
#'  will be not removed by that cutpoint.
#'@param pattern A pattern of file name used to identify images to be imported.
#'  For example, if `pattern = "im"` all images in the current working directory
#'  that the name matches the pattern (e.g., img1.-, image1.-, im2.-) will be
#'  imported as a list. Providing any number as pattern (e.g., `pattern = "1"`)
#'  will select images that are named as 1.-, 2.-, and so on. An error will be
#'  returned if the pattern matches any file that is not supported (e.g.,
#'  img1.pdf).
#'@param parallel If `TRUE` processes the images asynchronously (in parallel) in
#'  separate R sessions running in the background on the same machine. It may
#'  speed up the processing time, especially when `pattern` is used is informed.
#'  When `object_index` is informed, multiple sections will be used to extract
#'  the RGB values for each object in the image. This may significantly speed up
#'  processing time when an image has lots of objects (say >1000).
#'@param workers A positive numeric scalar or a function specifying the number
#'  of parallel processes that can be active at the same time. By default, the
#'  number of sections is set up to 50% of available cores.
#'@param watershed If `TRUE` (default) performs watershed-based object
#'  detection. This will detect objects even when they are touching one other.
#'  If `FALSE`, all pixels for each connected set of foreground pixels are set
#'  to a unique object. This is faster but is not able to segment touching
#'  objects.
#'@param veins Logical argument indicating whether vein features are computed.
#'  This will call [object_edge()] and applies the Sobel-Feldman Operator to
#'  detect edges. The result is the proportion of edges in relation to the
#'  entire area of the object(s) in the image. Note that **THIS WILL BE AN
#'  OPERATION ON AN IMAGE LEVEL, NOT OBJECT!**.
#'@param sigma_veins Gaussian kernel standard deviation used in the gaussian
#'  blur in the edge detection algorithm
#' @param ab_angles  Logical argument indicating whether apex and base angles
#'   should be computed. Defaults to `FALSE`. If `TRUE`, `poly_apex_base_angle()`
#'   are called and the base and apex angles are computed considering the 25th
#'   and 75th percentiles of the object height. These percentiles can be changed
#'   with the argument `ab_angles_percentiles`.
#' @param ab_angles_percentiles The percentiles indicating the heights of the
#'   object for which the angle should be computed (from the apex and the
#'   bottom). Defaults to c(0.25, 0.75), which means considering the 25th and
#'   75th percentiles of the object height.
#'@param haralick Logical value indicating whether Haralick features are
#'  computed. Defaults to `FALSE`.
#'@param har_nbins An integer indicating the number of bins using to compute the
#'  Haralick matrix. Defaults to 32. See Details
#'@param har_scales A integer vector indicating the number of scales to use to
#'  compute the Haralick features. See Details.
#'@param har_band The band to compute the Haralick features (1 = R, 2 = G, 3 =
#'  B). Defaults to 1.
#' @param pcv Computes the Perimeter Complexity Value? Defaults to `FALSE`.
#' @param pcv_niter An integer specifying the number of smoothing iterations for
#'   computing the  Perimeter Complexity Value. Defaults to 100.
#'@param resize Resize the image before processing? Defaults to `FALSE`. Use a
#'  numeric value of range 0-100 (proportion of the size of the original image).
#'@param trim Number of pixels removed from edges in the analysis. The edges of
#'  images are often shaded, which can affect image analysis. The edges of
#'  images can be removed by specifying the number of pixels. Defaults to
#'  `FALSE` (no trimmed edges).
#'@param fill_hull Fill holes in the binary image? Defaults to `FALSE`. This is
#'  useful to fill holes in objects that have portions with a color similar to
#'  the background. IMPORTANT: Objects touching each other can be combined into
#'  one single object, which may underestimate the number of objects in an
#'  image.
#'@param filter Performs median filtering in the binary image? See more at
#'  [image_filter()]. Defaults to `FALSE`. Use a positive integer to define the
#'  size of the median filtering. Larger values are effective at removing noise,
#'  but adversely affect edges.
#'@param invert Inverts the binary image if desired. This is useful to process
#'  images with a black background. Defaults to `FALSE`. If `reference = TRUE`
#'  is use, `invert` can be declared as a logical vector of length 2 (eg.,
#'  `invert = c(FALSE, TRUE`). In this case, the segmentation of objects and
#'  reference from the foreground using `back_fore_index` is performed using the
#'  default (not inverted), and the segmentation of objects from the reference
#'  is performed by inverting the selection (selecting pixels higher than the
#'  threshold).
#'@param object_size The size of the object. Used to automatically set up
#'  `tolerance` and `extension` parameters. One of the following. `"small"`
#'  (e.g, wheat grains), `"medium"` (e.g, soybean grains), `"large"`(e.g, peanut
#'  grains), and `"elarge"` (e.g, soybean pods)`.
#'@param index A character value specifying the target mode for conversion to
#'  binary image when `foreground` and `background` are not declared. Defaults
#'  to `"NB"` (normalized blue). See [image_index()] for more details. User can
#'  also calculate your own index using the bands names, e.g. `index = "R+B/G"`
#'@param object_index Defaults to `FALSE`. If an index is informed, the average
#'  value for each object is returned. It can be the R, G, and B values or any
#'  operation involving them, e.g., `object_index = "R/B"`. In this case, it
#'  will return for each object in the image, the average value of the R/B
#'  ratio. Use [pliman_indexes_eq()] to see the equations of available indexes.
#' @param pixel_level_index Return the indexes computed in `object_index` in the
#'   pixel level? Defaults to `FALSE` to avoid returning large data.frames.
#' @param return_mask Returns the mask for the analyzed image? Defaults to `FALSE`.
#'@param efourier Logical argument indicating if Elliptical Fourier should be
#'  computed for each object. This will call [efourier()] internally. It
#'  `efourier = TRUE` is used, both standard and normalized Fourier coefficients
#'  are returned.
#'@param nharm An integer indicating the number of harmonics to use. Defaults to
#'  10. For more details see [efourier()].
#'@param tolerance The minimum height of the object in the units of image
#'  intensity between its highest point (seed) and the point where it contacts
#'  another object (checked for every contact pixel). If the height is smaller
#'  than the tolerance, the object will be combined with one of its neighbors,
#'  which is the highest.
#'@param extension Radius of the neighborhood in pixels for the detection of
#'  neighboring objects. Higher value smooths out small objects.
#'@param lower_noise To prevent noise from affecting the image analysis, objects
#'  with lesser than 10% of the mean area of all objects are removed
#'  (`lower_noise = 0.1`). Increasing this value will remove larger noises (such
#'  as dust points), but can remove desired objects too. To define an explicit
#'  lower or upper size, use the `lower_size` and `upper_size` arguments.
#'@param lower_size,upper_size Lower and upper limits for size for the image
#'  analysis. Plant images often contain dirt and dust.  Upper limit is set to
#'  `NULL`, i.e., no upper limit used. One can set a known area or use
#'  `lower_limit = 0` to select all objects (not advised). Objects that matches
#'  the size of a given range of sizes can be selected by setting up the two
#'  arguments. For example, if `lower_size = 120` and `upper_size = 140`,
#'  objects with size greater than or equal 120 and less than or equal 140 will
#'  be considered.
#'@param topn_lower,topn_upper Select the top `n` objects based on its area.
#'  `topn_lower` selects the `n` elements with the smallest area whereas
#'  `topn_upper` selects the `n` objects with the largest area.
#'@param lower_eccent,upper_eccent,lower_circ,upper_circ Lower and upper limit
#'  for object eccentricity/circularity for the image analysis. Users may use
#'  these arguments to remove objects such as square papers for scale (low
#'  eccentricity) or cut petioles (high eccentricity) from the images. Defaults
#'  to `NULL` (i.e., no lower and upper limits).
#'@param randomize Randomize the lines before training the model?
#'@param nrows The number of lines to be used in training step. Defaults to
#'  2000.
#'@param plot Show image after processing?
#'@param show_original Show the count objects in the original image?
#'@param show_chull Show the convex hull around the objects? Defaults to
#'  `FALSE`.
#'@param show_contour Show a contour line around the objects? Defaults to
#'  `TRUE`.
#'@param contour_col,contour_size The color and size for the contour line around
#'  objects. Defaults to `contour_col = "red"` and `contour_size = 1`.
#'@param show_lw If `TRUE`, plots the length and width lines on each object
#'  calling [plot_lw()].
#'@param show_background Show the background? Defaults to `TRUE`. A white
#'  background is shown by default when `show_original = FALSE`.
#'@param show_segmentation Shows the object segmentation colored with random
#'  permutations. Defaults to `FALSE`.
#'@param col_foreground,col_background Foreground and background color after
#'  image processing. Defaults to `NULL`, in which `"black"`, and `"white"` are
#'  used, respectively.
#'@param marker,marker_col,marker_size The type, color and size of the object
#'  marker. Defaults to `NULL`, which plots the object id. Use `marker =
#'  "point"` to show a point in each object or `marker = FALSE` to omit object
#'  marker.
#'@param save_image Save the image after processing? The image is saved in the
#'  current working directory named as `proc_*` where `*` is the image name
#'  given in `img`.
#'@param prefix The prefix to be included in the processed images. Defaults to
#'  `"proc_"`.
#'@param dir_original,dir_processed The directory containing the original and
#'  processed images. Defaults to `NULL`. In this case, the function will search
#'  for the image `img` in the current working directory. After processing, when
#'  `save_image = TRUE`, the processed image will be also saved in such a
#'  directory. It can be either a full path, e.g., `"C:/Desktop/imgs"`, or a
#'  subfolder within the current working directory, e.g., `"/imgs"`.
#'@param verbose If `TRUE` (default) a summary is shown in the console.
#'@param known_area The known area of the template object.
#'@param ... Depends on the function:
#' * For [analyze_objects_iter()], further arguments passed on to
#'  [analyze_objects()].
#'@return `analyze_objects()` returns a list with the following objects:
#'  * `results` A data frame with the following variables for each object in the
#'  image:
#'     - `id`:  object identification.
#'
#'     - `x`,`y`:  x and y coordinates for the center of mass of the object.
#'     - `area`:  area of the object (in pixels).
#'
#'     - `area_ch`:  the area of the convex hull around object (in pixels).
#'     - `perimeter`: perimeter (in pixels).
#'
#'     - `radius_min`, `radius_mean`, and `radius_max`: The minimum, mean, and
#'  maximum radius (in pixels), respectively.
#'
#'     - `radius_sd`: standard deviation of the mean radius (in pixels).
#'
#'     - `diam_min`, `diam_mean`, and `diam_max`: The minimum, mean, and
#'  maximum diameter (in pixels), respectively.
#'
#'     - `major_axis`, `minor_axis`: elliptical fit for major and minor axes (in
#'  pixels).
#'
#'     - `caliper`: The longest distance between any two points on the margin
#'  of the object. See [poly_caliper()] for more details
#'
#'     - `length`, `width` The length and width of objects (in pixels). These
#'  measures are obtained as the range of x and y coordinates after aligning
#'  each object with [poly_align()].
#'
#'     - `radius_ratio`: radius ratio given by `radius_max / radius_min`.
#'
#'     - `theta`: object angle (in radians).
#'
#'     - `eccentricity`: elliptical eccentricity computed using the
#'  ratio of the eigen values (inertia axes of coordinates).
#'
#'     - `form_factor` (Wu et al., 2007):  the difference between a leaf and a
#'  circle. It is defined as `4*pi*A/P`, where A is the area and P is the
#'  perimeter of the object.
#'
#'     - `narrow_factor` (Wu et al., 2007): Narrow factor (`caliper / length`).
#'
#'     - `asp_ratio` (Wu et al., 2007): Aspect ratio (`length / width`).
#'
#'     - `rectangularity` (Wu et al., 2007): The similarity between a leaf and
#'  a rectangle (`length * width/ area`).
#'
#'     - `pd_ratio` (Wu et al., 2007): Ratio of perimeter to diameter
#'  (`perimeter / caliper`)
#'
#'     - `plw_ratio` (Wu et al., 2007): Perimeter ratio of length and width
#'  (`perimeter / (length + width)`)

#'     - `solidity`: object solidity given by `area / area_ch`.
#'
#'     - `convexity`: The convexity of the object computed using the ratio
#'    between the perimeter of the convex hull and the perimeter of the polygon.
#'
#'     - `elongation`: The elongation of the object computed as `1 - width /
#'    length`.
#'
#'     - `circularity`: The object circularity given by `perimeter ^ 2 / area`.
#'
#'     - `circularity_haralick`: The Haralick's circularity (CH), computed as
#'     `CH =  m/sd`, where `m` and `sd` are the mean and standard deviations
#'     from each pixels of the perimeter to the centroid of the object.
#'
#'     - `circularity_norm`: The normalized circularity (Cn), to be unity for a
#'     circle. This measure is computed as `Cn = perimeter ^ 2 / 4*pi*area` and
#'     is invariant under translation, rotation, scaling transformations, and
#'     dimensionless.
#'
#'     - `asm`: The angular second-moment feature.
#'
#'     - `con`: The contrast feature
#'
#'     - `cor`: Correlation measures the linear dependency of gray levels of
#'     neighboring pixels.
#'
#'     - `var`: The variance of gray levels pixels.
#'
#'     - `idm`: The Inverse Difference Moment (IDM), i.e., the local
#'     homogeneity.
#'
#'     - `sav`: The Sum Average.
#'
#'     - `sva`: The Sum Variance.
#'
#'     - `sen`: Sum Entropy.
#'
#'     - `dva`: Difference Variance.
#'
#'     - `den`: Difference Entropy
#'
#'     - `f12`: Difference Variance.
#'
#'     - `f13`: The angular second-moment feature.
#'
#'  * `statistics`: A data frame with the summary statistics for the area of the
#'  objects.
#'  * `count`: If `pattern` is used, shows the number of objects in each image.
#'  * `obj_rgb`: If `object_index` is used, returns the R, G, and B values
#'  for each pixel of each object.
#'  * `object_index`: If `object_index` is used, returns the index computed for
#'  each object.
#'
#'  * Elliptical Fourier Analysis: If `efourier = TRUE` is used, the following
#'     objects are returned.
#'     - `efourier`: The Fourier coefficients.  For more details see
#'        [efourier()].
#'     - `efourier_norm`: The normalized Fourier coefficients. For more details
#'        see [efourier_norm()].
#'     - `efourier_error`: The error between original data and  reconstructed
#'        outline. For more details see [efourier_error()].
#'     - `efourier_power`: The spectrum of harmonic Fourier power.
#'        For more details see [efourier_power()].
#'
#'  * `veins`: If `veins = TRUE` is used, returns, for each image, the
#'  proportion of veins (in fact the object edges) related to the total object(s)' area.
#'
#'  * `analyze_objects_iter()` returns a data.frame containing the features
#'  described in the `results` object of [analyze_objects()].
#'
#'  * `plot.anal_obj()` returns a `trellis` object containing the distribution
#'  of the pixels, optionally  for each object when `facet = TRUE` is used.
#'
#' @references
#' Claude, J. (2008) \emph{Morphometrics with R}, Use R! series,
#' Springer 316 pp.
#'
#' Gupta, S., Rosenthal, D. M., Stinchcombe, J. R., & Baucom, R. S. (2020). The
#' remarkable morphological diversity of leaf shape in sweet potato (Ipomoea
#' batatas): the influence of genetics, environment, and G×E. New Phytologist,
#' 225(5), 2183–2195. \doi{10.1111/NPH.16286}
#'
#' Haralick, R.M., K. Shanmugam, and I. Dinstein. 1973. Textural Features for Image
#' Classification. IEEE Transactions on Systems, Man, and Cybernetics SMC-3(6): 610–621.
#' \doi{10.1109/TSMC.1973.4309314}
#'
#' Kuhl, F. P., and Giardina, C. R. (1982). Elliptic Fourier features of a
#' closed contour. Computer Graphics and Image Processing 18, 236–258. doi:
#' \doi{10.1016/0146-664X(82)90034-X}
#'
#' Lee, Y., & Lim, W. (2017). Shoelace Formula: Connecting the Area of a Polygon
#' and the Vector Cross Product. The Mathematics Teacher, 110(8), 631–636.
#' \doi{10.5951/mathteacher.110.8.0631}
#'
#' Montero, R. S., Bribiesca, E., Santiago, R., & Bribiesca, E. (2009). State
#' of the Art of Compactness and Circularity Measures. International
#' Mathematical Forum, 4(27), 1305–1335.
#'
#' Chen, C.H., and P.S.P. Wang. 2005. Handbook of Pattern Recognition and
#' Computer Vision. 3rd ed. World Scientific.
#'
#' Wu, S. G., Bao, F. S., Xu, E. Y., Wang, Y.-X., Chang, Y.-F., and Xiang, Q.-L.
#' (2007). A Leaf Recognition Algorithm for Plant Classification Using
#' Probabilistic Neural Network. in 2007 IEEE International Symposium on Signal
#' Processing and Information Technology, 11–16.
#' \doi{10.1109/ISSPIT.2007.4458016}
#'
#' @export
#' @name analyze_objects
#' @importFrom  utils install.packages
#' @importFrom grDevices col2rgb dev.off jpeg png rgb hcl.colors
#' @importFrom graphics lines par points rect text hist
#' @importFrom stats aggregate binomial glm kmeans predict sd runif dist var density
#' @importFrom utils menu
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(pliman)
#' img <- image_pliman("soybean_touch.jpg")
#' obj <- analyze_objects(img)
#' obj$statistics
#'
#' ########################### Example 1 #########################
#' # Enumerate the objects in the original image
#' # Return the top-5 grains with the largest area
#' top <-
#'  analyze_objects(img,
#'                  marker = "id",
#'                  topn_upper = 5)
#' top$results
#'
#'
#' #' ########################### Example 1 #########################
#' # Correct the measures based on the area of the largest object
#' # note that since the reference object
#'
#' img <- image_pliman("flax_grains.jpg")
#' res <-
#'   analyze_objects(img,
#'                   index = "GRAY",
#'                   marker = "point",
#'                   show_contour = FALSE,
#'                   reference = TRUE,
#'                   reference_area = 6,
#'                   reference_larger = TRUE,
#'                   lower_noise = 0.3)
#' }
#'
analyze_objects <- function(img,
                            foreground = NULL,
                            background = NULL,
                            pick_palettes = FALSE,
                            reference = FALSE,
                            reference_area = NULL,
                            back_fore_index = "R/(G/B)",
                            fore_ref_index = "B-R",
                            reference_larger = FALSE,
                            reference_smaller = FALSE,
                            pattern = NULL,
                            parallel = FALSE,
                            workers = NULL,
                            watershed = TRUE,
                            veins = FALSE,
                            sigma_veins = 1,
                            ab_angles = FALSE,
                            ab_angles_percentiles = c(0.25, 0.75),
                            haralick = FALSE,
                            har_nbins = 32,
                            har_scales = 1,
                            har_band = 1,
                            pcv = FALSE,
                            pcv_niter = 100,
                            resize = FALSE,
                            trim = FALSE,
                            fill_hull = FALSE,
                            filter = FALSE,
                            invert = FALSE,
                            object_size = "medium",
                            index = "NB",
                            object_index = NULL,
                            pixel_level_index = FALSE,
                            return_mask = FALSE,
                            efourier = FALSE,
                            nharm = 10,
                            threshold = "Otsu",
                            k = 0.1,
                            windowsize = NULL,
                            tolerance = NULL,
                            extension = NULL,
                            lower_noise = 0.10,
                            lower_size = NULL,
                            upper_size = NULL,
                            topn_lower = NULL,
                            topn_upper = NULL,
                            lower_eccent = NULL,
                            upper_eccent = NULL,
                            lower_circ = NULL,
                            upper_circ = NULL,
                            randomize = TRUE,
                            nrows = 1000,
                            plot = TRUE,
                            show_original = TRUE,
                            show_chull = FALSE,
                            show_contour = TRUE,
                            contour_col = "red",
                            contour_size = 1,
                            show_lw = FALSE,
                            show_background = TRUE,
                            show_segmentation = FALSE,
                            col_foreground = NULL,
                            col_background = NULL,
                            marker = FALSE,
                            marker_col = NULL,
                            marker_size = NULL,
                            save_image = FALSE,
                            prefix = "proc_",
                            dir_original = NULL,
                            dir_processed = NULL,
                            verbose = TRUE){
  lower_noise <- ifelse(isTRUE(reference_larger), lower_noise * 3, lower_noise)
  if(!object_size %in% c("small", "medium", "large", "elarge")){
    stop("'object_size' must be one of 'small', 'medium', 'large', or 'elarge'")
  }
  if(!missing(img) & !missing(pattern)){
    stop("Only one of `img` or `pattern` arguments can be used.", call. = FALSE)
  }
  if(is.null(dir_original)){
    diretorio_original <- paste0("./")
  } else{
    diretorio_original <-
      ifelse(grepl("[/\\]", dir_original),
             dir_original,
             paste0("./", dir_original))
  }
  if(is.null(dir_processed)){
    diretorio_processada <- paste0("./")
  } else{
    diretorio_processada <-
      ifelse(grepl("[/\\]", dir_processed),
             dir_processed,
             paste0("./", dir_processed))
  }
  help_count <-
    function(img, foreground, background, pick_palettes, resize, fill_hull, threshold, filter,
             tolerance, extension, randomize, nrows, plot, show_original,
             show_background, marker, marker_col, marker_size, save_image,
             prefix, dir_original, dir_processed, verbose, col_background,
             col_foreground, lower_noise, ab_angles, ab_angles_percentiles, return_mask, pcv){
      if(is.character(img)){
        all_files <- sapply(list.files(diretorio_original), file_name)
        check_names_dir(img, all_files, diretorio_original)
        imag <- list.files(diretorio_original, pattern = paste0("^",img, "\\."))
        name_ori <- file_name(imag)
        extens_ori <- file_extension(imag)
        img <- image_import(paste(name_ori, ".", extens_ori, sep = ""), path = diretorio_original)
      } else{
        name_ori <- match.call()[[2]]
        extens_ori <- "png"
      }
      if(trim != FALSE){
        if(!is.numeric(trim)){
          stop("Argument `trim` must be numeric.", call. = FALSE)
        }
        img <- image_trim(img, trim)
      }
      if(resize != FALSE){
        if(!is.numeric(resize)){
          stop("Argument `resize` must be numeric.", call. = FALSE)
        }
        img <- image_resize(img, resize)
      }
      # when reference is not used
      if(isFALSE(reference)){
        if(isTRUE(pick_palettes)){
          if(interactive()){
            plot(img)
            message("Use the first mouse button to pick up BACKGROUND colors. Press Est to exit")
            background <- pick_palette(img,
                                       r = 5,
                                       verbose = FALSE,
                                       palette  = FALSE,
                                       plot = FALSE,
                                       col = "blue")
            message("Use the first mouse button to pick up FOREGROUND colors. Press Est to exit")
            foreground <- pick_palette(img,
                                       r = 5,
                                       verbose = FALSE,
                                       palette  = FALSE,
                                       plot = FALSE,
                                       col = "salmon")
          }
        }
        if(!is.null(foreground) && !is.null(background)){
          if(is.character(foreground)){
            all_files <- sapply(list.files(getwd()), file_name)
            imag <- list.files(getwd(), pattern = foreground)
            check_names_dir(foreground, all_files, getwd())
            name <- file_name(imag)
            extens <- file_extension(imag)
            foreground <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
          }
          if(is.character(background)){
            all_files <- sapply(list.files(getwd()), file_name)
            imag <- list.files(getwd(), pattern = background)
            check_names_dir(background, all_files, getwd())
            name <- file_name(imag)
            extens <- file_extension(imag)
            background <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
          }
          original <-
            data.frame(CODE = "img",
                       R = c(img@.Data[,,1]),
                       G = c(img@.Data[,,2]),
                       B = c(img@.Data[,,3]))
          foreground <-
            data.frame(CODE = "foreground",
                       R = c(foreground@.Data[,,1]),
                       G = c(foreground@.Data[,,2]),
                       B = c(foreground@.Data[,,3]))
          background <-
            data.frame(CODE = "background",
                       R = c(background@.Data[,,1]),
                       G = c(background@.Data[,,2]),
                       B = c(background@.Data[,,3]))
          back_fore <-
            transform(rbind(foreground[sample(1:nrow(foreground)),][1:nrows,],
                            background[sample(1:nrow(background)),][1:nrows,]),
                      Y = ifelse(CODE == "background", 0, 1))

          formula <- as.formula(paste("Y ~ ", back_fore_index))

          modelo1 <- suppressWarnings(glm(formula,
                                          family = binomial("logit"),
                                          data = back_fore))
          pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
          foreground_background <- matrix(pred1, ncol = dim(img)[[2]])
          if(is.numeric(filter) & filter > 1){
            foreground_background <- EBImage::medianFilter(foreground_background, size = filter)
          }
          ID <- c(foreground_background == 1)
          ID2 <- c(foreground_background == 0)
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(foreground_background)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(foreground_background),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(foreground_background)
          }

        } else{
          img2 <- help_binary(img,
                              index = index,
                              invert = invert,
                              fill_hull = fill_hull,
                              threshold = threshold,
                              k = k,
                              windowsize = windowsize,
                              filter = filter,
                              resize = FALSE)
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img2)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(img2),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(img2)
          }
          ID <- which(img2 == 1)
          ID2 <- which(img2 == 0)
        }
        if(isTRUE(fill_hull)){
          nmask <- EBImage::fillHull(nmask)
        }
        shape <- compute_measures(mask = nmask,
                                  img = img,
                                  haralick = haralick,
                                  har_nbins = har_nbins,
                                  har_scales = har_scales,
                                  har_band = har_band)
        object_contour <- shape$cont
        ch <- shape$ch
        shape <- shape$shape

      } else{
        # when reference is used
        if(is.null(reference_area)){
          stop("A known area must be declared when a template is used.", call. = FALSE)
        }
        if(isFALSE(reference_larger) & isFALSE(reference_smaller)){
          # segment back and fore
          if(!isFALSE(invert)){
            invert1 <- ifelse(length(invert) == 1, invert, invert[1])
          } else{
            invert1 <- FALSE
          }
          img_bf <-
            help_binary(img,
                        threshold = threshold,
                        index = back_fore_index,
                        filter = filter,
                        k = k,
                        windowsize = windowsize,
                        invert = invert1,
                        fill_hull = fill_hull)
          img3 <- img
          img3@.Data[,,1][which(img_bf != 1)] <- 2
          img3@.Data[,,2][which(img_bf != 1)] <- 2
          img3@.Data[,,3][which(img_bf != 1)] <- 2
          ID <-  which(img_bf == 1) # IDs for foreground
          ID2 <- which(img_bf == 0) # IDs for background
          # segment fore and ref
          if(!isFALSE(invert)){
            invert2 <- ifelse(length(invert) == 1, invert, invert[2])
          } else{
            invert2 <- FALSE
          }
          img4 <-
            help_binary(img3,
                        threshold = threshold,
                        index = fore_ref_index,
                        filter = filter,
                        k = k,
                        windowsize = windowsize,
                        invert = invert2)
          mask <- img_bf
          pix_ref <- which(img4 != 1)
          img@.Data[,,1][pix_ref] <- 1
          img@.Data[,,2][pix_ref] <- 0
          img@.Data[,,3][pix_ref] <- 0
          npix_ref <- length(pix_ref)
          mask[pix_ref] <- 0
          if(is.numeric(filter) & filter > 1){
            mask <- EBImage::medianFilter(mask, size = filter)
          }
          if(isTRUE(watershed)){
            parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
            res <- length(img)
            parms2 <- parms[parms$object_size == object_size,]
            rowid <-
              which(sapply(as.character(parms2$resolution), function(x) {
                eval(parse(text=x))}))
            ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
            tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
            nmask <- EBImage::watershed(EBImage::distmap(mask),
                                        tolerance = tol,
                                        ext = ext)
          } else{
            nmask <- EBImage::bwlabel(mask)
          }
          shape <- compute_measures(mask = nmask,
                                    img = img,
                                    haralick = haralick,
                                    har_nbins = har_nbins,
                                    har_scales = har_scales,
                                    har_band = har_band)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape
          if(isTRUE(show_lw)){
            shape_ori <- shape
          }
          # correct measures based on the area of the reference object
          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side^2
          shape$area_ch <- shape$area_ch * px_side^2
          shape[6:18] <- apply(shape[6:18], 2, function(x){
            x * px_side
          })
        } else{
          # correct the measures based on larger or smaller objects

          if(!is.null(foreground) && !is.null(background) | isTRUE(pick_palettes)){

            if(isTRUE(pick_palettes)){
              if(interactive()){
                plot(img)
                message("Use the first mouse button to pick up BACKGROUND colors. Press Est to exit")
                background <- pick_palette(img,
                                           r = 5,
                                           verbose = FALSE,
                                           palette  = FALSE,
                                           plot = FALSE,
                                           col = "blue")
                message("Use the first mouse button to pick up FOREGROUND colors. Press Est to exit")
                foreground <- pick_palette(img,
                                           r = 5,
                                           verbose = FALSE,
                                           palette  = FALSE,
                                           plot = FALSE,
                                           col = "salmon")
              }
            }

            if(is.character(foreground)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = foreground)
              check_names_dir(foreground, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              foreground <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            if(is.character(background)){
              all_files <- sapply(list.files(getwd()), file_name)
              imag <- list.files(getwd(), pattern = background)
              check_names_dir(background, all_files, getwd())
              name <- file_name(imag)
              extens <- file_extension(imag)
              background <- image_import(paste(getwd(), "/", name, ".", extens, sep = ""))
            }
            original <-
              data.frame(CODE = "img",
                         R = c(img@.Data[,,1]),
                         G = c(img@.Data[,,2]),
                         B = c(img@.Data[,,3]))
            foreground <-
              data.frame(CODE = "foreground",
                         R = c(foreground@.Data[,,1]),
                         G = c(foreground@.Data[,,2]),
                         B = c(foreground@.Data[,,3]))
            background <-
              data.frame(CODE = "background",
                         R = c(background@.Data[,,1]),
                         G = c(background@.Data[,,2]),
                         B = c(background@.Data[,,3]))
            back_fore <-
              transform(rbind(foreground[sample(1:nrow(foreground)),][1:nrows,],
                              background[sample(1:nrow(background)),][1:nrows,]),
                        Y = ifelse(CODE == "background", 0, 1))

            formula <- as.formula(paste("Y ~ ", back_fore_index))

            modelo1 <- suppressWarnings(glm(formula,
                                            family = binomial("logit"),
                                            data = back_fore))
            pred1 <- round(predict(modelo1, newdata = original, type="response"), 0)
            foreground_background <- matrix(pred1, ncol = dim(img)[[2]])
            if(is.numeric(filter) & filter > 1){
              foreground_background <- EBImage::medianFilter(foreground_background, size = filter)
            }
            ID <- c(foreground_background == 1)
            ID2 <- c(foreground_background == 0)
            if(isTRUE(watershed)){
              parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
              res <- length(foreground_background)
              parms2 <- parms[parms$object_size == object_size,]
              rowid <-
                which(sapply(as.character(parms2$resolution), function(x) {
                  eval(parse(text=x))}))
              ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
              tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
              nmask <- EBImage::watershed(EBImage::distmap(foreground_background),
                                          tolerance = tol,
                                          ext = ext)
            } else{
              nmask <- EBImage::bwlabel(foreground_background)
            }

          } else{


            mask <-
              help_binary(img,
                          threshold = threshold,
                          index = index,
                          k = k,
                          windowsize = windowsize,
                          filter = filter,
                          invert = invert,
                          fill_hull = fill_hull)
            ID <-  which(mask == 1) # IDs for foreground
            ID2 <- which(mask == 0) # IDs for background
            if(isTRUE(watershed)){
              parms <- read.csv(file=system.file("parameters.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
              res <- length(mask)
              parms2 <- parms[parms$object_size == object_size,]
              rowid <-
                which(sapply(as.character(parms2$resolution), function(x) {
                  eval(parse(text=x))}))
              ext <- ifelse(is.null(extension),  parms2[rowid, 3], extension)
              tol <- ifelse(is.null(tolerance), parms2[rowid, 4], tolerance)
              nmask <- EBImage::watershed(EBImage::distmap(mask),
                                          tolerance = tol,
                                          ext = ext)
            } else{
              nmask <- EBImage::bwlabel(mask)
            }

          }


          shape <- compute_measures(mask = nmask,
                                    img = img,
                                    haralick = haralick,
                                    har_nbins = har_nbins,
                                    har_scales = har_scales,
                                    har_band = har_band)
          object_contour <- shape$cont
          ch <- shape$ch
          shape <- shape$shape
          if(isTRUE(show_lw)){
            shape_ori <- shape
          }

          if(isTRUE(reference_larger)){
            id_ref <- which.max(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
          } else{
            shape <- shape[shape$area > mean(shape$area) * lower_noise, ]
            id_ref <- which.min(shape$area)
            npix_ref <- shape[id_ref, 4]
            shape <- shape[-id_ref,]
          }

          px_side <- sqrt(reference_area / npix_ref)
          shape$area <- shape$area * px_side ^ 2
          shape$area_ch <- shape$area_ch * px_side ^ 2
          shape[6:18] <- apply(shape[6:18], 2, function(x){
            x * px_side
          })
        }
      }



      if(!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & !is.null(topn_upper)){
        stop("Only one of 'lower_*' or 'topn_*' can be used.")
      }
      ifelse(!is.null(lower_size),
             shape <- shape[shape$area > lower_size, ],
             shape <- shape[shape$area > mean(shape$area) * lower_noise, ])
      if(!is.null(upper_size)){
        shape <- shape[shape$area < upper_size, ]
      }
      if(!is.null(topn_lower)){
        shape <- shape[order(shape$area),][1:topn_lower,]
      }
      if(!is.null(topn_upper)){
        shape <- shape[order(shape$area, decreasing = TRUE),][1:topn_upper,]
      }
      if(!is.null(lower_eccent)){
        shape <- shape[shape$eccentricity > lower_eccent, ]
      }
      if(!is.null(upper_eccent)){
        shape <- shape[shape$eccentricity < upper_eccent, ]
      }
      if(!is.null(lower_circ)){
        shape <- shape[shape$circularity > lower_circ, ]
      }
      if(!is.null(upper_circ)){
        shape <- shape[shape$circularity < upper_circ, ]
      }
      object_contour <- object_contour[as.character(shape$id)]
      ch <- ch[as.character(shape$id)]

      # check if fourier is computed
      if(isTRUE(efourier)){
        efr <-
          efourier(object_contour,
                   nharm = nharm)
        efer <- efourier_error(efr, plot = FALSE)$stats
        efpowwer <- efourier_power(efr, plot = FALSE)
        efpow <- efpowwer$cum_power
        min_harm <- efpowwer$min_harm
        efrn <- efourier_norm(efr)
        efr <- efourier_coefs(efr)
        names(efr)[1] <- "id"
        efrn <- efourier_coefs(efrn)
        names(efrn)[1] <- "id"
      } else{
        efr <- NULL
        efrn <- NULL
        efer <- NULL
        efpow <- NULL
        min_harm <- NULL
      }

      # check if angles should be computed
      if(isTRUE(ab_angles)){
        angles <- poly_apex_base_angle(object_contour, ab_angles_percentiles)
      } else{
        angles <- NULL
      }

      # check if veins is computed
      if(isTRUE(veins)){
        vein <- object_edge(img, sigma = sigma_veins, plot = FALSE)
        prop_veins <- data.frame(prop_veins = sum(vein) / sum(shape$area))
      } else{
        prop_veins <- NULL
      }

      # check if perimeter complexity value is computed
      if(isTRUE(pcv)){
        pcv <- poly_pcv(object_contour, niter = pcv_niter)
      } else{
        pcv <- NULL
      }

      if(!is.null(object_index)){
        if(!is.character(object_index)){
          stop("`object_index` must be a character.", call. = FALSE)
        }
        ind <- read.csv(file=system.file("indexes.csv", package = "pliman", mustWork = TRUE), header = T, sep = ";")
        if(any(object_index %in% ind$Index)){
          ind_formula <- ind[which(ind$Index %in% object_index), 2]
        } else{
          ind_formula <- object_index
        }
        ind_name <- object_index
        data_mask <- nmask@.Data
        obj_rgb <- object_rgb(img, data_mask)
        obj_rgb <- subset(obj_rgb, id %in% shape$id)
        obj_rgb <- cbind(obj_rgb, rgb_to_hsb(obj_rgb[, 2:4]))
        # indexes by id
        tmp <-
          by(obj_rgb,
             INDICES = obj_rgb$id,
             FUN = function(x){
               data.frame(
                 do.call(cbind,
                         lapply(seq_along(ind_formula), function(i){
                           data.frame(transform(x, index = eval(parse(text = ind_formula[i])))[,8])
                         })
                 )
               )
             }
          )
        tmp <-
          do.call(rbind,
                  lapply(tmp, data.frame)
          )
        colnames(tmp) <- ind_name
        obj_rgb <- cbind(obj_rgb, tmp)
        indexes <- data.frame(cbind(id = obj_rgb$id, tmp))
        colnames(indexes) <- c("id", ind_name)
        indexes <- aggregate(. ~ id, indexes, mean, na.rm = TRUE)
        if(isFALSE(pixel_level_index)){
          obj_rgb <- NULL
        }
      } else{
        obj_rgb <- NULL
        indexes <- NULL
      }
      if(isTRUE(return_mask)){
        mask <- nmask
      } else{
        mask <- NULL
      }
      stats <- data.frame(stat = c("n", "min_area", "mean_area", "max_area",
                                   "sd_area", "sum_area", "coverage"),
                          value = c(length(shape$area),
                                    min(shape$area),
                                    mean(shape$area),
                                    max(shape$area),
                                    sd(shape$area),
                                    sum(shape$area),
                                    sum(shape$coverage)))
      results <- list(results = shape,
                      statistics = stats,
                      object_rgb = obj_rgb,
                      object_index = indexes,
                      efourier = efr,
                      efourier_norm = efrn,
                      efourier_error = efer,
                      efourier_power = efpow,
                      efourier_minharm = min_harm,
                      veins = prop_veins,
                      angles = angles,
                      mask = mask,
                      pcv = pcv,
                      parms = list(index = index))
      class(results) <- "anal_obj"
      if(plot == TRUE | save_image == TRUE){
        backg <- !is.null(col_background)
        # color for background
        if (is.null(col_background)){
          col_background <- col2rgb("white") / 255
        } else{
          ifelse(is.character(col_background),
                 col_background <- col2rgb(col_background) / 255,
                 col_background <- col_background / 255)
        }
        # color for lesions
        if (is.null(col_foreground)){
          col_foreground <- col2rgb("black") / 255
        } else{
          ifelse(is.character(col_foreground),
                 col_foreground <- col2rgb(col_foreground) / 255,
                 col_foreground <- col_foreground / 255)
        }

        if(show_original == TRUE & show_segmentation == FALSE){
          im2 <- img
          if(backg){
            im3 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im3@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im3@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im3@.Data[,,3]==0)] <- col_background[3]
          }
        }
        if(show_original == TRUE & show_segmentation == TRUE){
          im2 <- EBImage::colorLabels(nmask)
          if(backg){
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
          } else{
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- img@.Data[,,1][which(im2@.Data[,,1]==0)]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- img@.Data[,,2][which(im2@.Data[,,2]==0)]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- img@.Data[,,3][which(im2@.Data[,,3]==0)]
          }
        }
        if(show_original == FALSE){
          if(show_segmentation == TRUE){
            im2 <- EBImage::colorLabels(nmask)
            im2@.Data[,,1][which(im2@.Data[,,1]==0)] <- col_background[1]
            im2@.Data[,,2][which(im2@.Data[,,2]==0)] <- col_background[2]
            im2@.Data[,,3][which(im2@.Data[,,3]==0)] <- col_background[3]
          } else{
            im2 <- img
            im2@.Data[,,1][ID] <- col_foreground[1]
            im2@.Data[,,2][ID] <- col_foreground[2]
            im2@.Data[,,3][ID] <- col_foreground[3]
            im2@.Data[,,1][ID2] <- col_background[1]
            im2@.Data[,,2][ID2] <- col_background[2]
            im2@.Data[,,3][ID2] <- col_background[3]
          }
        }
        show_mark <- ifelse(isFALSE(marker), FALSE, TRUE)
        marker <- ifelse(is.null(marker), "id", marker)
        if(!isFALSE(show_mark) & marker != "point" & !marker %in% colnames(shape)){
          warning("Accepted 'marker' are: {", paste(colnames(shape), collapse = ", "),
                  "}. Drawing the object id.", call. = FALSE)
          marker <- "id"
        }
        marker_col <- ifelse(is.null(marker_col), "white", marker_col)
        marker_size <- ifelse(is.null(marker_size), 0.75, marker_size)
        if(plot == TRUE){
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2],
                   shape[, 3],
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
            if(isTRUE(show_chull)){
              plot_contour(ch |> poly_close(), col = "black")
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour)  & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2],
                     shape[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          # plot length and width
          if(isTRUE(show_lw)){
            if(isTRUE(reference)){
              plot_lw(shape_ori)
            } else{
              plot_lw(results)
            }
          }
        }

        if(save_image == TRUE){
          if(dir.exists(diretorio_processada) == FALSE){
            dir.create(diretorio_processada, recursive = TRUE)
          }
          png(paste0(diretorio_processada, "/",
                     prefix,
                     name_ori, ".",
                     extens_ori),
              width = dim(im2@.Data)[1],
              height = dim(im2@.Data)[2])
          if(marker != "point"){
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              text(shape[, 2],
                   shape[, 3],
                   round(shape[, marker], 2),
                   col = marker_col,
                   cex = marker_size)
            }
          } else{
            plot(im2)
            if(isTRUE(show_contour) & isTRUE(show_original)){
              plot_contour(object_contour, col = contour_col, lwd = contour_size)
            }
            if(show_mark){
              points(shape[, 2],
                     shape[, 3],
                     col = marker_col,
                     pch = 16,
                     cex = marker_size)
            }
          }
          if(isTRUE(show_lw)){
            if(isTRUE(reference)){
              plot_lw(shape_ori)
            } else{
              plot_lw(results)
            }
          }
          dev.off()
        }
      }
      invisible(results)
    }

  if(missing(pattern)){
    help_count(img, foreground, background, pick_palettes, resize, fill_hull, threshold, filter,
               tolerance , extension, randomize, nrows, plot, show_original,
               show_background, marker, marker_col, marker_size, save_image, prefix,
               dir_original, dir_processed, verbose, col_background,
               col_foreground, lower_noise, ab_angles, ab_angles_percentiles, return_mask, pcv)
  } else{
    if(pattern %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
      pattern <- "^[0-9].*$"
    }
    plants <- list.files(pattern = pattern, diretorio_original)
    extensions <- as.character(sapply(plants, file_extension))
    names_plant <- as.character(sapply(plants, file_name))
    if(length(grep(pattern, names_plant)) == 0){
      stop(paste("Pattern '", pattern, "' not found in '",
                 paste(getwd(), sub(".", "", diretorio_original), sep = ""), "'", sep = ""),
           call. = FALSE)
    }
    if(!all(extensions %in% c("png", "jpeg", "jpg", "tiff", "PNG", "JPEG", "JPG", "TIFF"))){
      stop("Allowed extensions are .png, .jpeg, .jpg, .tiff")
    }
    if(parallel == TRUE){
      init_time <- Sys.time()
      nworkers <- ifelse(is.null(workers), trunc(detectCores()*.5), workers)
      cl <- parallel::makePSOCKcluster(nworkers)
      doParallel::registerDoParallel(cl)
      on.exit(stopCluster(cl))

      if(verbose == TRUE){
        message("Processing ", length(names_plant), " images in multiple sessions (",nworkers, "). Please, wait.")
      }
      ## declare alias for dopar command
      `%dopar%` <- foreach::`%dopar%`

      results <-
        foreach::foreach(i = seq_along(names_plant), .packages = c("pliman", "EBImage")) %dopar%{
          help_count(names_plant[i],
                     foreground, background, pick_palettes, resize, fill_hull, threshold, filter,
                     tolerance , extension, randomize, nrows, plot, show_original,
                     show_background, marker, marker_col, marker_size, save_image, prefix,
                     dir_original, dir_processed, verbose, col_background,
                     col_foreground, lower_noise, ab_angles, ab_angles_percentiles, return_mask, pcv)
        }

    } else{
      init_time <- Sys.time()
      pb <- progress(max = length(plants), style = 4)
      foo <- function(plants, ...){
        if(verbose == TRUE){
          run_progress(pb, ...)
        }
        help_count(img  = plants,
                   foreground, background, pick_palettes, resize, fill_hull, threshold, filter,
                   tolerance, extension, randomize, nrows, plot, show_original,
                   show_background, marker, marker_col, marker_size, save_image,
                   prefix, dir_original, dir_processed, verbose, col_background,
                   col_foreground, lower_noise, ab_angles, ab_angles_percentiles, return_mask, pcv)
      }
      results <-
        lapply(seq_along(names_plant), function(i){
          foo(names_plant[i],
              actual = i,
              text = paste("Processing image", names_plant[i]))
        })
    }

    ## bind the results
    names(results) <- names_plant

    stats <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["statistics"]],
                          id =  names(results[i]))[,c(3, 1, 2)]
              })
      )

    if(!is.null(object_index)){
      if(!is.null(results[[1]][["object_rgb"]])){
        obj_rgb <-
          do.call(rbind,
                  lapply(seq_along(results), function(i){
                    transform(results[[i]][["object_rgb"]],
                              img =  names(results[i]))
                  })
          )
        obj_rgb <- obj_rgb[, c(ncol(obj_rgb), 1:ncol(obj_rgb) - 1)]
      } else{
        obj_rgb <- NULL
      }
      object_index <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["object_index"]],
                            img =  names(results[i]))
                })
        )
      object_index <- object_index[, c(ncol(object_index), 1:ncol(object_index) - 1)]
    } else{
      obj_rgb <- NULL
      object_index <- NULL
    }



    if(!isFALSE(efourier)){
      efourier <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier"]],
                            img =  names(results[i]))
                })
        )
      efourier <- efourier[, c(ncol(efourier), 1:ncol(efourier)-1)]
      names(efourier)[2] <- "id"

      efourier_norm <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_norm"]],
                            img =  names(results[i]))
                })
        )
      efourier_norm <- efourier_norm[, c(ncol(efourier_norm), 1:ncol(efourier_norm)-1)]
      names(efourier_norm)[2] <- "id"


      efourier_error <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_error"]],
                            img =  names(results[i]))
                })
        )
      efourier_error <- efourier_error[, c(ncol(efourier_error), 1:ncol(efourier_error)-1)]
      names(efourier_error)[2] <- "id"

      efourier_power <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_power"]],
                            img =  names(results[i]))
                })
        )
      efourier_power <- efourier_power[, c(ncol(efourier_power), 1:ncol(efourier_power)-1)]
      names(efourier_power)[2] <- "id"

      efourier_minharm <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["efourier_minharm"]],
                            img =  names(results[i]))
                })
        )
      efourier_minharm <- efourier_minharm[, c(ncol(efourier_minharm), 1:ncol(efourier_minharm)-1)]
      names(efourier_minharm)[2] <- "id"

    } else{
      efourier <- NULL
      efourier_norm <- NULL
      efourier_error <- NULL
      efourier_power <- NULL
      efourier_minharm <- NULL
    }



    if(isTRUE(veins)){
      veins <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["veins"]],
                            img =  names(results[i]))
                })
        )

      veins <- veins[, c(ncol(veins), 1:ncol(veins)-1)]
    } else{
      veins <- NULL
    }



    if(isTRUE(ab_angles)){
      angles <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  transform(results[[i]][["angles"]],
                            img =  names(results[i]))
                })
        )

      angles <- angles[, c(ncol(angles), 1:ncol(angles)-1)]
    } else{
      angles <- NULL
    }

    if(isTRUE(pcv)){
      pcv <-
        do.call(rbind,
                lapply(seq_along(results), function(i){
                  data.frame(pcv = results[[i]][["pcv"]]) |>
                    transform(img =  names(results[i]))
                })
        )

      pcv <- pcv[, c("img", "pcv")]
    } else{
      pcv <- NULL
    }

    results <-
      do.call(rbind,
              lapply(seq_along(results), function(i){
                transform(results[[i]][["results"]],
                          img =  names(results[i]))
              })
      )

    if("img" %in% colnames(results)){
      results <- results[, c(ncol(results), 1:ncol(results) - 1)]
    }
    summ <- stats[stats$stat == "n", c(1, 3)]
    if(verbose == TRUE){
      names(summ) <- c("Image", "Objects")
      cat("--------------------------------------------\n")
      print(summ, row.names = FALSE)
      cat("--------------------------------------------\n")
      message("Done!")
      message("Elapsed time: ", sec_to_hms(as.numeric(difftime(Sys.time(),  init_time, units = "secs"))))

    }

    invisible(
      structure(
        list(statistics = stats,
             count = summ,
             results = results,
             obj_rgb = obj_rgb,
             object_index = object_index,
             efourier = efourier,
             efourier_norm = efourier_norm,
             efourier_error = efourier_error,
             efourier_minharm = efourier_minharm,
             veins = veins,
             angles = angles,
             pcv = pcv),
        class = "anal_obj_ls"
      )
    )
  }
}


#' @name analyze_objects
#' @param x An object of class `anal_obj`.
#' @param which Which to plot. Either 'measure' (object measures) or 'index'
#'   (object index). Defaults to `"measure"`.
#' @param measure The measure to plot. Defaults to `"area"`.
#' @param type The type of plot. Either `"hist"` or `"density"`. Partial matches
#'   are recognized.
#' @method plot anal_obj
#' @export
#'
#' @examples
#' \donttest{
#' library(pliman)
#'
#' img <- image_pliman("soy_green.jpg")
#' # Segment the foreground (grains) using the normalized blue index (NB, default)
#' # Shows the average value of the blue index in each object
#'
#' rgb <-
#'    analyze_objects(img,
#'                    marker = "id",
#'                    object_index = "B",
#'                    pixel_level_index = TRUE)
#' # density of area
#' plot(rgb)
#'
#' # histogram of perimeter
#' plot(rgb, measure = "perimeter", type = "histogram") # or 'hist'
#'
#' # density of the blue (B) index
#' plot(rgb, which = "index")
#' }
plot.anal_obj <- function(x,
                          which = "measure",
                          measure = "area",
                          type = c("density", "histogram"),
                          ...){
  if(!which %in% c("measure", "index")){
    stop("'which' must be one of 'measure' or 'index'", call. = FALSE)
  }
  if(which == "measure"){
    nam <- colnames(x$results)
    if(!measure %in% nam){
      stop("Measure '", measure, "' not available in 'x'. Try one of the '",
           paste0(nam, collapse = ", "), call. = FALSE)
    }
    temp <- x$results[[measure]]
    types <- c("density", "histogram")
    matches <- grepl(type[1], types)
    type <- types[matches]
    if(type == "histogram"){
      hist(temp,  xlab = paste(measure), main = NA, col = "cyan")
    } else{
      density_data <- density(temp)  # Calculate the density for the column
      plot(density_data, col = "red", main = NA, lwd = 2, xlab = paste(measure), ylab = "Density")  # Create the density plot
      points(x = temp, y = rep(0, length(temp)), col = "red")
    }
  } else{
    rgb <- x$object_rgb
    if(is.null(rgb)){
      stop("RGB values not found. Use `object_index` in the function `analyze_objects()`.\nHave you accidentally missed the argument `pixel_level_index = TRUE`?", call. = FALSE)
    }
    plot(density(rgb$R),
         main = NA,
         col = "red",
         lwd = 2,
         xlim = c(min(rgb$R, rgb$G, rgb$B), max(rgb$R, rgb$G, rgb$B)),
         ylim = c(0, max(density(rgb$R)$y, density(rgb$G)$y, density(rgb$B)$y)),
         xlab = "Pixel value",
         ylab = "Density")

    # Add the density curves for G and B
    lines(density(rgb$G), col = "green", lwd = 2)
    lines(density(rgb$B), col = "blue", lwd = 2)
    # Add a legend
    legend("topright", legend = c("R", "G", "B"), col = c("red", "green", "blue"), lty = 1,
           lwd = 2)

  }
}



#' @export
#' @name analyze_objects
analyze_objects_iter <- function(pattern,
                                 known_area,
                                 verbose = TRUE,
                                 ...){
  if (interactive()) {
    imgs <- list.files(pattern = pattern)
    measures <- list()
    for (i in 1:length(imgs)) {
      tmp <-
        analyze_objects(img = file_name(imgs[[i]]),
                        marker = "id",
                        ...)
      object <- as.numeric(readline("Known object: "))
      tmp <- get_measures(tmp,
                          id = object,
                          measure = area ~ known_area,
                          verbose = verbose)
      tmp$img <- file_name(imgs[[i]])
      tmp <- tmp[, c(ncol(tmp), 1:ncol(tmp) - 1)]
      measures[[i]] <- tmp
    }
    do.call(rbind, lapply(measures, function(x){x}))
  }
}


