# pliman (development version)
## New functions
* `analyze_objects_iter()` to execute an interactive section of `analyze_objects()`.
* `measure_disease_byl()` to measure disease severity '`by` `l`eaf' in an image with several leaves.
* `object_split()` to split multiples objects of an image into a list of images.
* `pca()`, `plot.pca()`, `get_biplot()` as helper functions to perform Principal Component Analysis.
* `rownames_to_column()`, `column_to_rownames()`, `separate_col()`, `round_cols()` as helper functions to manipulate data.

* A set of `poly_*()` function to analyze polygons. All of them are based on a set of coordinate points describing the edge of the object(s). See ?`utils_polygon` for more details.

* `get_na_here()` and `set_na_here()` to deal with working directories.

* `apply_fun_to_imgs()` to apply a function (or functions) to a set of images stored in the working directory.

* `make_brush()`, `make_mask()`, and `image_segment_mask()` to create masks and segment images based on such a mask.

## New features
* Several measures were added in `analyze_objects()`. The function now wraps some `poly_*()` functions to compute shape measures such as width, length, elongation, and circularity. This improvement was at cost of a slight increase in computation time.

* Haralick's features that quantify pixel texture for image objects were included.

* `analyze_objects()`, `measure_disease()`, and `measure_disease_byl()` have now a `filter` argument that applies a median filtering in the binary mask. This is useful to reduce the noise in the segmentation of objects.


## Minor changes
* `get_measures()` now remove known objects from the results when using the `id` argument.
* The right-hand of argument `measure` in `get_measures()` now accepts a numeric object stored in the global environment.
* `analyze_objects()` now returns the objects `object_index` and `object_rgb` when the argument `pattern` is used. Thanks to [João Paulo Oliveira Ribeiro](https://www.researchgate.net/profile/Joao-Paulo-Ribeiro) for alerting me regarding this issue.
* New argument `reference` in `analyze_objects()` to adjust measures using a reference object in the image.
* Argument `object_index` in `analyze_objects()` now recognizes the names of built-in indexes (see ?`pliman_indexes()`).


# pliman 1.1.0
## New functions
* `measure_disease_iter()` to measure disease in an interactive section.
* `pick_count()` to count objects in an image manually.
* `pick_palette()` to create an image palette by picking up color point(s) from the image
* `pick_rgb()` to pick up the RGB values from selected point(s) in the image.
* `summary_index()` to summary the index either between and within objects.
* `pliman` now exports the foward-pipe operator `%>%`. Code from [poorman](https://github.com/nathaneastwood/poorman/blob/master/R/pipe.R) package.

## Minor changes
* Deprecated functions in the last version (`count_objects()`, `image_show()`, `leaf_area()`, `objects_rgb()`, `prop_segmented()`, and `symptomatic_area()`) were removed.
* Use Bootstrap 5 from [pkgdown 2.0.0](https://pkgdown.r-lib.org/news/index.html#bootstrap-2-0-0) in the package site.


# pliman 1.0.0
## New functions
* `analyze_objects()` is now used as the main function to compute the number and shape of objects.
* `measure_disease()` is now used as the main function to perform phytopatometry studies. The function can compute symptomatic area, as well as the number and shape of lesions.
* `image_segment_iter()` is used to performs iterative image segmentation.
* `conv_hull()`, `poly_area()`, `poly_mass()`, `poly_spline()`, `plot_contour()`, and `plot_ellipse()` as utilities for analyzing polygons.
* `dpi()` to compute the resolution (dots per inch) of an image.
* `tune_tolerance()` for tunning the `tolerance` parameter.

## Deprecated functions.
* `objects_rgb()` will be depracated in the future. Now, to compute an index for each object use the `object_index` argument in  `analyze_objects()`, for example, `analyze_objects(object_index = "B")`.
* `leaf_area()` will be depracated in the future. Now, combine `analyze_objects()` with `get_measures()` to obtain the area and shape of objects (leaves).
* `prop_segmented()` is now deprecated in favour of `image_segment_iter()`.
* `count_lesions()` is now deprecated. Now, to compute the number and shape of lesions, use the argument `show_features = TRUE` in `measure_disease()`.
* `image_show()` is now deprecated in favour of `plot()`.
   
## Minor improvements
* Include `fill_hull` argument in `symptomatic_area()` and `count_lesions()`
* Improve `image_contrast()` function to avoid error regarding image resolution.
* New argument `subfolder` in `image_export()` to export an image to a subfolder.
* Now `EBImage` installation is checked when pliman is installed.
* `image_pliman()` now returns the image object instead of the path to the image. So, it is not necessarily to call it within `image_import()`.

# pliman 0.3.0
## New functions
* `image_autocrop()` for automatic image cropping.
* `image_filter()` to perform median-based filtering.
* `image_contrast()` to improve contrast by performing adaptive histogram equalization
* `object_coord()` to get the object coordinates and (optionally) draw a bounding rectangle around multiple objects in an image.
* `object_id()` to get the object identification in an image.
* `object_isolate()` to isolate an object from an image.
* `prop_segmented()` to perform (iterative) image segmentation with pixels proportion.

## Minor improvements
* New argument `filter` in `count_objects()` and `prop_segmented()`.

# pliman 0.2.0
* Includes the suggestions given by the CRAN team in the first submission

# pliman 0.1.0

* The first version of `pliman` package submitted to CRAN.
