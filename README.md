
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pliman <img src="man/figures/logo_pliman.svg" align="right" height="140/"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version-ago/pliman)](https://CRAN.R-project.org/package=pliman)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental-1)
![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-week/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-day/pliman?color=orange)](https://r-pkg.org/pkg/pliman)
[![DOI](https://zenodo.org/badge/352844585.svg)](https://zenodo.org/badge/latestdoi/352844585)
<!-- badges: end -->

`pliman` (**pl**ant **im**age **an**alysis) is designed to analyze plant
images, especially related to leaf analysis. You provide color palettes,
tell `pliman` what each one represents, and it takes care of the
details. Image indexes can also be used to segment images. The package
will help you to:

- Measure leaf area;
- Measure disease severity;
- Count the number of lesions;
- Obtain the shape of lesions;
- Produce Standard Area Diagrams;
- Count objects in an image;
- Get object features (area, perimeter, radius, circularity,
  eccentricity, solidity);
- Get the RGB values for each object in an image;
- Get the object coordinates;
- Get the object contours;
- Get convex hulls;
- Isolate objects;
- Plot object measures;
- Perform Elliptical Fourier Analysis;
- 

`pliman` also provides useful functions for image
[transformation](https://tiagoolivoto.github.io/pliman/reference/utils_transform.html),
[binarization](https://tiagoolivoto.github.io/pliman/reference/image_binary.html),
[segmentation](https://tiagoolivoto.github.io/pliman/reference/image_segment.html),
and
[resolution](https://tiagoolivoto.github.io/pliman/reference/utils_dpi.html).
Please, visit the
[Examples](https://tiagoolivoto.github.io/pliman/index.html) page in
`pliman` website for a detailed documentation of each function.

# Installation

Install the latest stable version of `pliman` from
[CRAN](https://CRAN.R-project.org/package=pliman) with:

``` r
install.packages("pliman")
```

The development version of `pliman` can be installed from
[GitHub](https://github.com/TiagoOlivoto/pliman) with:

``` r
devtools::install_github("TiagoOlivoto/pliman")

# To build the HTML vignette use
devtools::install_github("TiagoOlivoto/pliman", build_vignettes = TRUE)
```

*Note*: If you are a Windows user, you should also first download and
install the latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

# Basic usage

The two main functions of pliman are `analyze objects()` and
`measure disease()`, which can be used to compute object features and
measure disease severity, respectively.

## Analyze objects

The function `analyze_objects()` can be used to analyze objects such as
leaves, grains, pods, and pollen in an image. The following example
counts and computes several features of soybean grains of an image with
30 grains.

``` r
library(pliman)
# |==========================================================|
# | Tools for Plant Image Analysis (pliman 1.2.0)            |
# | Author: Tiago Olivoto                                    |
# | Type 'citation('pliman')' to know how to cite pliman     |
# | Type 'vignette('pliman_start')' for a short tutorial     |
# | Visit 'http://bit.ly/pkg_pliman' for a complete tutorial |
# |==========================================================|
img <-image_pliman("soybean_touch.jpg", plot = TRUE)
soy <- analyze_objects(img, marker = "id")
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
str(soy$results)
# 'data.frame': 30 obs. of  47 variables:
#  $ id                  : num  1 2 3 4 5 6 7 8 9 10 ...
#  $ x                   : num  246 538 238 345 278 ...
#  $ y                   : num  510 402 340 106 261 ...
#  $ area                : num  2275 2299 2311 2443 2161 ...
#  $ area_ch             : num  2302 2270 2289 2413 2123 ...
#  $ perimeter           : num  184 178 181 186 172 ...
#  $ radius_mean         : num  26.5 26.6 26.7 27.5 25.8 ...
#  $ radius_min          : num  22.9 25 24 24.3 24.2 ...
#  $ radius_max          : num  29.4 28.7 29.4 30.5 28 ...
#  $ radius_sd           : num  1.369 0.945 1.24 1.757 0.803 ...
#  $ diam_mean           : num  52.9 53.2 53.4 55 51.5 ...
#  $ diam_min            : num  45.9 49.9 48 48.7 48.5 ...
#  $ diam_max            : num  58.9 57.4 58.9 61.1 56.1 ...
#  $ major_axis          : num  56.1 56.5 57.5 60.9 54.1 ...
#  $ minor_axis          : num  51.9 51.9 51.3 51.1 50.9 ...
#  $ caliper             : num  57.3 56.9 57.7 61 54.4 ...
#  $ length              : num  56.6 56.5 57.2 61 54 ...
#  $ width               : num  51.6 52.4 52 51 50.6 ...
#  $ radius_ratio        : num  1.28 1.15 1.23 1.25 1.16 ...
#  $ theta               : num  -0.889 -0.841 -0.565 -0.993 -0.218 ...
#  $ eccentricity        : num  0.893 0.853 0.816 0.754 0.896 ...
#  $ form_factor         : num  0.849 0.91 0.886 0.891 0.921 ...
#  $ narrow_factor       : num  1.01 1.01 1.01 1 1.01 ...
#  $ asp_ratio           : num  1.1 1.08 1.1 1.2 1.07 ...
#  $ rectangularity      : num  1.28 1.29 1.29 1.27 1.26 ...
#  $ pd_ratio            : num  3.2 3.13 3.14 3.04 3.16 ...
#  $ plw_ratio           : num  1.7 1.64 1.66 1.66 1.64 ...
#  $ solidity            : num  0.988 1.013 1.01 1.012 1.018 ...
#  $ convexity           : num  0.886 0.881 0.911 0.915 0.898 ...
#  $ elongation          : num  0.0884 0.0728 0.0911 0.1638 0.0642 ...
#  $ circularity         : num  14.8 13.8 14.2 14.1 13.6 ...
#  $ circularity_haralick: num  19.3 28.1 21.5 15.7 32.1 ...
#  $ circularity_norm    : num  1.22 1.14 1.17 1.16 1.13 ...
#  $ coverage            : num  0.00426 0.0043 0.00432 0.00457 0.00404 ...
#  $ asm                 : num  0.0679 0.0886 0.0848 0.0869 0.0865 ...
#  $ con                 : num  0.527 0.603 0.702 0.73 0.834 ...
#  $ cor                 : num  0.961 0.959 0.95 0.945 0.944 ...
#  $ var                 : num  7.68 8.37 7.98 7.58 8.46 ...
#  $ idm                 : num  0.81 0.807 0.794 0.801 0.795 ...
#  $ sav                 : num  40.5 40.7 40.7 40.6 40.9 ...
#  $ sva                 : num  1566 1591 1590 1583 1606 ...
#  $ sen                 : num  1.24 1.18 1.17 1.17 1.2 ...
#  $ ent                 : num  1.4 1.36 1.36 1.36 1.4 ...
#  $ dva                 : num  0.527 0.603 0.702 0.73 0.834 ...
#  $ den                 : num  0.347 0.364 0.386 0.385 0.392 ...
#  $ f12                 : num  0.569 0.56 0.535 0.544 0.542 ...
#  $ f13                 : num  0.82 0.807 0.793 0.798 0.804 ...
```

## Disease severity

``` r
img <-image_pliman("sev_leaf.jpg")
healthy <-image_pliman("sev_healthy.jpg")
symptoms <-image_pliman("sev_sympt.jpg")
background <-image_pliman("sev_back.jpg")
image_combine(img, healthy, symptoms, background, ncol = 4)
```

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

``` r
sev <- 
  measure_disease(img = img,
                  img_healthy = healthy,
                  img_symptoms = symptoms,
                  img_background = background)
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

``` r

sev$severity
#    healthy symptomatic
# 1 89.20609    10.79391
```

`pliman` takes the advantage of several powerful functions from [EBImage
package](https://bioconductor.org/packages/release/bioc/html/EBImage.html).
Thanks to Andrzej Oleś and collaborators for the impressive job done
with EBImage!

# Citation

``` r
citation("pliman")

Please, support this project by citing it in your publications!

  Olivoto, T.(2022). Lights, camera, pliman! An R package for plant
  image analysis. Methods Ecol Evol. 13:789-798
  doi:10.1111/2041-210X.13803

A BibTeX entry for LaTeX users is

  @Article{Olivoto2022,
    author = {Tiago Olivoto},
    title = {Lights, camera, pliman! An R package for plant image analysis},
    journal = {Methods in Ecology and Evolution},
    volume = {13},
    number = {4},
    pages = {789-798},
    year = {2022},
    doi = {10.1111/2041-210X.13803},
  }
```

# Getting help

- If you encounter a clear bug, please file a minimal reproducible
  example on [github](https://github.com/TiagoOlivoto/pliman/issues).
  The package [reprex](https://reprex.tidyverse.org/) may help you with
  that.

- Suggestions and criticisms to improve the quality and usability of the
  package are welcome!

# Code of Conduct

Please note that the pliman project is released with a [Contributor Code
of Conduct](https://tiagoolivoto.github.io/pliman/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

<div align="center">

<a href='https://www.free-website-hit-counter.com'><img src='https://www.free-website-hit-counter.com/c.php?d=9&id=144207&s=2' border='0' alt='Free Website Hit Counter'></a><br/><small><a href='https://www.free-website-hit-counter.com' title="Free Website Hit Counter">Free
website hit counter</a></small>

</div>
