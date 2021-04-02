
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plaiman <img src="man/figures/logo_plaiman.png" align="right" height="140/"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version-ago/metan)](https://CRAN.R-project.org/package=metan)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
![Downloads](http://cranlogs.r-pkg.org/badges/metan) ![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/metan)
[<img src="https://zenodo.org/badge/130062661.svg" alt="DOI" width="186"/>](https://zenodo.org/badge/latestdoi/130062661)

<!-- badges: end -->

`plaiman` (**pla**nt **im**age **an**anlysis) is designed to analyze
plant images, especially related to leaf analysis. You provide color
palettes, tell `plaiman` what each one represents, and it takes care of
the details. The package will help you to:

-   Measure leaf area
-   Measure disease severity
-   Count the number of lesions

# Installation

The development version of `plaiman` can be installed from
[GitHub](https://github.com/TiagoOlivoto/metan) with:

``` r
devtools::install_github("TiagoOlivoto/metan")

# To build the HTML vignette use
devtools::install_github("TiagoOlivoto/metan", build_vignettes = TRUE)
```

*Note*: If you are a Windows user, you should also first download and
install the latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

For the latest release notes on this development version, see the [NEWS
file](https://tiagoolivoto.github.io/metan/news/index.html).

# Usage

``` r
library(plaiman)
# |===========================================|
# | Tools for Plant Image Analysis (plaiman)  |
# | Author: Tiago Olivoto                     |
# |===========================================|
img <- import_image(system.file("tmp_images/sev2.png", package = "plaiman"))
healthy <- import_image(system.file("tmp_images/sev_healthy.png", package = "plaiman"))
symptoms <- import_image(system.file("tmp_images/sev_sympt.png", package = "plaiman"))
background <- import_image(system.file("tmp_images/sev_back.png", package = "plaiman"))
combine_images(img, healthy, symptoms, background, ncol = 4)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

``` r
symptomatic_area(img = img,
                 img_healthy = healthy,
                 img_symptoms = symptoms,
                 img_background = background,
                 show_image = TRUE)
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

    #    healthy symptomatic
    # 1 88.67028    11.32972

# Getting help

-   If you encounter a clear bug, please file a minimal reproducible
    example on [github](https://github.com/TiagoOlivoto/plaiman/issues)

-   Suggestions and criticisms to improve the quality and usability of
    the package are welcome!
