* **Package**: pliman
* **Version**: 2.1.0
* **Date**: 2023/10/14
* **Dependencies**
   - {pliman} depends on {EBImage} that is not on CRAN. I check this internally and provide a guide to install it. 
* **Reverse dependencies**: None
* **Release summary**
   - This is a path release (2.1.0) that fixes some bugs and includes new functionalities. In this version, I removed the dependency on the raster package in favor of more recent packages, namely, terra, stars, and sf.

   
* **R CMD check results (run with --run-donttest --as-cran)**
0 errors | 0 warnings | 1 notes
checking CRAN incoming feasibility ... [13s] NOTE 
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>' (I think the following note doesn't matter)

* **devtools::check()**
0 errors | 0 warnings | 0 notes
