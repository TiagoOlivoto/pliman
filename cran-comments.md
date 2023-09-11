* **Package**: pliman
* **Version**: 2.0.1
* **Date**: 2023/09/11
* **Dependencies**
   - {pliman} depends on {EBImage} that is not on CRAN. I check this internally and provide a guide to install it. 
* **Reverse dependencies**: None
* **Release summary**
   - This is a minor release (2.0.1) that fixes the Note *Package has a VignetteBuilder field but no prebuilt vignette index.*
   

I removed the `VignetteBuilder` field from DESCRIPTION file since there is no longer a vignete build in the package.
   
   
* **R CMD check results (run with --run-donttest --as-cran)**
0 errors | 0 warnings | 1 notes
checking CRAN incoming feasibility ... [13s] NOTE 
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>' (I think the following note doesn't matter)

* **devtools::check()**
0 errors | 0 warnings | 0 notes
