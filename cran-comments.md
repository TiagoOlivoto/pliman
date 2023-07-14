* **Package**: pliman
* **Version**: 2.0.1
* **Date**: 2023/07/14
* **Dependencies**
   - Pliman depends on {EBImage} that is not on CRAN. I check this internally and provide a guide to install it. 
* **Reverse dependencies**: None
* **Release summary**
   - This is a major release (2.0.0) that includes several new features and breaking changes. Please, see `NEWS.md` for additional information.
   
   In this resubmission, I removed `CXX_STD=CXX11` from `src/Makevars` and `src/Makevars.win` to fix the note "checking C++ specification ... NOTE - Specified C++11: please drop specification unless essential"
   
* **R CMD check results (run with --run-donttest --as-cran)**
0 errors | 0 warnings | 1 notes
checking CRAN incoming feasibility ... [13s] NOTE 
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>' (I think the following note doesn't matter)

* **devtools::check()**
0 errors | 0 warnings | 0 notes
