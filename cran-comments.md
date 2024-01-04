* **Package**: pliman
* **Version**: 2.1.0
* **Date**: 2023/10/14
* **Dependencies**
   - {pliman} depends on {EBImage} that is not on CRAN. I check this internally and provide a guide to install it. 
* **Reverse dependencies**: None
* **Release summary**
   - This is a path release (2.2.0) that fixes some bugs and includes the `mosaic_*()` and `shapefile_*()` family of functions. In this version, I removed the dependency on stars package.

   
* **R CMD check results (run with --run-donttest --as-cran)**
0 errors | 0 warnings | 2 notes


```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.


```
* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''
```

As noted in [R-hub issue #560](https://github.com/r-hub/rhub/issues/560), this seems to be an Rhub issue and so can likely be ignored. 



* **devtools::check()**
0 errors | 0 warnings | 0 notes
