## Release summary
This is the first submission of `pliman` (0.1.0) to CRAN

## Test environments
* local R installation, R 4.0.5
* ubuntu 16.04 (on travis-ci), R 4.0.5
* win-builder (devel)

## R CMD check results (run with --run-donttest --as-cran)
0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... Note_to_CRAN_maintainers Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'

   - I think the following note doesn't matter.

* checking top-level files ... NOTE Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
   - I can put these in `.Rbuildignore` but CRAN supports `NEWS.md` now.
