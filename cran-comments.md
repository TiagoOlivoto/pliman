## Release summary
This is a resubmission of the `pliman` (0.1.0) that includes the suggestions given by Julia Haider in the first submission of the package to CRAN.

## Suggestions

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* Implemented. Please, see the DESCRIPTION file.


Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      image_pallete.Rd: \value
      utils_file.Rd: \value
      
* Done! I have included \value to all the .Rd files.

I hope that fixing these two issues makes the package suitable for publication. If any changes still need to be made, please, let me known.
Thank you for the time in reviewing this new submission

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
