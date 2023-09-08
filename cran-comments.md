* **Package**: pliman
* **Version**: 2.0.1
* **Date**: 2023/09/08
* **Dependencies**
   - {pliman} depends on {EBImage} that is not on CRAN. I check this internally and provide a guide to install it. 
* **Reverse dependencies**: None
* **Release summary**
   - This is a minor release (2.0.1) that fixes the error in r-devel <https://cran.r-project.org/web/checks/check_results_pliman.html> (Error in eval(substitute(list(...)), `_data`, parent.frame()): object 'X1' not found).

I just changed "`transform(x, X1 = X1 - mean(X1)...`" with `X1 = x[,1] - mean(x[,1])...`. See the changes at https://github.com/TiagoOlivoto/pliman/commit/c58eb891d80c0632c4c55c8c07516020ef4df4f5#diff-0828b0c0e0293ce78954f344d4d9b805e6f2bda776b29b943866880e0cf2a1f4R454
   
   
* **R CMD check results (run with --run-donttest --as-cran)**
0 errors | 0 warnings | 1 notes
checking CRAN incoming feasibility ... [13s] NOTE 
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>' (I think the following note doesn't matter)

* **devtools::check()**
0 errors | 0 warnings | 0 notes
