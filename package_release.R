# create a github repo
usethis::create_from_github("https://github.com/TiagoOlivoto/pliman.git",
                            destdir = "E:/Desktop")

# Ignore files
usethis::use_build_ignore("logo_pliman.ai")
usethis::use_git_ignore("logo_pliman.ai")

# Build package
devtools::install_github("TiagoOlivoto/pliman")
devtools::install_github("TiagoOlivoto/pliman", build_vignettes = TRUE)
devtools::build_vignettes()
devtools::build_manual()
devtools::build_readme()
devtools::document()
mathjaxr::preview_rd("leaf_area")
setwd("E:/Desktop/pliman")
devtools::check()
devtools::install(quick = TRUE)
devtools::check(run_dont_test = TRUE)
usethis::use_release_issue("1.13.0")
usethis::use_version('minor')
# major
# minor
# patch

#Prepare for release:
#Check current CRAN check results
#Polish NEWS
devtools::build_readme()
urlchecker::url_check()
devtools::check(remote = TRUE, manual = TRUE)
devtools::check_win_devel()
rhub::check_for_cran()
revdepcheck::revdep_check(num_workers = 4)

#Update cran-comments.md
#Review pkgdown reference index for, e.g., missing topics
#Draft blog post

#Submit to CRAN:
usethis::use_version('minor')
devtools::submit_cran()
#Approve email
#Wait for CRAN...

#Accepted
usethis::use_github_release()
usethis::use_dev_version()

#Finish blog post
#Tweet
#Add link to blog post in pkgdown news menu



# Site ####
pkgdown::build_site()
pkgdown::build_news()
pkgdown::build_articles()
pkgdown::preview_site()
pkgdown::deploy_site_github()


# Check status ####
foghorn::cran_incoming(pkg = "metan")
foghorn::summary_cran_results(pkg = "metan")
foghorn::cran_results(pkg = "metan")
foghorn::cran_details(pkg = "metan")


# Install RTools ######
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")


# Restart R #######
Sys.which("make")

# Set environment variables ####
Sys.setenv(PATH = paste(Sys.getenv("PATH"), ";C:\\Program Files\\qpdf-10.0.1\\bin"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), ";C:\\Program Files\\RStudio\\bin\\pandoc"))
Sys.which(Sys.getenv("R_QPDF", "qpdf"))
Sys.getenv("PATH")


# install all packages #####
lib_loc <- "D:/Documents/R/win-library/4.0"
to_install <- unname(installed.packages(lib.loc = lib_loc)[, "Package"])
to_install
install.packages(pkgs = to_install)


# Format description  ####
desc::desc_normalize(file = ".")

# find and replace
xfun::gsub_dir(dir = "R",
               pattern = "pliman",
               replacement  = "pliman")

