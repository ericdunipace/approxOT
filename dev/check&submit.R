# runs checks 
usethis::use_release_issue("1.1")

# cran comments
usethis::use_cran_comments()

# url checks
urlchecker::url_check()

# command check local
devtools::check(remote = TRUE, manual = TRUE)

# command checks, remote. Depends on RcppCGAL being on CRAN
devtools::check_win_devel(quiet=TRUE)
devtools::check_win_release(quiet = TRUE)
devtools::check_win_oldrelease(quiet = TRUE)

out <- rhub::check_for_cran(show_status = FALSE)

# reverse dependency
# run if no rev dep check devtools::install_github('r-lib/revdepcheck')
revdepcheck::revdep_reset()
revdepcheck::revdep_check(num_workers = 4,
                          timeout = as.difftime(240, units = "mins"))


# update cran comments
usethis::use_cran_comments()

# to submit
devtools::submit_cran()
