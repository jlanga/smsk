#!/usr/bin/env Rscript

# Based on https://orfe.princeton.edu/help/r-packages

library_path = Sys.getenv("R_LIBS_USER")

dir.create(
    library_path,
    showWarnings = FALSE,
    recursive = TRUE
)


install.packages(
    pkgs= c(
        "ggplot2"
    ),
    lib = library_path,
    repos = "https://cran.rstudio.com"
)
