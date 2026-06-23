# install_packages.R
# Run this script once before the course to install all required packages.
# Usage: source("install_packages.R")

required_packages <- c(
  # Spectroscopy and chemometrics
  "prospectr",  # spectral pre-processing and sample selection
  "resemble",   # memory-based learning and local models
  "pls",        # PLS and PCR regression
  # Machine learning utilities
  "caret",
  # Data wrangling and visualisation
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "patchwork",  # combining ggplot2 panels
  # Reporting
  "rmarkdown",
  "knitr"
)

# Install any packages that are not yet present
to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(to_install) > 0) {
  message("Installing missing packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}

# Verify successful installation
loaded <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (all(loaded)) {
  message("All packages loaded successfully. You are ready for the course!")
} else {
  warning(
    "The following packages could not be loaded: ",
    paste(names(loaded)[!loaded], collapse = ", ")
  )
}
