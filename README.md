# CAC 2026 — Open-source chemometrics for real-world NIR handheld spectroscopy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Course overview

This repository contains all materials for the **CAC 2026** short course:

> **Open-source chemometrics for real-world NIR handheld spectroscopy**

The course is taught entirely in **R** and covers the full analytical workflow — from raw spectral data collected with a handheld NIR instrument through to validated calibration models ready for field deployment.

---

## Learning objectives

By the end of this course participants will be able to:

1. Load, inspect and visualise NIR spectra from handheld instruments in R.
2. Apply standard spectral pre-processing methods (SNV, MSC, Savitzky–Golay derivatives).
3. Use principal component analysis (PCA) for exploratory analysis and outlier detection.
4. Build and optimise partial least squares (PLS) regression calibration models.
5. Assess model performance through cross-validation and independent test-set validation.
6. Apply memory-based / local modelling strategies suited to heterogeneous field data.
7. Report and communicate results using reproducible R Markdown documents.

---

## Prerequisites

| Requirement | Details |
|---|---|
| R | ≥ 4.2 |
| RStudio (recommended) | ≥ 2023.03 |
| Basic R knowledge | Familiarity with data frames, basic plots and functions |
| NIR data | Either bring your own or use the example datasets provided |

---

## Repository structure

```
CAC-2026-open-source-chemometrics/
├── README.md                        # This file
├── CAC-2026-chemometrics.Rproj      # RStudio project file
├── install_packages.R               # One-time package installation script
├── data/
│   └── README.md                    # Dataset descriptions and provenance
└── notebooks/
    ├── 01_introduction.Rmd          # R environment and NIR spectroscopy basics
    ├── 02_data_loading.Rmd          # Loading and visualising spectra
    ├── 03_preprocessing.Rmd         # Spectral pre-processing
    ├── 04_exploratory_analysis.Rmd  # PCA and outlier detection
    ├── 05_pls_calibration.Rmd       # PLS calibration model building
    └── 06_model_validation.Rmd      # Cross-validation and external validation
```

---

## Getting started

### 1 — Clone or download the repository

```bash
git clone https://github.com/l-ramirez-lopez/CAC-2026-open-source-chemometrics.git
```

Or download the ZIP archive from the green **Code** button on GitHub.

### 2 — Open the RStudio project

Double-click `CAC-2026-chemometrics.Rproj` to open the project in RStudio.  
Working inside the project keeps file paths consistent across machines.

### 3 — Install required packages

Run the installation script once before the course:

```r
source("install_packages.R")
```

### 4 — Work through the notebooks in order

Open each `.Rmd` file in `notebooks/` in RStudio and click **Knit** (or run
chunks interactively with Ctrl+Enter / Cmd+Enter).

---

## Key R packages used

| Package | Purpose |
|---|---|
| [`prospectr`](https://CRAN.R-project.org/package=prospectr) | Spectral pre-processing and sample selection |
| [`resemble`](https://CRAN.R-project.org/package=resemble) | Memory-based learning and local PLS models |
| [`pls`](https://CRAN.R-project.org/package=pls) | PLS and PCR regression |
| [`ggplot2`](https://CRAN.R-project.org/package=ggplot2) | Publication-quality graphics |
| [`dplyr`](https://CRAN.R-project.org/package=dplyr) | Data manipulation |
| [`tidyr`](https://CRAN.R-project.org/package=tidyr) | Data reshaping |
| [`caret`](https://CRAN.R-project.org/package=caret) | Model training utilities |

---

## Instructor

**Leonardo Ramirez-Lopez**  
Author of the `prospectr` and `resemble` R packages.  
Research in applied spectroscopy and chemometrics for agriculture, food and environment.

---

## License

This course material is released under the [MIT License](LICENSE).  
You are free to reuse and adapt it with attribution.
