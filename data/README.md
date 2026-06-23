# Datasets

This directory holds the spectral and reference datasets used throughout the course.

## Datasets included

### `NIRsoil` (bundled with `prospectr`)

The `NIRsoil` dataset is shipped with the **prospectr** package and does not need
to be downloaded separately. It contains diffuse reflectance NIR spectra (1100–2498 nm,
2 nm steps) of 825 soil samples along with reference measurements for:

- `Nt`  — total nitrogen (%)
- `Ciso` — isotopic carbon (‰)
- `CEC` — cation exchange capacity (cmol/kg)
- `moisture` — gravimetric soil moisture (%)

Access it in R with:

```r
library(prospectr)
data(NIRsoil)
```

### `shootout` (bundled with `pls`)

The classic NIR shootout dataset (2002 IDRC Shootout) is included in the **pls**
package. It contains NIR spectra of pharmaceutical tablets and reference tablet
hardness values.

```r
library(pls)
data(NIR)
```

### Bring your own data

If you have spectra from your own handheld NIR instrument, you can load them
using `read.csv()` or `readr::read_csv()`. The expected format is:

| sample_id | 1000 | 1002 | … | 2500 | reference_value |
|-----------|------|------|---|------|-----------------|
| S001      | 0.41 | 0.42 | … | 0.55 | 12.3            |

Place your CSV or Excel files in this directory and adjust the file paths in the
notebooks accordingly.

## Data provenance

| Dataset | Source | Licence |
|---------|--------|---------|
| `NIRsoil` | Ramirez-Lopez et al. (2013), bundled in `prospectr` | GPL-2 |
| `NIR` (shootout) | 2002 IDRC Shootout, bundled in `pls` | GPL-2 |
