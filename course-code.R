library(proximetricsR)

mlocal_cal <- proxiscout_read_data(
  "data/local-samples-spain/local-spectra-cal.xlsx", 
  references_file = "data/local-samples-spain/local-properties-cal.xlsx"
)

mlocal_val <- proxiscout_read_data(
  "data/local-samples-spain/local-spectra-val.xlsx", 
  references_file = "data/local-samples-spain/local-properties-val.xlsx"
)

mlibrary <- proxiscout_read_data(
  "data/spectral-library-usa/spectral-library-spectra.xlsx",
  references_file = "data/spectral-library-usa/spectral-library-properties.xlsx"
)
