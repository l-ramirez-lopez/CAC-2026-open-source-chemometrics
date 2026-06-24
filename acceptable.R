## Packages
library("proximetricsR")
library("tidyverse")
library("curl")
library("qs2") # >=0.25.5

## Separate files
soil_url <- "https://storage.googleapis.com/soilspec4gg-public/neospectra_soillab_v1.2.csv.gz"
site_url <- "https://storage.googleapis.com/soilspec4gg-public/neospectra_soilsite_v1.2.csv.gz"
nir_url  <- "https://storage.googleapis.com/soilspec4gg-public/neospectra_nir_v1.2.csv.gz"
mir_url  <- "https://storage.googleapis.com/soilspec4gg-public/neospectra_mir_v1.2.csv.gz"

soil <- readr::read_csv(soil_url)
# site <- readr::read_csv(site_url)
nir  <- readr::read_csv(nir_url)



spc_lib <- as.matrix(nir[, grep("^scan_nir.", colnames(nir))])
colnames(spc_lib) <- gsub("^scan_nir.|_ref", "", colnames(spc_lib))
colnames(spc_lib) <- 10000000 / as.numeric(colnames(spc_lib))

lib_aggregated_spc <- aggregate(
  spc_lib, 
  by = list(id.sample_local_c = nir$id.sample_local_c),
  FUN = mean
)

lib_spc <- as.matrix(lib_aggregated_spc[, -1])



lib_merged <- merge(
  lib_aggregated_spc,
  soil,
  by.x = "id.sample_local_c",
  by.y = "id.sample_local_c"
)

lib_data <- lib_merged[, -grep("^[0-9]{3}", colnames(lib_merged))]
lib_data$spc <- as.matrix(lib_merged[, grep("^[0-9]{3}", colnames(lib_merged))])

spectral_library <- lib_data
spectral_library$spc <- log10(1 / spectral_library$spc)
colnames(spectral_library$spc) <- 10000000 / as.numeric(colnames(spectral_library$spc))


resample_recipe <- preprocess_recipe(
  prep_resample(grid = "proxiscout"),
  device = "unspecified"
)

spectral_library$spc <- process(
  spectral_library$spc, 
  resample_recipe
)


# TALK ABOUT THE STRUCTURE OF THE SPECTRAL DATA IN R
# spc IS A MATRIX OF SPECTRA ...
# send the PAPER and PUT a Screenshot of it
# for gesearch... the example Leo is a bad szudent and in competitions he alys makes his team lose

library(resemble)
library(proximetricsR)
library(readxl)


local_data <- proxiscout_read_data("data/local-soil-spectra.xlsx", "data/local-soil-properties.xlsx")
local_data_spc <- proxiscout_read_data("data/local-soil-spectra.xlsx")
local_data_properties <- as.data.frame(read_excel("data/local-soil-properties.xlsx"))


local_aggregated_spc <- aggregate(
  local_data_spc$spc, 
  by = list(Group_Name = local_data_spc$Group_Name),
  FUN = mean
)

local_spc <- as.matrix(local_aggregated_spc[, -1])


local_data_properties$SampleName <- as.character(local_data_properties$SampleName)

local_merged <- merge(
  local_aggregated_spc,
  local_data_properties,
  by.x = "Group_Name",
  by.y = "SampleName"
)

local_data <- local_merged[, -grep("^[0-9]{3}", colnames(local_merged))]
local_data$spc <- as.matrix(local_merged[, grep("^[0-9]{3}", colnames(local_merged))])

# convert to absorbance
local_data$spc <- log10(1 / local_data$spc)
local_data$OC <- as.numeric(local_data$OC)




## the preprocessing recipe to cundcut the search
my_preprocessing_search <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_derivative(m = 2, w = 11, p = 2, algorithm = "savitzky-golay"),
  prep_snv(),
  device = "proxiscout"
)

colnames(spectral_library$spc)
colnames(local_data$spc)


pp = process(local_data$spc, my_preprocessing_search)
matplot(
  as.numeric(colnames(pp)),
  t(pp), 
  lty = 1, col = rgb(0, 0, 0, 0.2),
  type = "l"
)




md <- dissimilarity(
  process(local_data$spc, my_preprocessing_search),
  apply(process(local_data$spc, my_preprocessing_search), 2, median),
  diss_method = diss_correlation(ws = 51, center = T, scale = T)
)

threshold <- 0.5

local_data <- local_data[md$dissimilarity < threshold, ]



ex2 <- search_neighbors(
  Xr = process(spectral_library$spc, my_preprocessing_search), 
  Xu = process(local_data$spc, my_preprocessing_search),
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(40),
    scale = FALSE,
    return_projection = TRUE
  ),
  Yr =  spectral_library$oc_usda.c729_w.pct,
  neighbors = neighbors_k(50)
)

length(ex2$unique_neighbors)



stopCluster(cl)
registerDoSEQ()



# Parallel processing
library(doParallel)
n_cores <- 12
cl <- makeCluster(n_cores)
registerDoParallel(cl)

spectral_library <- spectral_library[!is.na(spectral_library$oc_usda.c729_w.pct), ]

set.seed(8011)
rsamples <- sample(nrow(local_data), 30)
local_y <- local_data$OC
local_y[-rsamples] <- NA



gs <- gesearch(
  Xr = process(spectral_library$spc, my_preprocessing_search)[ex2$unique_neighbors, ], 
  Yr = spectral_library$oc_usda.c729_w.pct[ex2$unique_neighbors],
  Xu = process(local_data$spc, my_preprocessing_search), 
  Yu = local_y,
  k = 50, b = 150, retain = 0.97,
  target_size = 50,
  fit_method = fit_pls(ncomp = 10, method = "mpls", scale = FALSE),
  optimization = c("reconstruction", "similarity", "response"),
  control = gesearch_control(retain_by = "probability"),
  seed = 1410, 
  intermediate_models = TRUE, 
  pchunks = 5
)

stopCluster(cl)
registerDoSEQ()

aa <- predict(
  gs, 
  newdata = process(local_data$spc, my_preprocessing_search)[-rsamples, ], 
  what = "all_generations"
)
cor(aa[[34]][[1]], as.numeric(local_data$OC[-rsamples]), use = "complete.obs")^2


# saveRDS(gs, "gesearch-results-soil-lucas2009.rds")

set.seed(20)
pls_mod <- model(
  Xr = process(local_data$spc, my_preprocessing_search), 
  Yr = local_data$OM,
  fit_method = fit_pls(ncomp = 10, scale = FALSE),
  control = model_control(validation_type = "lgo", number = 10)
)

