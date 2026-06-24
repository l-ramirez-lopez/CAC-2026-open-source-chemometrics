# FIRST TRAINING CONDUCTED IN AFRICA (NAME THE COUNTRY) AND PUT PICS [NO PUBLICALY AVAILABLE YET]


# TALK ABOUT THE STRUCTURE OF THE SPECTRAL DATA IN R
# spc IS A MATRIX OF SPECTRA ...
# send the PAPER and PUT a Screenshot of it
# for gesearch... the example Leo is a bad szudent and in competitions he alys makes his team lose

library(resemble)
library(proximetricsR)
library(readxl)


local_data <- proxiscout_read_data("data/local-soil-spectra.xlsx", "data/local-soil-properties.xlsx")
local_data_spc <- proxiscout_read_data("data/local-soil-spectra.xlsx")
local_data_properties <- as.data.frame(read_excel("data/local-soil-properties2.xlsx"))


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




##### SPECTRAL LIBRARY ######
spectral_library <- readRDS("data/spectral-library-soil-lucas2009.rds")

# convert the wavelengths to wavenumbers
library_wavn <- 10000000 / as.numeric(colnames(spectral_library$spc))
colnames(spectral_library$spc) <- library_wavn


resample_recipe <- preprocess_recipe(
  prep_resample(grid = "proxiscout"),
  device = "unspecified"
)

spectral_library$spc <- process(
  spectral_library$spc, 
  resample_recipe
)


which(!colnames(local_data$spc) %in% colnames(spectral_library$spc))


local_data$spc <- local_data$spc[, colnames(spectral_library$spc)]

dim(local_data$spc)
dim(spectral_library$spc)




recipe_00 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_derivative(m = 2, w = 11, p = 2, algorithm = "savitzky-golay"),
  prep_snv(),
  device = "proxiscout"
)

recipe_01 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_snv(),
  prep_derivative(m = 1, w = 9, p = 1, algorithm = "savitzky-golay"),
  device = "proxiscout"
)


recipe_02 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_snv(),
  prep_derivative(m = 1, w = 13, p = 1, algorithm = "savitzky-golay"),
  device = "proxiscout"
)

recipe_03 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_derivative(m = 2, w = 7, p = 2, algorithm = "savitzky-golay"),
  prep_detrend(p = 2),
  device = "proxiscout"
)


local_data$Sand <- as.numeric(local_data$SaNA)
local_data$Clay <- as.numeric(local_data$Clay)
my_model <- calibrate_models(
  list(Sand ~ spc),
  data = local_data,
  preprocess = list(recipe_00, recipe_01, recipe_02, recipe_03),
  method = list(fit_plsr(ncomp = 7, type = "modified")),
  control = calibration_control(
    validation_type = "kfold",
    number = 5,
    seed = 42, 
    remove_outliers = 3
  ),
  verbose = FALSE
)

my_model

to_remove_local <- unique(unlist(my_model$final_models[[1]]$detected_outliers$model_1))



spectral_library <- spectral_library[-to_remove_local, ]

## the preprocessing recipe to cundcut the search
my_preprocessing_search <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_derivative(m = 2, w = 15, p = 2, algorithm = "savitzky-golay"),
  prep_snv(),
  device = "proxiscout"
)



ex2 <- search_neighbors(
  Xr = process(spectral_library$spc, my_preprocessing_search), 
  Xu = process(local_data$spc, my_preprocessing_search),
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(40),
    scale = FALSE,
    return_projection = TRUE
  ),
  Yr =  spectral_library$OC,
  neighbors = neighbors_k(50)
)

length(ex2$unique_neighbors)













# Parallel processing
library(doParallel)
n_cores <- 12
cl <- makeCluster(n_cores)
registerDoParallel(cl)

spectral_library <- spectral_library[!is.na(spectral_library$Sand), ]

set.seed(8011)
rsamples <- sample(nrow(local_data), 50)


gs <- gesearch(
  Xr = process(spectral_library$spc, my_preprocessing_search)[ex2$unique_neighbors, ], 
  Yr = spectral_library$Sand[ex2$unique_neighbors],
  Xu = process(local_data$spc, my_preprocessing_search)[rsamples, ], 
  Yu = local_data$Sand[rsamples],
  k = 75, b = 75, retain = 0.95,
  target_size = 120,
  fit_method = fit_pls(ncomp = 15, method = "simpls", scale = FALSE),
  optimization = c("reconstruction", "response"),
  control = gesearch_control(retain_by = "probability"),
  seed = 1410, 
  intermediate_models = TRUE, 
  pchunks = 5
)

stopCluster(cl)
registerDoSEQ()

saveRDS(gs, "gesearch-results-soil-lucas2009.rds")

set.seed(20)
pls_mod <- model(
  Xr = process(local_data$spc, my_preprocessing_search), 
  Yr = local_data$OM,
  fit_method = fit_pls(ncomp = 10, scale = FALSE),
  control = model_control(validation_type = "lgo", number = 10)
)




aa <- predict(
  gs, 
  newdata = process(local_data$spc, my_preprocessing_search)[-rsamples, ], 
  what = "all_generations"
)
cor(aa[[20]][[1]], as.numeric(local_data$Sand[-rsamples]), use = "complete.obs")^2





