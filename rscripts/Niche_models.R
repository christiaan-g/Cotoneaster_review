################################################################################
#
# Cotoneaster: climatic niche comparison and model transferability
#
# Authors : Luke Potgieter & Christiaan Gildenhuys
# Affiliation: Centre for Invasion Biology (C·I·B), Department of Botany
#              and Zoology, Stellenbosch University, South Africa
# Year    : 2026
#
# Associated manuscript:
#   Potgieter LJ & Gildenhuys C. A global review of Cotoneaster (Rosaceae)
#   invasions: distribution, drivers, and projected spread.
#   Biological Invasions (submitted 2026).
#
# Description:
#   Runs end-to-end from the cleaned GBIF occurrence file produced by the
#   main pipeline (gbif_final_species.rds). Performs climatic niche
#   comparisons (COUE), native-range SDMs with spatially blocked cross-
#   validation, geographic transfer evaluation against geographically
#   constrained null models, MESS extrapolation diagnostics, and a
#   reciprocal transfer test (naturalised -> native).
#
# Focal species: Cotoneaster acutifolius, C. divaricatus, C. horizontalis
#   (other taxa excluded due to insufficient native-range records for
#    reliable model calibration)
#
# Data requirements (not included in this repository):
#   1. gbif_final_species.rds  -- cleaned occurrence data from main pipeline
#   2. CHELSA V2.1 bioclimatic layers (.tif)
#      https://chelsa-climate.org/bioclim/
#   3. Natural Earth 50 m land polygon shapefile (ne_50m_land.zip)
#      https://www.naturalearthdata.com/downloads/50m-physical-vectors/
#
# R package dependencies:
#   tidyverse, sf, terra, ade4, ecospat, usdm, maxnet, blockCV,
#   precrec, predicts, rnaturalearth, patchwork, tidytext
#   flexsdm  (used via namespace: flexsdm::calib_area — install separately)
#
# CONTENTS
#   0.  Setup, paths, packages
#   1.  Occurrence preparation and spatial thinning (+ sensitivity)
#   2.  Climate predictors: load, mask, collinearity reduction
#   3.  Accessible areas and background sampling
#   4.  COUE niche comparison (+ intersection sweep, equivalency, similarity)
#   5.  Native-only models: fit, project, threshold
#   6.  Model evaluation: blocked CV, geographic transfer, binary performance
#   7.  Null model comparison
#   8.  Regional breakdown of transfer performance
#   9.  Extrapolation diagnostics (MESS, aggregate and per variable)
#   10. Response curves
#   11. Reciprocal test (naturalised -> native), full and sample-matched
#   12. Figures
#   13. Tables
#   14. Save objects and session info
#
# NOTES ON DECISIONS MADE (see manuscript Methods):
#   - Thinning on an equal-area grid, not the raster grid, so that native
#     (c. 30°N) and naturalised (c. 50°N) records are filtered at the same
#     true distance. 5 km used; sensitivity at 1/2/5/10 km reported.
#   - CHELSA bio layers carry interpolated values over ocean and MUST be
#     masked to land before any extraction or background sampling.
#   - bio6 substituted for bio1 and bio4 on a priori ecological grounds
#     (not a VIF result; reported as such in the manuscript).
#   - Feature classes and regularisation held constant across species so
#     that performance differences are not confounded with model configuration.
#   - Null models evaluated identically to observed models (restricted to the
#     naturalised accessible area); evaluating nulls globally inflates them.
#
################################################################################

################################################################################



################################################################################

# 0. SETUP

################################################################################


library(tidyverse)
library(sf)
library(terra)
library(ade4)
library(ecospat)
library(usdm)
library(maxnet)
library(blockCV)
library(precrec)
library(predicts)
library(rnaturalearth)
library(patchwork)
library(tidytext)


set.seed(20260721)


# ---- Paths
#
# Set the three paths below to match your local file system.
# All output sub-directories (figures/, tables/) are created automatically.
#
# chelsa_dir must contain the CHELSA V2.1 bio*.tif files.
# land_shp   must point to the ne_50m_land.shp file from Natural Earth.
# out_main   must already contain gbif_final_species.rds (from main pipeline).
#            -------------------------------------------------------------------

out_main   <- "path/to/output"          # <-- SET THIS
chelsa_dir <- "path/to/CHELSA_V2.1"     # <-- SET THIS

land_shp   <- "path/to/ne_50m_land/ne_50m_land.shp"  # <-- SET THIS


fig_dir <- file.path(out_main, "figures")
tab_dir <- file.path(out_main, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Analysis constants ------------------------------------------------------

focal_species <- c("Cotoneaster acutifolius","Cotoneaster divaricatus", "Cotoneaster horizontalis")


THIN_KM <- 5 # occurrence filtering distance
BUFFER_KM <- 200 # accessible area buffer
N_BG <- 10000 # background points per range
AGG_FACT <- 20 # aggregation factor for global projection (~20 km)
FC <- "lq" # maxnet feature classes
REG <- 2 # maxnet regularisation multiplier
N_NULL <- 99 # null model replicates
N_MATCH <- 50 # sample-matched reciprocal replicates
N_PERM <- 1000 # equivalency / similarity permutations


# Treatment of records disputed between WCVP (native) and GloNAF
(alien):

# "drop" | "native" | "naturalised"

DISPUTED_TREATMENT <- "drop"


# ---- Helper: maxnet prediction wrapper for terra::predict --------------------

# maxnet's clamping code uses matrix-subscript assignment, which tibbles reject,
# so newdata must always be a plain data.frame. predict.maxnet also returns an
# n x 1 matrix, hence as.vector where a vector is needed.


maxnet_fun <- function(model, data, ...) {
  predict(model, newdata = data, type = "cloglog", clamp = TRUE)
}

maxnet_fun_noclamp <- function(model, data, ...) {
  predict(model, newdata = data, type = "cloglog", clamp = FALSE)
}


pred_pts <- function(model, df, clamp = TRUE) {
  as.vector(predict(model, newdata = as.data.frame(df),
                    type = "cloglog", clamp = clamp))
}


################################################################################

# 1. OCCURRENCE PREPARATION AND THINNING

################################################################################


gbif_final_sp <- readRDS(file.path(out_main, "gbif_final_species.rds"))

occ_all <-
  gbif_final_sp %>%
  filter(wcvp_accepted_name %in% focal_species,
         !is.na(longitude), !is.na(latitude),
         nativeness %in% c("native", "nonnative", "disputed")) %>%
  mutate(range = case_when(
    nativeness == "native" ~ "native",
    nativeness == "nonnative" ~ "naturalised",
    nativeness == "disputed" & DISPUTED_TREATMENT == "native" ~ "native",
    nativeness == "disputed" & DISPUTED_TREATMENT == "naturalised" ~
      "naturalised",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(range)) %>%
  select(species = wcvp_accepted_name, range, longitude, latitude, year)


occ_all %>% count(species, range)


# ---- Equal-area thinning -----------------------------------------------------

# Projecting to cylindrical equal-area before binning means one record per
# d x d km cell everywhere, rather than one per d x d degree cell, which would
# filter the Chinese native records far more loosely than the European ones.

thin_grid <- function(df, d_km) {
  p <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform("+proj=cea +lat_ts=30 +datum=WGS84")
  xy <- st_coordinates(p)
  key <- paste(floor(xy[, 1] / (d_km * 1000)),
               floor(xy[, 2] / (d_km * 1000)))
  df[!duplicated(key), ]
}


# ---- Sensitivity to filtering distance ---------------------------------------

thin_sensitivity <-
  occ_all %>%
  group_by(species, range) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(.x$species[1], .x$range[1], sep = "_")))%>%
  
  map_dfr(function(df) {
    
    tibble(d_km = c(1, 2, 5, 10)) %>%
      mutate(n = map_int(d_km, ~ nrow(thin_grid(df, .x))))
    
  }, .id = "set") %>%
  
  pivot_wider(names_from = d_km, values_from = n, names_prefix = "km_")


thin_sensitivity

write_csv(thin_sensitivity, file.path(tab_dir,"tableS2_thinning_sensitivity.csv"))


# ---- Apply chosen distance ---------------------------------------------------

  occ_thin <-
  occ_all %>%
  group_by(species, range) %>%
  group_modify(~ thin_grid(.x, THIN_KM)) %>%
  ungroup()


occ_thin %>% count(species, range)


saveRDS(occ_thin, file.path(out_main, paste0("occ_thinned_", THIN_KM, "km.rds")))


################################################################################

# 2. CLIMATE PREDICTORS

################################################################################


# ---- Load 19 CHELSA bioclim layers -------------------------------------------

bio_files <- list.files(chelsa_dir,
                        pattern = "^CHELSA_bio\\d+_1981-2010_V\\.2\\.1\\.tif$",
                        full.names = TRUE)

stopifnot(length(bio_files) == 19)


predictors <- rast(bio_files)

names(predictors) <- paste0("bio", str_extract(basename(bio_files), "(?<=bio)\\d+"))

predictors <- predictors[[order(as.numeric(str_remove(names(predictors), "bio")))]]

names(predictors)


# Confirm scaling: terra applies CHELSA's embedded scale/offset, so bio1 should
# be in degrees C (roughly -55 to +35). If it returns 0-3500, apply
# value * 0.1 - 273.15 to the temperature layers.

global(predictors[["bio1"]], "range", na.rm = TRUE)


# ---- Mask to land ------------------------------------------------------------

# CRITICAL: CHELSA bio layers carry interpolated values over the ocean. Without
# masking, background points fall in the sea and the presence/background
# contrast becomes partly land-vs-water.

land <- st_read(land_shp, quiet = TRUE)

land_v <- vect(land)

predictors <- mask(predictors, land_v)

global(is.na(predictors[[1]]), "sum", na.rm = FALSE) # should now be large

# ---- Collinearity reduction --------------------------------------------------

occ_vect <- vect(occ_thin, geom = c("longitude", "latitude"), crs = "EPSG:4326")

env_at_occ <- terra::extract(predictors, occ_vect, ID = FALSE) %>%
  drop_na()

vif_all <- usdm::vifstep(env_at_occ, th = 5)

vif_all

# Stepwise VIF retains: bio1 bio3 bio4 bio8 bio9 bio15 bio18 bio19
# No cold-extreme variable survives. bio6 is substituted for bio1 and bio4 on
# a priori grounds (winter cold expected to limit temperate montane shrubs).
# This is an ecological decision, NOT a VIF result - report it as such.


keep_vars <- c("bio6", "bio3", "bio8", "bio9", "bio15", "bio18", "bio19")

usdm::vif(env_at_occ[, keep_vars]) # confirm all < 5


predictors <- predictors[[keep_vars]]

names(predictors)


# ---- Coarse layers for global projection -------------------------------------

# Fitting uses full-resolution extractions; only mapping is coarsened.

pred_coarse <- terra::aggregate(predictors, fact = AGG_FACT, fun = "mean",
                                na.rm = TRUE)

pred_coarse <- mask(pred_coarse, land_v)

res(pred_coarse)


################################################################################

# 3. ACCESSIBLE AREAS AND BACKGROUND

################################################################################

make_ca <- function(df, buffer_km = BUFFER_KM) {
  flexsdm::calib_area(data = df, x = "longitude", y = "latitude",
                      method = c("buffer", width = buffer_km * 1000),
                      crs = "EPSG:4326")
}


calib_areas <-
  occ_thin %>%
  group_by(species, range) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(.x$species[1], .x$range[1], sep = "_")))%>%
  map(make_ca)

map_dbl(calib_areas, ~ expanse(.x, unit = "km")) # report in Methods


# ---- Background sampling -----------------------------------------------------

# Random within each accessible area. Sampling the polygon directly avoids
# masking a global raster (which exhausts memory for the naturalised areas);
# points falling off-raster are discarded and redrawn.


draw_background <- function(ca, n = N_BG, oversample = 1.6) {
  
  repeat {
    pts <- terra::spatSample(ca, size = ceiling(n * oversample), method =
                               "random")
    xy <- terra::crds(pts) %>% as_tibble() %>%
      rename(longitude = x, latitude = y)
    
    ok <- terra::extract(predictors[[1]],
                         vect(xy, geom = c("longitude", "latitude"),
                              crs = "EPSG:4326"), ID = FALSE)[[1]]
    xy <- xy[!is.na(ok), ]
    
    if (nrow(xy) >= n) return(xy[sample(nrow(xy), n), ])
    oversample <- oversample * 1.5
  }
}


backgrounds <- map(calib_areas, draw_background)

map_int(backgrounds, nrow)

# Verify none are marine (should all be zero)

map_int(backgrounds, function(b) {
  v <- terra::extract(predictors[[1]],
                      vect(b, geom = c("longitude", "latitude"),
                           crs = "EPSG:4326"), ID = FALSE)[[1]]
  sum(is.na(v))
})

################################################################################

# 4. COUE NICHE COMPARISON

################################################################################

ex <- function(d) {
  terra::extract(predictors,
                 vect(d, geom = c("longitude", "latitude"), crs = "EPSG:4326"),
                 ID = FALSE) %>% drop_na()
}

run_coue <- function(sp, R = 100, n_reps = 100) {
  
  nat <- occ_thin %>% filter(species == sp, range == "native")
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  
  bg_nat <- backgrounds[[paste(sp, "native", sep = "_")]]
  bg_inv <- backgrounds[[paste(sp, "naturalised", sep = "_")]]
  
  env_nat <- ex(nat)
  env_inv <- ex(inv)
  env_bg_nat <- ex(bg_nat)
  env_bg_inv <- ex(bg_inv)
  
  # PCA calibrated on background AND occurrences of both ranges. Occurrences
  # must be included or a record outside the background's range yields NA
  # density and ecospat.grid.clim.dyn fails.
  
  pca_cal <- dudi.pca(rbind(env_bg_nat, env_bg_inv, env_nat, env_inv),
                      center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
  
  scores_bg_nat <- suprow(pca_cal, env_bg_nat)$li
  scores_bg_inv <- suprow(pca_cal, env_bg_inv)$li
  scores_nat <- suprow(pca_cal, env_nat)$li
  scores_inv <- suprow(pca_cal, env_inv)$li
  scores_glob <- rbind(scores_bg_nat, scores_bg_inv, scores_nat, scores_inv)
  
  grid_nat <- ecospat.grid.clim.dyn(scores_glob, scores_bg_nat, scores_nat, R = R)
  grid_inv <- ecospat.grid.clim.dyn(scores_glob, scores_bg_inv, scores_inv, R = R)
  
  overlap <- ecospat.niche.overlap(grid_nat, grid_inv, cor = TRUE)
  eq <- ecospat.niche.equivalency.test(grid_nat, grid_inv, rep = n_reps, overlap.alternative = "lower")
  sim <- ecospat.niche.similarity.test(grid_nat, grid_inv, rep = n_reps, overlap.alternative = "higher", rand.type = 2)
  dyn <- ecospat.niche.dyn.index(grid_nat, grid_inv, intersection = 0.1)
  
  list(species = sp, pca = pca_cal,
       grid_native = grid_nat, grid_naturalised = grid_inv,
       overlap = overlap, equivalency = eq, similarity = sim,
       dynamic = dyn$dynamic.index.w,
       n_native = nrow(env_nat), n_naturalised = nrow(env_inv))
}

coue <- map(focal_species, run_coue) %>% set_names(focal_species)


# ---- Equivalency and similarity at full replication --------------------------

# These permute the existing grids, so they are far cheaper than refitting.

rerun_tests <- function(res, n_reps = N_PERM) {
  
  list(
    equivalency = ecospat.niche.equivalency.test(res$grid_native,
                                                 res$grid_naturalised,
                                                 rep = n_reps,
                                                 overlap.alternative = "lower"),
    
    similarity = ecospat.niche.similarity.test(res$grid_native,
                                               res$grid_naturalised,
                                               rep = n_reps,
                                               overlap.alternative = "higher",
                                               rand.type = 2)
  )
}


tests_full <- map(coue, rerun_tests)

for (sp in focal_species) {
  coue[[sp]]$equivalency <- tests_full[[sp]]$equivalency
  coue[[sp]]$similarity <- tests_full[[sp]]$similarity
}


map_dbl(coue, ~ .x$overlap$D)
map_dbl(coue, ~ .x$equivalency$p.D) # report 0.000999 as p < 0.001
map_dbl(coue, ~ .x$similarity$p.D)


# ---- Intersection sweep ------------------------------------------------------

# Expansion is highly sensitive to how shared climate space is delimited.

sweep_dyn <- function(res, sp) {
  map_dfr(
    c(0, 0.05, 0.1, 0.15, 0.2, 0.25),
    ~ ecospat.niche.dyn.index(res$grid_native, res$grid_naturalised, intersection = .x)$dynamic.index.w %>%
      as_tibble_row() %>% 
      mutate(intersection = .x)
  ) %>%
    mutate(species = sp) %>%
    select(species, intersection, expansion, stability, unfilling)
}


coue_sweep <- map_dfr(focal_species, ~ sweep_dyn(coue[[.x]], .x))

coue_sweep

write_csv(coue_sweep, file.path(tab_dir,"tableS3_intersection_sweep.csv"))

# ---- Axis interpretation -----------------------------------------------------

axis_var <- map_dfr(focal_species, function(sp) {
  e <- coue[[sp]]$pca$eig
  tibble(species = sp, pc1 = 100 * e[1] / sum(e), pc2 = 100 * e[2] / sum(e))
})


axis_meaning <- map_dfr(focal_species, function(sp) {
  co <- coue[[sp]]$pca$co
  top <- function(k) paste(rownames(co)[order(-abs(co[, k]))][1:2],
                           collapse = ", ")
  tibble(species = sp, pc1_top = top(1), pc2_top = top(2))
})


pca_loadings <- map_dfr(focal_species, function(sp) {
  coue[[sp]]$pca$co %>%
    as_tibble(rownames = "variable") %>%
    rename(PC1 = Comp1, PC2 = Comp2) %>%
    mutate(species = sp, .before = 1)
})


axis_var; axis_meaning

write_csv(pca_loadings, file.path(tab_dir, "tableS4_pca_loadings.csv"))


################################################################################

# 5. NATIVE-ONLY MODELS

################################################################################


fit_native_only <- function(sp) {
  pres <- occ_thin %>% filter(species == sp, range == "native")
  bg <- backgrounds[[paste(sp, "native", sep = "_")]]
  d <- bind_rows(pres %>% select(longitude, latitude) %>% mutate(pr_ab = 1),
                 bg %>% select(longitude, latitude) %>% mutate(pr_ab = 0))
  env <- terra::extract(predictors,
                        vect(d, geom = c("longitude", "latitude"), crs = "EPSG:4326"),
                        ID = FALSE)
  d <- as.data.frame(bind_cols(d, env) %>% drop_na())
  
  # lq only: at 119-343 presences, hinge features are unstable (glmnet fails to
  # converge for C. horizontalis) and prone to overfitting. Settings held
  # constant across species so performance is comparable.
  
  m <- maxnet(p = d$pr_ab,
              data = d %>% select(all_of(names(predictors))),
              f = maxnet.formula(d$pr_ab, d %>% select(all_of(names(predictors))),
                                 classes = FC),
              regmult = REG)
  
  list(model = m, data = d, n_pres = sum(d$pr_ab == 1))
}


native_fits <- map(focal_species, fit_native_only) %>%
set_names(focal_species)

map_int(native_fits, "n_pres")

native_proj <- map(native_fits, 
                   ~ mask(terra::predict(pred_coarse,.x$model,
                                         fun = maxnet_fun, na.rm = TRUE), land_v))


for (sp in focal_species) {
  writeRaster(native_proj[[sp]],
              file.path(fig_dir, paste0("raster_nativeonly_",
                                        str_replace_all(sp, " ", "_"), ".tif")),
              overwrite = TRUE)
}


# ---- Thresholds --------------------------------------------------------------

# Derived from NATIVE training presences only, so naturalised records remain
# unseen and the transfer test stays independent.

thresholds <- map_dfr(focal_species, function(sp) {
  d <- native_fits[[sp]]$data
  
  ps <- terra::extract(native_proj[[sp]],
                       vect(d %>% filter(pr_ab == 1),
                            geom = c("longitude", "latitude"), 
                            crs = "EPSG:4326"),
                       ID = FALSE)[[1]]
  
  bgv <- terra::extract(native_proj[[sp]],
                        vect(d %>% filter(pr_ab == 0),
                             geom = c("longitude", "latitude"), 
                             crs = "EPSG:4326"),
                        ID = FALSE)[[1]]
  
  ps <- ps[!is.na(ps)]; bgv <- bgv[!is.na(bgv)]
  
  cand <- seq(0, 1, by = 0.005)
  sss <- map_dbl(cand, ~ mean(ps >= .x) + mean(bgv < .x))
  
  tibble(species = sp,
         mtp = min(ps),
         p10 = quantile(ps, 0.10),
         maxSSS = cand[which.max(sss)])
})

thresholds

write_csv(thresholds, file.path(tab_dir, "tableS5_thresholds.csv"))


################################################################################

# 6. MODEL EVALUATION

################################################################################


# ---- 6a. Spatially blocked cross-validation ----------------------------------

# Random CV inflates performance on spatially autocorrelated occurrence data.
# Block size from the autocorrelation range of the predictors themselves; using
# autocorrelation in pr_ab would partly reflect background sampling instead.


ac <- cv_spatial_autocor(r = pred_coarse, num_sample = 5000,
                         plot = FALSE, progress = FALSE)

block_size_m <- ac$range
block_size_m


run_block_cv <- function(sp, k = 5) {
  
  d <- native_fits[[sp]]$data
  d_sf <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
  
  folds <- cv_spatial(x = d_sf, column = "pr_ab", k = k, size = block_size_m,
                      selection = "random", iteration = 100,
                      progress = FALSE, plot = FALSE, report = FALSE)
  
  fid <- folds$folds_ids
  
  map_dfr(seq_len(k), function(i) {
    train <- d[fid != i, ]
    test <- d[fid == i, ]
    if (sum(train$pr_ab == 1) < 20 || sum(test$pr_ab == 1) < 5) {
      return(tibble(species = sp, fold = i, auc = NA_real_, boyce = NA_real_,
                    n_test_pres = sum(test$pr_ab == 1)))
    }
    
    trd <- as.data.frame(train %>% select(all_of(names(predictors))))
    m <- maxnet(p = train$pr_ab, data = trd,
                f = maxnet.formula(train$pr_ab, trd, classes = FC), regmult = REG)
    p_test <- pred_pts(m, test %>% select(all_of(names(predictors))))
    
    a <- precrec::auc(precrec::evalmod(scores = p_test, labels = test$pr_ab))
    b <- ecospat.boyce(fit = p_test, obs = p_test[test$pr_ab == 1], PEplot = FALSE)
    
    tibble(species = sp, fold = i,
           auc = a$aucs[a$curvetypes == "ROC"], boyce = b$cor,
           n_test_pres = sum(test$pr_ab == 1))
    
  })
}

block_cv <- map_dfr(focal_species, run_block_cv)

block_cv_summary <-
  block_cv %>%
  group_by(species) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            mean_boyce = mean(boyce, na.rm = TRUE),
            sd_boyce = sd(boyce, na.rm = TRUE),
            folds_used = sum(!is.na(auc)), .groups = "drop")

block_cv_summary

write_csv(block_cv, file.path(tab_dir, "tableS6_blockcv_folds.csv"))

# ---- 6b. Geographic transfer -------------------------------------------------

# Boyce computed over the projection RESTRICTED to the naturalised accessible
# area. The global version partly rewards correctly excluding regions no
# Cotoneaster could occupy, and inflates the statistic.

transfer_eval <- map_dfr(focal_species, function(sp) {
  ca <- calib_areas[[paste(sp, "naturalised", sep = "_")]]
  r <- mask(crop(native_proj[[sp]], ca), ca)
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  suit <- terra::extract(r, vect(inv, geom = c("longitude", "latitude"),
                                 crs = "EPSG:4326"), ID = FALSE)[[1]]
  suit <- suit[!is.na(suit)]
  
  tibble(species = sp, n_eval = length(suit),
         boyce_restricted = ecospat.boyce(values(r, na.rm = TRUE), suit,
                                          PEplot = FALSE)$cor)
})

transfer_eval

# ---- 6c. Binary performance --------------------------------------------------

binary_eval <- map_dfr(focal_species, function(sp) {
  th <- thresholds %>% filter(species == sp) %>% pull(p10)
  b <- native_proj[[sp]] >= th
  a <- cellSize(b, unit = "km")
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  hit <- terra::extract(b, vect(inv, geom = c("longitude", "latitude"),
                                crs = "EPSG:4326"), ID = FALSE)[[1]]
  
  tibble(species = sp, threshold = th,
         pct_land = 100 * as.numeric(global(a * b, "sum", na.rm = TRUE)) /
           as.numeric(global(a * !is.na(b), "sum", na.rm = TRUE)),
         pct_naturalised_captured = 100 * mean(hit, na.rm = TRUE))
})


binary_eval %>% mutate(enrichment = pct_naturalised_captured / pct_land)


# ---- 6d. Threshold sensitivity ----------------------------------------------

threshold_sensitivity <-
  map_dfr(c("mtp", "p10", "maxSSS"), function(rule) {
    
    map_dfr(focal_species, function(sp) {
      th <- thresholds %>% filter(species == sp) %>% pull(!!rule)
      b <- native_proj[[sp]] >= th
      a <- cellSize(b, unit = "km")
      inv <- occ_thin %>% filter(species == sp, range == "naturalised")
      hit <- terra::extract(b, vect(inv, geom = c("longitude", "latitude"),
                                    crs = "EPSG:4326"), ID = FALSE)[[1]]
      
      tibble(rule = rule, species = sp, threshold = th,
             pct_land = 100 * as.numeric(global(a * b, "sum", na.rm = TRUE)) /
               as.numeric(global(a * !is.na(b), "sum", na.rm = TRUE)),
             pct_captured = 100 * mean(hit, na.rm = TRUE))
    })
  })

threshold_sensitivity

write_csv(threshold_sensitivity, file.path(tab_dir, "tableS7_threshold_sensitivity.csv"))


# ---- 6e. Clamping sensitivity ------------------------------------------------

native_proj_noclamp <- map(native_fits, ~ mask(

terra::predict(pred_coarse, .x$model, fun = maxnet_fun_noclamp, na.rm = TRUE), land_v))


clamp_comparison <- 
  map_dfr(focal_species, function(sp) {
    
  th <- thresholds %>% filter(species == sp) %>% pull(p10)
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  iv <- vect(inv, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  b_cl <- native_proj[[sp]] >= th
  b_nc <- native_proj_noclamp[[sp]] >= th
  a <- cellSize(b_cl, unit = "km")
  land_km <- as.numeric(global(a * !is.na(b_cl), "sum", na.rm = TRUE))
  
  tibble(species = sp,
         pct_land_clamped = 100 * as.numeric(global(a * b_cl, "sum", na.rm = TRUE)) / land_km,
         pct_land_unclamped = 100 * as.numeric(global(a * b_nc, "sum", na.rm = TRUE)) / land_km,
         captured_clamped = 100 * mean(terra::extract(b_cl, iv, ID = FALSE)[[1]], na.rm = TRUE),
         captured_unclamped = 100 * mean(terra::extract(b_nc, iv, ID = FALSE)[[1]], na.rm = TRUE))
})

clamp_comparison

write_csv(clamp_comparison, file.path(tab_dir, "tableS8_clamping.csv"))

################################################################################

# 7. NULL MODEL COMPARISON

################################################################################

# Boyce has no fixed null expectation. Nulls draw presence sets of the same size
# from within the NATIVE accessible area, so any performance reflects the signal
# in the specific occurrence locations rather than the geometry of the region.

# CRITICAL: nulls must be evaluated exactly as observed models are - restricted
# to the naturalised accessible area. Evaluating nulls over the full global
# projection rewards correct rejection of remote regions and inflates them.

run_null <- function(sp, n_null = N_NULL) {
  
  d_obs <- native_fits[[sp]]$data
  n_pres <- sum(d_obs$pr_ab == 1)
  ca <- calib_areas[[paste(sp, "native", sep = "_")]]
  ca_inv <- calib_areas[[paste(sp, "naturalised", sep = "_")]]
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  inv_v <- vect(inv, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  bg <- d_obs %>% filter(pr_ab == 0)
  
  map_dfr(seq_len(n_null), function(i) {
    
    pts <- terra::spatSample(ca, size = ceiling(n_pres * 2), method = "random")
    xy <- terra::crds(pts) %>% as_tibble() %>% rename(longitude = x, latitude = y)
    env <- terra::extract(predictors,
                          vect(xy, geom = c("longitude", "latitude"),
                               crs = "EPSG:4326"), ID = FALSE)
    keep <- complete.cases(env)
    xy <- xy[keep, ]; env <- env[keep, ]
    
    if (nrow(xy) < n_pres) {
      return(tibble(species = sp, rep = i, boyce = NA_real_, n_pres = n_pres))
    }
    
    take <- sample(nrow(xy), n_pres)
    
    d <- bind_rows(
      bind_cols(xy[take, ], env[take, ]) %>% mutate(pr_ab = 1),
      bg %>% select(longitude, latitude, all_of(names(predictors)), pr_ab))
    
    dd <- as.data.frame(d %>% select(all_of(names(predictors))))
    m <- try(maxnet(p = d$pr_ab, data = dd,
                    f = maxnet.formula(d$pr_ab, dd, classes = FC),
                    regmult = REG), silent = TRUE)
    if (inherits(m, "try-error")) {
      return(tibble(species = sp, rep = i, boyce = NA_real_, n_pres = n_pres))
    }
    
    r <- mask(terra::predict(pred_coarse, m, fun = maxnet_fun, na.rm = TRUE), land_v)
    r_eval <- mask(crop(r, ca_inv), ca_inv)
    suit <- terra::extract(r_eval, inv_v, ID = FALSE)[[1]]
    suit <- suit[!is.na(suit)]
    b <- try(ecospat.boyce(fit = values(r_eval, na.rm = TRUE), obs = suit,
                           PEplot = FALSE)$cor, silent = TRUE)
    
    tibble(species = sp, rep = i,
           boyce = if (inherits(b, "try-error")) NA_real_ else b,
           n_pres = n_pres)
  })
}

null_results <- map_dfr(focal_species, ~ run_null(.x, n_null = N_NULL))

null_results %>% distinct(species, n_pres) # confirm thinning used

null_pct <-
  null_results %>%
  left_join(transfer_eval %>% select(species, obs = boyce_restricted), by = "species") %>%
  group_by(species) %>%
  summarise(null_boyce_median = median(boyce, na.rm = TRUE),
            pct_nulls_below_obs = 100 * mean(boyce < obs, na.rm = TRUE),
            n_null = sum(!is.na(boyce)), .groups = "drop")

null_pct

write_csv(null_results %>% select(species, rep, boyce),
          file.path(tab_dir, "tableS9_null_reps.csv"))


################################################################################

# 8. REGIONAL BREAKDOWN

################################################################################

# A global capture rate can conceal complete failure in a whole region.

capture_by_region <-
  map_dfr(focal_species, function(sp) {
    th <- thresholds %>% filter(species == sp) %>% pull(p10)
    b <- native_proj[[sp]] >= th
    
    inv <- occ_thin %>% 
      filter(species == sp, range == "naturalised") %>%
      left_join(gbif_final_sp %>% select(longitude, latitude, tdwg3) %>% distinct(),
                by = c("longitude", "latitude"))
    
    hit <- terra::extract(b,vect(inv, geom = c("longitude", "latitude"),
                                 crs = "EPSG:4326"), ID = FALSE)[[1]]
    inv %>% 
      mutate(hit = hit) %>%
      group_by(tdwg3) %>%
      summarise(n = n(), pct_captured = 100 * mean(hit, na.rm = TRUE), .groups
                = "drop") %>%
      mutate(species = sp)
  })

capture_by_region %>% 
  filter(n >= 20) %>% 
  arrange(species, pct_captured) %>% 
  print(n = 60)

write_csv(capture_by_region, file.path(tab_dir, "tableS10_capture_by_region.csv"))

# ---- Climate of failing units ------------------------------------------------

poor_units <- capture_by_region %>% filter(pct_captured < 50, n >= 20)

poor_climate <- map_dfr(seq_len(nrow(poor_units)), function(i) {
  sp <- poor_units$species[i]; u <- poor_units$tdwg3[i]
  
  pts <- occ_thin %>% 
    filter(species == sp, range == "naturalised") %>%
    left_join(gbif_final_sp %>% select(longitude, latitude, tdwg3) %>% distinct(),
              by = c("longitude", "latitude")) %>% 
    filter(tdwg3 == u)
  
  e <- terra::extract(predictors, vect(pts, geom = c("longitude","latitude"),
                                       crs = "EPSG:4326"), ID = FALSE)
  
  tibble(species = sp, tdwg3 = u, n = nrow(pts),
         bio19 = median(e$bio19, na.rm = TRUE),
         bio15 = median(e$bio15, na.rm = TRUE),
         bio6 = median(e$bio6, na.rm = TRUE))
})

poor_climate

write_csv(poor_climate, file.path(tab_dir, "tableS11_failing_unit_climate.csv"))


# Native reference medians for comparison

native_climate <- map_dfr(focal_species, function(sp) {
  as.data.frame(native_fits[[sp]]$data %>% filter(pr_ab == 1) %>%
                  select(all_of(names(predictors)))) %>%
    summarise(across(everything(), median)) %>%
    mutate(species = sp, .before = 1)
})

native_climate


################################################################################

# 9. EXTRAPOLATION DIAGNOSTICS

################################################################################

# MESS takes the MINIMUM across variables, so a single novel dimension drives
# the aggregate. The per-variable breakdown is essential context.

mess_maps <- map(focal_species, function(sp) {
  train <- as.data.frame(native_fits[[sp]]$data %>% filter(pr_ab == 1) %>%
                           select(all_of(names(predictors))))
  mask(mess(x = pred_coarse, v = train, full = FALSE), land_v)
  
}) %>% set_names(focal_species)


for (sp in focal_species) {
  writeRaster(mess_maps[[sp]],
              file.path(fig_dir, paste0("raster_mess_", str_replace_all(sp, " ", "_"), ".tif")),
              overwrite = TRUE)
}


mess_at_occ <- map_dfr(focal_species, function(sp) {
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  v <- terra::extract(mess_maps[[sp]],vect(inv, geom = c("longitude", "latitude"),
                                           crs = "EPSG:4326"), ID = FALSE)[[1]]
  v <- v[!is.na(v)]
  tibble(species = sp, n = length(v), pct_novel = 100 * mean(v < 0),
         median_mess = median(v), min_mess = min(v))
})


mess_by_var <- map_dfr(focal_species, function(sp) {
  train <- as.data.frame(native_fits[[sp]]$data %>% 
                           filter(pr_ab == 1) %>% 
                           select(all_of(names(predictors))))
  
  mf <- mask(mess(x = pred_coarse, v = train, full = TRUE), land_v)
  inv <- occ_thin %>% filter(species == sp, range == "naturalised")
  
  terra::extract(mf, vect(inv, geom = c("longitude", "latitude"),
                          crs = "EPSG:4326"), ID = FALSE) %>%
    summarise(across(everything(), ~ 100 * mean(.x < 0, na.rm = TRUE))) %>%
    mutate(species = sp)
})

mess_at_occ; mess_by_var

write_csv(mess_at_occ, file.path(tab_dir, "tableS12_mess_aggregate.csv"))
write_csv(mess_by_var, file.path(tab_dir, "tableS13_mess_by_variable.csv"))

################################################################################

# 10. RESPONSE CURVES

################################################################################

# All other variables held at their means. With lq features there are no

# interaction terms, so this is a fair representation of the fitted model.

response_data <-
  map_dfr(focal_species, function(sp) {
    d <- as.data.frame(native_fits[[sp]]$data %>%
                         select(all_of(names(predictors))))
    m <- native_fits[[sp]]$model
    means <- map_dbl(d, mean)
    
    map_dfr(names(predictors), function(v) {
      rng <- range(d[[v]], na.rm = TRUE)
      grid <- tibble(!!v := seq(rng[1], rng[2], length.out = 100))
      for (other in setdiff(names(predictors), v)) 
        
        grid[[other]] <- means[[other]]
      
      grid <- as.data.frame(grid[, names(predictors)])
      
      tibble(species = sp, variable = v,
             value = seq(rng[1], rng[2], length.out = 100),
             suitability = pred_pts(m, grid, clamp = FALSE))
    })
  })

write_csv(response_data, file.path(tab_dir, "tableS14_response_curves.csv"))


# Native bio19 distribution - shows how far the training data are censored

map_dfr(focal_species, function(sp) {
  tr <- as.data.frame(native_fits[[sp]]$data %>% filter(pr_ab == 1))
  
  tibble(species = sp,
         bio19_q50 = quantile(tr$bio19, 0.50),
         bio19_q90 = quantile(tr$bio19, 0.90),
         bio19_q95 = quantile(tr$bio19, 0.95),
         bio19_max = max(tr$bio19))
})

################################################################################

# 11. RECIPROCAL TEST (naturalised -> native)

################################################################################

# If the native range is merely under-documented, the far better sampled
# naturalised range should contain the native climate within its envelope, and
# transfer should succeed in this direction. Symmetric failure would instead
# indicate genuine divergence.

fit_naturalised_only <- function(sp) {
  pres <- occ_thin %>% filter(species == sp, range == "naturalised")
  bg <- backgrounds[[paste(sp, "naturalised", sep = "_")]]
  d <- bind_rows(pres %>% select(longitude, latitude) %>% mutate(pr_ab = 1),
                 bg %>% select(longitude, latitude) %>% mutate(pr_ab = 0))
  
  env <- terra::extract(predictors,
                        vect(d, geom = c("longitude", "latitude"), crs = "EPSG:4326"),
                        ID = FALSE)
  
  d <- as.data.frame(bind_cols(d, env) %>% drop_na())
  dd <- d %>% select(all_of(names(predictors)))
  
  m <- maxnet(p = d$pr_ab, data = dd,
              f = maxnet.formula(d$pr_ab, dd, classes = FC), regmult = REG)
  
  list(model = m, data = d, n_pres = sum(d$pr_ab == 1))
}


natzd_fits <- map(focal_species, fit_naturalised_only) %>%
  set_names(focal_species)

natzd_proj <- map(natzd_fits, ~ mask(terra::predict(pred_coarse,.x$model,
                                                    fun = maxnet_fun, na.rm = TRUE),
                                     land_v))

reciprocal_eval <- map_dfr(focal_species, function(sp) {
  d <- natzd_fits[[sp]]$data
  ca_nat <- calib_areas[[paste(sp, "native", sep = "_")]]
  r_eval <- mask(crop(natzd_proj[[sp]], ca_nat), ca_nat)
  ps <- terra::extract(natzd_proj[[sp]],
                       vect(d %>% filter(pr_ab == 1),
                            geom = c("longitude", "latitude"), crs = "EPSG:4326"),
                       ID = FALSE)[[1]]
  
  th <- quantile(ps, 0.10, na.rm = TRUE)
  nat <- occ_thin %>% filter(species == sp, range == "native")
  suit <- terra::extract(r_eval, vect(nat, geom = c("longitude", "latitude"),
                                      crs = "EPSG:4326"), ID = FALSE)[[1]]
  suit <- suit[!is.na(suit)]
  a <- cellSize(natzd_proj[[sp]], unit = "km")
  
  tibble(species = sp, direction = "naturalised -> native",
         n_train = natzd_fits[[sp]]$n_pres, n_eval = length(suit),
         threshold_p10 = th,
         boyce_restricted = ecospat.boyce(values(r_eval, na.rm = TRUE), suit, PEplot = FALSE)$cor,
         pct_land = 100 * as.numeric(global(a * (natzd_proj[[sp]] >= th), "sum", na.rm = TRUE)) /
           as.numeric(global(a * !is.na(natzd_proj[[sp]]), "sum", na.rm = TRUE)),
         pct_captured = 100 * mean(suit >= th))
})

reciprocal_eval

# ---- Sample-matched reciprocal -----------------------------------------------

# Naturalised models train on 20-40x more presences, so raw superiority could

# be sample size alone. Matching to native n controls for this.

fit_naturalised_matched <- function(sp, n_reps = N_MATCH) {
  d_full <- natzd_fits[[sp]]$data
  n_target <- occ_thin %>% 
    filter(species == sp, range == "native") %>%
    nrow()
  
  pres_idx <- which(d_full$pr_ab == 1)
  bg_idx <- which(d_full$pr_ab == 0)
  nat <- occ_thin %>% filter(species == sp, range == "native")
  nat_v <- vect(nat, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  ca_nat <- calib_areas[[paste(sp, "native", sep = "_")]]
  
  
  map_dfr(seq_len(n_reps), function(i) {
    d <- d_full[c(sample(pres_idx, n_target), bg_idx), ]
    dd <- d %>% select(all_of(names(predictors)))
    m <- try(maxnet(p = d$pr_ab, data = dd,
                    f = maxnet.formula(d$pr_ab, dd, classes = FC),
                    regmult = REG), silent = TRUE)
    if (inherits(m, "try-error")) {
      return(tibble(species = sp, rep = i, boyce = NA_real_, pct_captured =
                      NA_real_))
    }

    r <- mask(terra::predict(pred_coarse, m, fun = maxnet_fun, na.rm = TRUE), land_v)
    r_eval <- mask(crop(r, ca_nat), ca_nat)
    ps <- terra::extract(r, vect(d %>% filter(pr_ab == 1),
                                 geom = c("longitude", "latitude"),
                                 crs = "EPSG:4326"), ID = FALSE)[[1]]
    
    th <- quantile(ps, 0.10, na.rm = TRUE)
    suit <- terra::extract(r_eval, nat_v, ID = FALSE)[[1]]
    suit <- suit[!is.na(suit)]
    b <- try(ecospat.boyce(values(r_eval, na.rm = TRUE), suit, PEplot =
                             FALSE)$cor, silent = TRUE)
    
    tibble(species = sp, rep = i,
           boyce = if (inherits(b, "try-error")) NA_real_ else b,
           pct_captured = 100 * mean(suit >= th, na.rm = TRUE))
  })
}


matched_reps <- map_dfr(focal_species, fit_naturalised_matched)


matched_summary <-
  matched_reps %>%
  group_by(species) %>%
  summarise(n_reps = sum(!is.na(boyce)),
            boyce_median = median(boyce, na.rm = TRUE),
            boyce_q05 = quantile(boyce, 0.05, na.rm = TRUE),
            boyce_q95 = quantile(boyce, 0.95, na.rm = TRUE),
            capture_median = median(pct_captured, na.rm = TRUE), .groups = "drop")

matched_summary


direction_comparison <-
  bind_rows(
    transfer_eval %>% 
      select(species, boyce = boyce_restricted) %>%
      left_join(binary_eval %>% select(species,
                                       pct_captured = pct_naturalised_captured,
                                       pct_land), by = "species") %>%
      mutate(direction = "native -> naturalised", model = "full sample"),
    
    reciprocal_eval %>% 
      select(species, boyce = boyce_restricted, pct_captured, pct_land) %>%
      mutate(direction = "naturalised -> native", model = "full sample"),
    
    matched_summary %>%
      select(species, boyce = boyce_median, pct_captured = capture_median) %>%
      mutate(direction = "naturalised -> native", model = "sample-matched", pct_land = NA_real_)
  ) %>%
  
  mutate(enrichment = pct_captured / pct_land) %>%
  select(species, direction, model, boyce, pct_captured, pct_land, enrichment) %>%
  arrange(species, direction, model)



direction_comparison

write_csv(direction_comparison, file.path(tab_dir, "table3_direction_comparison.csv"))

write_csv(matched_reps, file.path(tab_dir, "tableS15_reciprocal_matched_reps.csv"))



# ---- bio19 response under naturalised calibration ----------------------------

# The forward models learned a monotonic preference for dry winters because the
# native range is censored on bio19. The naturalised range is not censored, so
# its fitted response should admit the dry native climate.


resp_reciprocal <- map_dfr(focal_species, function(sp) {
  d <- natzd_fits[[sp]]$data %>% select(all_of(names(predictors)))
  means <- map_dbl(d, mean)
  rng <- range(d$bio19, na.rm = TRUE)
  grid <- tibble(bio19 = seq(rng[1], rng[2], length.out = 100))
  
  for (v in setdiff(names(predictors), "bio19")) grid[[v]] <- means[[v]]
  
  tibble(species = sp, direction = "naturalised-fitted",
         value = seq(rng[1], rng[2], length.out = 100),
         suitability = pred_pts(natzd_fits[[sp]]$model,
                                as.data.frame(grid[, names(predictors)]),
                                clamp = FALSE))
})


################################################################################

# 12. FIGURES

################################################################################

coast <- ne_coastline(scale = 50, returnclass = "sf")


# ---- Fig A: COUE niche panels ------------------------------------------------

# dyn raster codes (verified empirically against the density grids):
# 0 = neither range
# 1, 2 = unfilling (native only; analogous / non-analogous)
# 3 = stability (both)
# 4, 5 = expansion (naturalised only; non-analogous / analogous)



coue_to_df <- function(res, sp, intersection = 0.1) {
  d <- ecospat.niche.dyn.index(res$grid_native, res$grid_naturalised,
                               intersection = intersection)
  dyn_df <- as.data.frame(d$dyn, xy = TRUE, na.rm = FALSE)
  names(dyn_df)[3] <- "cat"
  
  dyn_df %>%
    mutate(category = case_when(cat %in% c(1, 2) ~ "Unfilling",
                                cat == 3 ~ "Stability",
                                cat %in% c(4, 5) ~ "Expansion",
                                TRUE ~ NA_character_),
           category = factor(category, levels = c("Expansion", "Stability", "Unfilling")),
           species = sp) %>%
    filter(!is.na(category))
}

coue_outline_df <- function(res, sp) {
  bind_rows(
    as.data.frame(res$grid_native$Z, xy = TRUE, na.rm = FALSE) %>%
      rename(value = 3) %>% mutate(range = "Native"),
    as.data.frame(res$grid_naturalised$Z, xy = TRUE, na.rm = FALSE) %>%
      rename(value = 3) %>% mutate(range = "Naturalised")
  ) %>% mutate(species = sp)
}


build_coue_panel <- function(res, sp, panel_letter, av) {
  cats <- coue_to_df(res, sp)
  outl <- coue_outline_df(res, sp)
  dyn <- ecospat.niche.dyn.index(res$grid_native, res$grid_naturalised,
                                 intersection = 0.1)$dynamic.index.w
  
  ggplot() +
    geom_raster(data = cats, aes(x, y, fill = category)) +
    geom_contour(data = outl, aes(x, y, z = value, colour = range),
                 breaks = 1e-6, linewidth = 0.4) +
    scale_fill_manual(values = c(Expansion = "#E8A0A0", Stability = "#9C9CE8",
                                 Unfilling = "#A0D8A0"), drop = FALSE, name = NULL) +
    scale_colour_manual(values = c(Native = "#1E8449", Naturalised =
                                     "#B03A2E"),
                        name = "Background extent") +
    coord_equal() +
    labs(title = bquote(bold(.(panel_letter)) ~ italic(.(sp))),
         subtitle = sprintf("E = %.2f S = %.2f U = %.2f D = %.3f",
                            dyn["expansion"], dyn["stability"],
                            dyn["unfilling"], res$overlap$D),
         x = sprintf("PC1 (%.1f%%)", av$pc1),
         y = sprintf("PC2 (%.1f%%)", av$pc2)) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 7, colour = "grey25"),
          legend.key.size = unit(0.35, "cm"))
}


panels <- map(seq_along(focal_species), function(i) {
  sp <- focal_species[i]
  build_coue_panel(coue[[sp]], sp, paste0("(", letters[i], ")"),
                   axis_var %>% filter(species == sp))
})

figA <- wrap_plots(panels, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figA_coue.png"), figA, width = 11, height = 4.2, dpi = 300)

# ---- Fig B: binary suitability with naturalised records ----------------------

map_binary <- function(sp, th) {
  b <- native_proj[[sp]] >= th
  df <- as.data.frame(b, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "suitable"
  df$suitable <- as.integer(df$suitable) # logical -> 0/1 for the manual scale
  pts <- occ_thin %>% filter(species == sp, range == "naturalised")
  
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = factor(suitable))) +
    geom_sf(data = coast, colour = "grey30", linewidth = 0.12) +
    geom_point(data = pts, aes(longitude, latitude), shape = 21,
               colour = "black", fill = NA, size = 0.45, stroke = 0.18, alpha = 0.6) +
    scale_fill_manual(values = c("0" = "grey92", "1" = "#B03A2E"),
                      labels = c("Unsuitable", "Suitable"), name = NULL) +
    coord_sf(expand = FALSE, ylim = c(-58, 84)) +
    labs(title = bquote(italic(.(sp))), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 6))
}


for (sp in focal_species) {
  th <- thresholds %>% filter(species == sp) %>% pull(p10)
  ggsave(file.path(fig_dir, paste0("figB_binary_", str_replace_all(sp, " ", "_"), ".png")),
         map_binary(sp, th), width = 9, height = 4.5, dpi = 300)
}



# ---- Fig C: bio19 response, native-fitted vs naturalised-fitted --------------

# The key mechanistic figure: naturalised-fitted models admit the dry native

# climate, while native-fitted models exclude the wet naturalised climate.


p_recip <- map(focal_species, function(sp) {
  tr_nat <- occ_thin %>% filter(species == sp, range == "native")
  env_nat <- terra::extract(predictors,
                            vect(tr_nat, geom = c("longitude", "latitude"),
                                 crs = "EPSG:4326"), ID = FALSE)
  bind_rows(
    response_data %>% 
      filter(variable == "bio19", species == sp) %>%
      mutate(direction = "native-fitted") %>% 
      select(value, suitability, direction),
    
    resp_reciprocal %>% 
      filter(species == sp) %>% 
      select(value, suitability, direction)
  ) %>%
    
    ggplot(aes(value, suitability, colour = direction)) +
    geom_line(linewidth = 0.6) +
    geom_rug(data = tibble(value = env_nat$bio19), aes(x = value),
             inherit.aes = FALSE, alpha = 0.3, colour = "#1E8449") +
    scale_colour_manual(values = c("native-fitted" = "#1E8449","naturalised-fitted" = "#B03A2E"), name = NULL) +
    labs(title = bquote(italic(.(sp))), x = "bio19 (mm)", y = "Predicted suitability") +
    theme_classic(base_size = 9)
  
})

ggsave(file.path(fig_dir, "figC_bio19_reciprocal.png"), 
       wrap_plots(p_recip, ncol = 3, guides = "collect") &
         theme(legend.position = "bottom"), width = 10, height = 3.4, dpi = 300)


# ---- Supplementary figures ---------------------------------------------------


# All response curves

p_resp <- response_data %>%
  mutate(sp_short = str_replace(species, "Cotoneaster ", "C. ")) %>%
  ggplot(aes(value, suitability, colour = sp_short)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "Predictor value", y = "Predicted suitability") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(face = "italic"))

ggsave(file.path(fig_dir, "figS_response_curves.png"), p_resp, width = 9, height = 5, dpi = 300)


# Null distributions

p_null <- null_results %>%
  left_join(transfer_eval %>% select(species, obs = boyce_restricted), by = "species") %>%
  mutate(sp_short = str_replace(species, "Cotoneaster ", "C. ")) %>%
  
  ggplot(aes(boyce)) +
  geom_histogram(bins = 20, fill = "grey75", colour = "white") +
  geom_vline(aes(xintercept = obs), colour = "#B03A2E", linewidth = 0.8) +
  facet_wrap(~ sp_short, scales = "free_y") +
  labs(x = "Boyce index", y = "Null replicates",
       caption = "Red line: observed value from the fitted native-range model") +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "italic"))

ggsave(file.path(fig_dir, "figS_null_distributions.png"), p_null, 
       width = 8, height = 3, dpi = 300)



# Capture by region

p_capture <- capture_by_region %>%
  filter(n >= 20) %>%
  mutate(sp_short = str_replace(species, "Cotoneaster ", "C. "),
         tdwg3 = tidytext::reorder_within(tdwg3, pct_captured, species)) %>%
  
  ggplot(aes(pct_captured, tdwg3, size = n)) +
  geom_point(colour = "#2C3E50", alpha = 0.8) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey60") +
  facet_wrap(~ sp_short, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_size_continuous(range = c(1, 5), name = "Records") +
  labs(x = "Naturalised records captured (%)", y = "WGSRPD level-3 unit") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "italic"))

ggsave(file.path(fig_dir, "figS_capture_by_region.png"), p_capture,
       width = 10, height = 5, dpi = 300)



# MESS maps

for (sp in focal_species) {
  df <- as.data.frame(mess_maps[[sp]], xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "mess"
  df$mess_c <- pmax(pmin(df$mess, 50), -50)
  pts <- occ_thin %>% filter(species == sp, range == "naturalised")
  
  p <- ggplot() +
    geom_raster(data = df, aes(x, y, fill = mess_c)) +
    geom_sf(data = coast, colour = "grey20", linewidth = 0.12) +
    geom_point(data = pts, aes(longitude, latitude), shape = 21, colour = "black",
               fill = NA, size = 0.4, stroke = 0.15, alpha = 0.5) +
    scale_fill_gradient2(low = "#B03A2E", mid = "grey95", high = "#1B4F72",
                         midpoint = 0, name = "MESS") +
    coord_sf(expand = FALSE, ylim = c(-58, 84)) +
    labs(title = bquote(italic(.(sp))), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 5))
  
  ggsave(file.path(fig_dir, paste0("figS_mess_", str_replace_all(sp, " ", "_"), ".png")),
         p, width = 9, height = 4.5, dpi = 300)
}


################################################################################

# 13. MAIN TABLE

################################################################################


table_A <-
  tibble(species = focal_species) %>%
  left_join(
    map_dfr(focal_species, function(sp) {
      r <- coue[[sp]]
      dyn <- ecospat.niche.dyn.index(r$grid_native, r$grid_naturalised,
                                     intersection = 0.1)$dynamic.index.w
      tibble(species = sp,
             n_native = r$n_native, n_naturalised = r$n_naturalised,
             pc12_var = 100 * sum(r$pca$eig[1:2]) / sum(r$pca$eig),
             schoener_D = r$overlap$D,
             expansion = dyn["expansion"], stability = dyn["stability"],
             unfilling = dyn["unfilling"],
             equivalency_p = r$equivalency$p.D, similarity_p = r$similarity$p.D)
    }), by = "species") %>%
  
  left_join(block_cv_summary %>%
              select(species, cv_auc = mean_auc, cv_auc_sd = sd_auc,
                     cv_boyce = mean_boyce, cv_boyce_sd = sd_boyce), by = "species") %>%
  left_join(transfer_eval %>% select(species, n_eval, boyce_restricted), by = "species") %>%
  left_join(null_pct, by = "species") %>%
  left_join(binary_eval %>%
              select(species, threshold_p10 = threshold,
                     pct_land_suitable = pct_land,
                     pct_captured = pct_naturalised_captured), by = "species") %>%
  
  mutate(enrichment = pct_captured / pct_land_suitable) %>%
  left_join(mess_at_occ %>% select(species, pct_novel), by = "species")%>%
  left_join(mess_by_var %>% select(species, bio6_novel = bio6_mess,
                                   bio15_novel = bio15_mess), by = "species") %>%
  left_join(reciprocal_eval %>% select(species, recip_boyce = boyce_restricted),
            by = "species") %>%
  left_join(matched_summary %>% select(species, recip_boyce_matched = boyce_median),
            by = "species")

table_A_fmt <-
  table_A %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)),
         species = str_replace(species, "Cotoneaster ", "C. ")) %>%
  pivot_longer(-species, names_to = "metric") %>%
  pivot_wider(names_from = species, values_from = value)

table_A_fmt %>% print(n = 40)

write_csv(table_A, file.path(tab_dir, "tableA_full.csv"))
write_csv(table_A_fmt, file.path(tab_dir, "tableA_formatted.csv"))

################################################################################

# 14. SAVE

################################################################################

saveRDS(list(
  occ_all = occ_all, occ_thin = occ_thin, thin_sensitivity = thin_sensitivity,
  keep_vars = keep_vars, calib_areas_expanse = map_dbl(calib_areas, ~ expanse(.x, unit = "km")),
  backgrounds = backgrounds,
  coue = coue, coue_sweep = coue_sweep, axis_var = axis_var,
  axis_meaning = axis_meaning, pca_loadings = pca_loadings,
  native_fits = native_fits, thresholds = thresholds,
  block_cv = block_cv, block_cv_summary = block_cv_summary,
  transfer_eval = transfer_eval, binary_eval = binary_eval,
  threshold_sensitivity = threshold_sensitivity, clamp_comparison = clamp_comparison,
  null_results = null_results, null_pct = null_pct,
  capture_by_region = capture_by_region, poor_climate = poor_climate,
  native_climate = native_climate,
  mess_at_occ = mess_at_occ, mess_by_var = mess_by_var,
  response_data = response_data, resp_reciprocal = resp_reciprocal,
  natzd_fits = natzd_fits, reciprocal_eval = reciprocal_eval,
  matched_reps = matched_reps, matched_summary = matched_summary,
  direction_comparison = direction_comparison,
  table_A = table_A
), file.path(out_main, "suitability_analysis_objects.rds"))


writeLines(capture.output(sessionInfo()),
           file.path(out_main, "sessionInfo_suitability.txt"))


################################################################################

# END

################################################################################
