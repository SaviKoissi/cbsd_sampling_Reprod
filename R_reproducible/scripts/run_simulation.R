# # run_surveillance.R
# library(terra)
# library(sf)
# library(readxl)
# library(parallel)
# 
# source("R/initialization.R")
# source("R/kernel.R")
# source("R/landscape_model.R")
# 
# # --- LOAD REAL DATA (REQUIRED) ---
# cassava <- rast("data/cassavaMap/cassava_nigeria.tif")
# whitefly <- rast("data/whitefly/whitefly.tif")
# 
# # Coordinates of grid cells (REQUIRED)
# coords <- terra::xyFromCell(cassava, 1:ncell(cassava))
# 
# # Create a common template
# template <- rast(
#   ext(cassava),
#   resolution = 0.01,   # ~1 km
#   crs = crs(cassava)
# )
# 
# # Aggregate cassava
# cassava_agg <- aggregate(cassava, fact = round(res(cassava)[1] / 0.01), fun = mean)
# 
# # 
# # cassava_res <- project(cassava_res, "EPSG:3857")
# # whitefly_res <- project(whitefly_res, "EPSG:3857")
# # 
# # coords <- terra::xyFromCell(cassava_res, 1:ncell(cassava_res))
# # coords <- coords[grid$cell_id, ]
# res(cassava)
# 
# # Resample both onto template
# cassava_res <- resample(cassava_agg, template, method = "bilinear")
# whitefly_res <- resample(whitefly, template, method = "bilinear")
# 
# # Check alignment (CRITICAL)
# ncell(cassava_res)
# ncell(whitefly_res)
# #compareGeom(cassava_res, whitefly_res)
# 
# # --- INITIALIZE ---
# grid <- create_grid(cassava_res, whitefly_res)
# init <- initialize_infection(grid)
# 
# # --- PARAMETERS (FROM PAPER / ABC — REQUIRED) ---
# params <- list(alpha = 2, beta = 0.5)  # ⚠️ placeholder
# 
# # --- RUN ---
# 
# sim <- simulate_landscape(
#   init$grid,
#   coords,
#   years = 30,
#   params = params
# )
# 
# saveRDS(sim, "outputs/simulation.rds")

# run_simulation.R

library(terra)

source("R/initialization.R")
source("R/kernel.R")
source("R/landscape_model.R")

#---------------------------------------------------
# Load rasters
#---------------------------------------------------
cassava <- rast("data/cassavaMap/cassava_nigeria.tif")

whitefly <- rast("data/whitefly/whitefly.tif")

#---------------------------------------------------
# Common template
#---------------------------------------------------
template <- rast(
  ext(cassava),
  resolution = 0.01,
  crs = crs(cassava)
)

#---------------------------------------------------
# Resample rasters
#---------------------------------------------------
cassava_res <- resample(
  cassava,
  template,
  method = "bilinear"
)

whitefly_res <- resample(
  whitefly,
  template,
  method = "bilinear"
)

#---------------------------------------------------
# Create grid
#---------------------------------------------------
grid <- create_grid(
  cassava_res,
  whitefly_res
)

#---------------------------------------------------
# CRITICAL FIX:
# coordinates ONLY for VALID CELLS
#---------------------------------------------------
coords <- terra::xyFromCell(
  cassava_res,
  grid$cell_id
)

#---------------------------------------------------
# Initialize epidemic
#---------------------------------------------------
init <- initialize_infection(grid)

#---------------------------------------------------
# Parameters
#---------------------------------------------------
params <- list(
  alpha = 2,
  beta = 0.01,
  max_dist = 5000,
  infection_increment = 0.05
)

#---------------------------------------------------
# Run simulation
#---------------------------------------------------
sim <- simulate_landscape(
  grid = init$grid,
  coords = coords,
  years = 30,
  params = params
)

#---------------------------------------------------
# Save output
#---------------------------------------------------
saveRDS(
  sim,
  "outputs/simulation.rds"
)
