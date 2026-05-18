#===============================================================
#        This code is inspired by Ferris et al. 2022 and       #
# written by Savi Koissi (2026), the entire pipeline is availa-#
# on our github repo                                           #
#==============================================================#



# This is to generate country and year specific host landscape file. 
# In this file in local is a .Tiff file but with the binary env we generated an rds
library(terra)
library(sf)
library(readxl)
library(dplyr)

base_dir <- "/home/savi/project/SurveyCBSD/SurveyCBSD/R_reproducible"
cassava_tif_path <- file.path(base_dir, "data/cassavaMap/CassavaMap_Prod_v1.tif")

print("Attempting native Base-R binary recovery of the base cassava map...")

# 1. Read raw binary files directly bypassing GDAL headers
file_size <- file.info(cassava_tif_path)$size
con <- file(cassava_tif_path, "rb")
raw_bytes <- readBin(con, what = "numeric", n = file_size, size = 4)
close(con)

# 2. Reconstruct the spatial matrix layout manually
total_elements <- length(raw_bytes)
cols <- 1201
rows <- floor(total_elements / cols)
clean_matrix <- matrix(raw_bytes[1:(rows * cols)], nrow = rows, ncol = cols)
clean_matrix[is.nan(clean_matrix)] <- 0

# 3. Rebuild a clean SpatRaster in-memory
cassava <- rast(clean_matrix)
ext(cassava) <- c(2.6685, 14.6788, 4.2730, 13.8944) # Nigeria bounding box
crs(cassava) <- "EPSG:4326"

print("Base cassava map successfully reconstructed via binary parsing!")

# --- 4. Read standard vector/excel data safely ---
whitefly <- read_excel(file.path(base_dir, "data/whitefly/dataWhiteflyNigeria_2022.xlsx"))
nigeria  <- vect(file.path(base_dir, "data/state_boundaries/nga_admin1.shp"))
if (crs(nigeria) == "") crs(nigeria) <- "EPSG:4326"

# --- 5. Spatial Processing ---
cassava_cropped <- crop(cassava, nigeria)
cassava_masked  <- mask(cassava_cropped, nigeria)

# --- 6. Aggregate Points ---
whitefly_agg <- whitefly %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    whitefly_mean = mean(Total_Whitefly_Count, na.rm = TRUE),
    disease_prev  = mean(disease, na.rm = TRUE),
    .groups = "drop"
  )

whitefly_sf   <- st_as_sf(whitefly_agg, coords = c("Longitude", "Latitude"), crs = 4326)
whitefly_vect <- vect(whitefly_sf)

# --- 7. Define Grid Template ---
template <- rast(
  ext(cassava_masked),
  resolution = 0.01,
  crs = "EPSG:4326"
)

# --- 8. Rasterize Vector Layers ---
whitefly_raster <- rasterize(whitefly_vect, template, field = "whitefly_mean", fun = mean, background = NA)
disease_raster  <- rasterize(whitefly_vect, template, field = "disease_prev", fun = mean, background = NA)

# --- 9. Resample Cassava to match model grid ---
cassava_resampled <- resample(cassava_masked, template, method = "bilinear")

# --- 10. Focal Smoothing ---
w_matrix <- matrix(1, 5, 5)
whitefly_smooth <- focal(whitefly_raster, w = w_matrix, fun = mean, na.rm = TRUE)
disease_smooth  <- focal(disease_raster,  w = w_matrix, fun = mean, na.rm = TRUE)

# --- 11. SAVE AS NATIVE R MATRICES (Bypassing GDAL Drivers Completely) ---
print("Converting spatial layers to matrices to bypass writeRaster...")
whitefly_mat <- as.matrix(whitefly_smooth, wide=TRUE)
disease_mat  <- as.matrix(disease_smooth, wide=TRUE)
cassava_mat  <- as.matrix(cassava_resampled, wide=TRUE)

saveRDS(whitefly_mat, file.path(base_dir, "data/whitefly/whitefly.rds"))
saveRDS(disease_mat,  file.path(base_dir, "data/whitefly/disease.rds"))
saveRDS(cassava_mat,  file.path(base_dir, "data/cassavaMap/cassava_nigeria.rds"))

print("SUCCESS: All spatial layers saved as clean .rds matrices!")
