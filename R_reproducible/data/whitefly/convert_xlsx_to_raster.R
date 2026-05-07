# Convert excel_file to raster_file 

library(terra)
library(sf)
library(readxl)
library(dplyr)
# read data

whitefly <- read_excel("data/whitefly/dataWhiteflyNigeria_2022.xlsx")
cassava <- rast("data/cassavaMap/CassavaMap_Prod_v1.tif")

nigeria <- vect("data/state_boundaries/nga_admin1.shp")

cassava <- crop(cassava, nigeria)
cassava <- mask(cassava, nigeria)

whitefly_agg <- whitefly %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    whitefly_mean = mean(Total_Whitefly_Count, na.rm = TRUE),
    disease_prev  = mean(disease, na.rm = TRUE),
    .groups = "drop"
  )


whitefly_sf <- st_as_sf(
  whitefly_agg,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

whitefly_vect <- vect(whitefly_sf)

template <- template <- rast(
  ext(cassava),
  resolution = 0.01,
  crs = crs(cassava)
)

# cassava_res <- resample(cassava, template, method = "bilinear")
# whitefly_res <- resample(whitefly, template, method = "bilinear")

whitefly_raster <- rasterize(
  whitefly_vect,
  template,
  field = "whitefly_mean",
  fun = mean,
  background = NA
)

disease_raster <- rasterize(
  whitefly_vect,
  template,
  field = "disease_prev",
  fun = mean,
  background = NA
)

whitefly_smooth <- focal(
  whitefly_raster,
  w = matrix(1, 5, 5),
  fun = mean,
  na.rm = TRUE
)

disease_smooth <- focal(
  disease_raster, 
  w = matrix(1, 5, 5), 
  fun = mean,
  na.rm = TRUE
)

writeRaster(whitefly_smooth, "data/whitefly/whitefly.tif", overwrite=TRUE)

writeRaster(disease_smooth, "data/whitefly/disease.tif", overwrite=TRUE)

writeRaster(cassava,"data/cassavaMap/cassava_nigeria.tif", overwrite=TRUE)

