
#attach_states.R
attach_states <- function(grid, raster_ref, states_sf) {
  
  # STEP 1: extract coordinates from SAME raster used for grid
  coords <- terra::xyFromCell(raster_ref, grid$cell_id)
  
  if (nrow(coords) != nrow(grid)) {
    stop("Coordinate mismatch between grid and raster_ref.")
  }
  
  # STEP 2: build sf object with correct CRS
  grid_sf <- sf::st_as_sf(
    data.frame(
      cell_id = grid$cell_id,
      x = coords[,1],
      y = coords[,2]
    ),
    coords = c("x", "y"),
    crs = sf::st_crs(raster_ref)
  )
  
  # STEP 3: align CRS of states
  states_sf <- sf::st_transform(states_sf, sf::st_crs(grid_sf))
  
  # STEP 4: spatial join (KEEP ALL GRID CELLS)
  joined <- sf::st_join(grid_sf, states_sf, left = TRUE)
  
  if (nrow(joined) == 0) {
    stop("Spatial join failed: no overlap between grid and states.")
  }
  
  if (!"adm1_name" %in% names(joined)) {
    stop("State join failed: adm1_name missing.")
  }
  
  # STEP 5: safe extraction
  state_vec <- as.character(joined$adm1_name)
  state_vec[is.na(state_vec)] <- "UNKNOWN"
  
  grid$state_id <- state_vec
  
  grid
}