#spatial_join_states.R

#========================================================
# SAFE STATE ASSIGNMENT (NO ROW ORDER DEPENDENCY)
#========================================================

attach_states <- function(grid, cassava_raster, states_sf) {
  
  coords <- terra::xyFromCell(cassava_raster, grid$cell_id)
  
  grid_sf <- sf::st_as_sf(
    data.frame(
      cell_id = grid$cell_id,
      x = coords[,1],
      y = coords[,2]
    ),
    coords = c("x", "y"),
    crs = sf::st_crs(states_sf)
  )
  
  grid_sf <- sf::st_join(grid_sf, states_sf["adm1_name"])
  
  # SAFE MERGE BACK USING cell_id (NOT ROW ORDER)
  state_map <- data.frame(
    cell_id = grid_sf$cell_id,
    state_id = as.character(grid_sf$adm1_name)
  )
  
  grid <- merge(grid, state_map, by = "cell_id", all.x = TRUE)
  
  # CRITICAL FIX: prevent NA state leakage
  grid$state_id[is.na(grid$state_id)] <- "UNKNOWN"
  
  return(grid)
}