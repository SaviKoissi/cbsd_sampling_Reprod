# grid_construction.R

#========================================================
# GRID CONSTRUCTION (FIXED: NA PROPAGATION ROOT FIX)
#========================================================

create_grid <- function(cassava, whitefly) {
  
  if (!terra::compareGeom(cassava, whitefly)) {
    stop("Cassava and whitefly rasters are not aligned.")
  }
  
  cassava_vals  <- terra::values(cassava)
  whitefly_vals <- terra::values(whitefly)
  
  valid <- !is.na(cassava_vals) & !is.na(whitefly_vals)
  
  grid <- data.frame(
    row_id = seq_len(sum(valid)),
    cell_id = which(valid),
    cassava_density = cassava_vals[valid],
    whitefly_density = whitefly_vals[valid],
    infected_prop = 0
  )
  
  # SAFE DEFAULT STATE (prevents NA cascade)
  grid$state_id <- NA_character_
  
  return(grid)
}
