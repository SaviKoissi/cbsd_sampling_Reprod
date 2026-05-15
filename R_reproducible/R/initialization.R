# initialization.R

# create_grid.R

create_grid <- function(cassava, whitefly) {
  
  if (!terra::compareGeom(cassava, whitefly)) {
    stop("rasters are not aligned")
  }
  
  cassava_vals <- terra::values(cassava)
  whitefly_vals <- terra::values(whitefly)
  
  valid <-
    !is.na(cassava_vals) &
    !is.na(whitefly_vals)
  
  grid <- data.frame(
    row_id = seq_len(sum(valid)),
    cell_id = which(valid),
    cassava_density = cassava_vals[valid],
    whitefly_density = whitefly_vals[valid],
    infected_prop = 0
  )
  
  grid
}
# 
# #initialization.R
# 
# initialize_infection <- function(grid) {
#   
#   probs <- grid$cassava_density
#   
#   probs[is.na(probs)] <- 0
#   probs[probs < 0] <- 0
#   
#   if (all(probs == 0)) {
#     stop("No valid cassava density for initialization.")
#   }
#   
#   probs <- probs / sum(probs)
#   
#   start_cell <- sample(seq_len(nrow(grid)), 1, prob = probs)
#   
#   grid$infected_prop[start_cell] <- runif(1, 0.05, 0.15)
#   
#   list(grid = grid, start_cell = start_cell)
# }

# initialization.R

initialize_infection <- function(grid) {
  
  probs <- grid$cassava_density
  
  probs[is.na(probs)] <- 0
  probs[probs < 0] <- 0
  
  probs <- probs / sum(probs)
  
  start_cell <- sample(
    seq_len(nrow(grid)),
    size = 1,
    prob = probs
  )
  
  # initialize only ONE location
  grid$infected_prop[start_cell] <- 0.1
  
  list(
    grid = grid,
    start_cell = start_cell
  )
}