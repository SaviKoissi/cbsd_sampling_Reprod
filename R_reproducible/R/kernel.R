#kernel.R

# Power-law dispersal kernel
power_kernel <- function(distance, alpha, beta) {
  beta * (distance + 1)^(-alpha)
}


# # Compute infection pressure across grid
# compute_infection_pressure <- function(grid, coords, alpha, beta) {
#   dist_mat <- as.matrix(dist(coords))
#   kernel_mat <- power_kernel(dist_mat, alpha, beta)
#   
#   as.vector(kernel_mat %*% grid$infected_prop)
# }
# 
# # infection pressure = focal convolution
# kernel_matrix <- focalMat(cassava_res, d = 5000, type = "Gauss")
# 
# pressure <- focal(
#   infected_raster,
#   w = kernel_matrix,
#   fun = sum,
#   na.rm = TRUE
# )


compute_infection_pressure <- function(grid, coords, alpha, beta, max_dist = 10000) {
  
  n <- nrow(grid)
  pressure <- numeric(n)
  
  for (i in seq_len(n)) {
    
    dx <- coords[,1] - coords[i,1]
    dy <- coords[,2] - coords[i,2]
    d  <- sqrt(dx^2 + dy^2)
    
    neighbors <- which(d < max_dist)
    
    if (length(neighbors) == 0) next
    
    k <- power_kernel(d[neighbors], alpha, beta)
    
    pressure[i] <- sum(
      k *
        grid$infected_prop[neighbors] *
        grid$whitefly_density[neighbors]
    )
  }
  
  pressure
}