# # landscape_model.R
# 
# library(parallel)
# 
# compute_infection_pressure_local <- function(grid, coords, alpha, beta, max_dist = 10000) {
#   
#   n <- nrow(grid)
#   
#   # detect cores
#   n_cores <- max(1, detectCores() - 1)
#   cl <- makeCluster(n_cores)
#   
#   # export needed variables to workers
#   clusterExport(cl, varlist = c(
#     "grid", "coords", "alpha", "beta", "max_dist", "power_kernel"
#   ), envir = environment())
#   
#   # split indices into chunks
#   chunks <- split(1:n, cut(1:n, n_cores, labels = FALSE))
#   
#   # parallel computation
#   results <- parLapply(cl, chunks, function(idx) {
#     
#     pressure_chunk <- numeric(length(idx))
#     
#     for (j in seq_along(idx)) {
#       
#       i <- idx[j]
#       
#       dx <- coords[,1] - coords[i,1]
#       dy <- coords[,2] - coords[i,2]
#       d  <- sqrt(dx^2 + dy^2)
#       
#       neighbors <- which(d < max_dist)
#       
#       if (length(neighbors) == 0) {
#         pressure_chunk[j] <- 0
#         next
#       }
#       
#       k <- power_kernel(d[neighbors], alpha, beta)
#       
#       val <- k *
#         grid$infected_prop[neighbors] *
#         grid$whitefly_density[neighbors]
#       
#       # safety against NA
#       val[is.na(val)] <- 0
#       
#       pressure_chunk[j] <- sum(val)
#     }
#     
#     pressure_chunk
#   })
#   
#   stopCluster(cl)
#   
#   # combine results
#   pressure <- unlist(results)
#   
#   # safety cleanup
#   pressure[is.na(pressure)] <- 0
#   
#   pressure
# }
# 
# 
# 
# 
# simulate_landscape <- function(grid, coords, years, params) {
#   
#   results <- vector("list", years)
#   
#   start_time <- Sys.time()
#   
#   for (t in seq_len(years)) {
#     
#     cat("Year:", t, "/", years, "\n")
#     
#     t0 <- Sys.time()
#     
#     pressure <- compute_infection_pressure_local(
#       grid, coords,
#       params$alpha,
#       params$beta,
#       max_dist = 5000
#     )
#     
#     # clean probabilities
#     pressure[is.na(pressure)] <- 0
#     pressure <- pmin(pmax(pressure, 0), 1)
#     
#     new_infections <- rbinom(
#       nrow(grid),
#       size = 1,
#       prob = pressure
#     )
#     
#     # smoother infection growth (avoid instant saturation)
#     grid$infected_prop <- pmin(
#       grid$infected_prop + new_infections * 0.05,
#       1
#     )
#     
#     results[[t]] <- grid
#     
#     cat("  Done in",
#         round(difftime(Sys.time(), t0, units = "secs"), 2),
#         "sec\n")
#   }
#   
#   cat("Total time:",
#       round(difftime(Sys.time(), start_time, units = "mins"), 2),
#       "minutes\n")
#   
#   results
# }
# 

# landscape_model.R

library(parallel)

#---------------------------------------------------
# Compute local infection pressure
#---------------------------------------------------
compute_infection_pressure_local <- function(
    grid,
    coords,
    alpha,
    beta,
    max_dist = 5000
) {
  
  n <- nrow(grid)
  
  pressure <- numeric(n)
  
  for (i in seq_len(n)) {
    
    # distances from focal cell
    dx <- coords[,1] - coords[i,1]
    dy <- coords[,2] - coords[i,2]
    
    d <- sqrt(dx^2 + dy^2)
    
    neighbors <- which(d > 0 & d < max_dist)
    
    if (length(neighbors) == 0) {
      pressure[i] <- 0
      next
    }
    
    # dispersal kernel
    k <- power_kernel(
      distance = d[neighbors],
      alpha = alpha,
      beta = beta
    )
    
    # infection contribution
    contribution <-
      k *
      grid$infected_prop[neighbors] *
      grid$whitefly_density[neighbors] *
      grid$cassava_density[neighbors]
    
    contribution[is.na(contribution)] <- 0
    
    pressure[i] <- sum(contribution)
  }
  
  pressure
}

#---------------------------------------------------
# Simulate landscape epidemic
#---------------------------------------------------
simulate_landscape <- function(
    grid,
    coords,
    years,
    params
) {
  
  # CRITICAL SAFETY CHECK
  if (nrow(grid) != nrow(coords)) {
    stop("grid and coords are not aligned")
  }
  
  results <- vector("list", years)
  
  for (t in seq_len(years)) {
    
    cat("Year:", t, "/", years, "\n")
    
    pressure <- compute_infection_pressure_local(
      grid = grid,
      coords = coords,
      alpha = params$alpha,
      beta = params$beta,
      max_dist = params$max_dist
    )
    
    # normalize pressure
    pressure <- pressure / max(pressure, na.rm = TRUE)
    
    pressure[is.na(pressure)] <- 0
    
    # probability of NEW infection
    new_infection_prob <-
      (1 - grid$infected_prop) * pressure
    
    new_infections <- rbinom(
      n = nrow(grid),
      size = 1,
      prob = pmin(new_infection_prob, 1)
    )
    
    # gradual epidemic growth
    grid$infected_prop <- pmin(
      grid$infected_prop +
        new_infections * params$infection_increment,
      1
    )
    
    results[[t]] <- grid
  }
  
  results
}
