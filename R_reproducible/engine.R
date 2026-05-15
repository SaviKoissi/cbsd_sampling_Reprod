#engine.R

library(terra)
library(sf)
library(dplyr)
library(parallel)

# --- DISPERSAL PHYSICS ---
power_kernel <- function(distance, alpha, beta) {
  beta * (distance + 1)^(-alpha)
}

compute_pressure <- function(grid, coords, alpha, beta, max_dist) {
  n <- nrow(grid)
  pressure <- numeric(n)
  # Vectorized distance calculation for speed
  for (i in seq_len(n)) {
    d <- sqrt((coords[,1] - coords[i,1])^2 + (coords[,2] - coords[i,2])^2)
    nb <- which(d > 0 & d < max_dist)
    if (length(nb) == 0) next
    
    k <- power_kernel(d[nb], alpha, beta)
    # Generic interaction: Kernel * Infected * Vector * Host
    val <- k * grid$infected_prop[nb] * grid$vector_density[nb] * grid$host_density[nb]
    pressure[i] <- sum(val, na.rm = TRUE)
  }
  return(pressure)
}

# --- SAMPLING STRATEGY CONTROLLER ---
# This is the 'Safe Sampler' we built, abstracted for any weights
run_sampler <- function(pool, n, weights = NULL) {
  if (length(pool) == 0) return(integer(0))
  n_pick <- min(n, length(pool))
  
  if (is.null(weights) || sum(weights, na.rm = TRUE) == 0) {
    return(sample(pool, n_pick))
  }
  
  w <- weights
  w[is.na(w)] <- 0
  pos_idx <- which(w > 0)
  
  if (length(pos_idx) >= n_pick) {
    return(sample(pool, n_pick, prob = w))
  } else {
    # Not enough high-prob cells? Pick all, then fill randomly
    picked_pos <- pool[pos_idx]
    remainder <- setdiff(pool, picked_pos)
    fill_n <- n_pick - length(picked_pos)
    return(c(picked_pos, sample(remainder, min(fill_n, length(remainder)))))
  }
}