## Sampling

#samplin.R 

#--------------------------------------------
# Select cells within 1 km of roads (FIXED)
#--------------------------------------------


select_accessible_cells <- function(grid_sf, roads) {
  
  within_dist <- sf::st_is_within_distance(grid_sf, roads, dist = 1000)
  accessible <- lengths(within_dist) > 0
  
  which(accessible)
}

#--------------------------------------------
# Subsampling fields for variability (OK)
#--------------------------------------------
resample_fields <- function(fields_pool, n_fields, n_rep = 500) {
  replicate(n_rep, {
    sample(fields_pool, size = n_fields, replace = TRUE)
  }, simplify = FALSE)
}

#--------------------------------------------
# Safe sampling helper (NEW - IMPORTANT)
#--------------------------------------------
safe_sample <- function(x, n, prob = NULL) {
  
  if (length(x) == 0) return(integer(0))
  
  if (length(x) <= n) {
    return(x)
  } else {
    return(sample(x, n, prob = prob))
  }
}

#--------------------------------------------
# Helper: filter eligible cells by state
#--------------------------------------------
filter_cells_by_state <- function(cells, grid, allowed_states) {
  cells[grid$state_id[cells] %in% allowed_states]
}

#--------------------------------------------
# 1. Baseline
#--------------------------------------------
strategy_baseline <- function(cells, n) {
  safe_sample(cells, n)
}


#--------------------------------------------
# 2. Single-detection
#--------------------------------------------
strategy_single_detection <- function(cells, grid, n, state_status) {
  
  allowed_states <- names(state_status)[state_status != "detected"]
  eligible_cells <- filter_cells_by_state(cells, grid, allowed_states)
  
  safe_sample(eligible_cells, n)
}

#--------------------------------------------
# 3. Survey-adjacent
#--------------------------------------------
strategy_survey_adjacent <- function(cells, grid, n, state_status, adjacency_list) {
  
  infected_states <- names(state_status)[state_status == "detected"]
  not_detected_states <- names(state_status)[state_status == "not_detected"]
  
  adjacent_states <- unique(unlist(adjacency_list[infected_states]))
  
  allowed_states <- intersect(not_detected_states, adjacent_states)
  eligible_cells <- filter_cells_by_state(cells, grid, allowed_states)
  
  # ⚠️ fallback (important)
  if (length(eligible_cells) == 0) {
    return(safe_sample(cells, n))  # fallback to baseline
  }
  
  safe_sample(eligible_cells, n)
}

#--------------------------------------------
# 4. Host-density
#--------------------------------------------
# strategy_host_density <- function(cells, grid, n) {
#   
#   probs <- grid$cassava_density[cells]
#   probs <- probs / sum(probs, na.rm = TRUE)
#   
#   safe_sample(cells, n, prob = probs)
# }

strategy_host_density <- function(cells, grid, n) {
  
  probs <- grid$cassava_density[cells]
  probs[is.na(probs)] <- 0
  
  if (sum(probs) == 0) {
    return(safe_sample(cells, n))
  }
  
  probs <- probs / sum(probs)
  
  safe_sample(cells, n, prob = probs)
}

#--------------------------------------------
# 5. Host-density + single-detection
#--------------------------------------------
# strategy_host_density_single_detection <- function(cells, grid, n, state_status) {
#   
#   allowed_states <- names(state_status)[state_status != "detected"]
#   eligible_cells <- filter_cells_by_state(cells, grid, allowed_states)
#   
#   if (length(eligible_cells) == 0) return(integer(0))
#   
#   probs <- grid$cassava_density[eligible_cells]
#   probs <- probs / sum(probs, na.rm = TRUE)
#   
#   safe_sample(eligible_cells, n, prob = probs)
# }

strategy_host_density_single_detection <- function(cells, grid, n, state_status) {
  
  allowed_states <- names(state_status)[state_status != "detected"]
  eligible_cells <- filter_cells_by_state(cells, grid, allowed_states)
  
  if (length(eligible_cells) == 0) {
    return(safe_sample(cells, n))
  }
  
  probs <- grid$cassava_density[eligible_cells]
  probs[is.na(probs)] <- 0
  
  if (sum(probs) == 0) {
    return(safe_sample(eligible_cells, n))
  }
  
  probs <- probs / sum(probs)
  
  safe_sample(eligible_cells, n, prob = probs)
}

#--------------------------------------------
# Update detection status
#--------------------------------------------
update_state_detection <- function(state_status, detected_cells, grid) {
  
  if (length(detected_cells) == 0) return(state_status)
  
  detected_states <- unique(grid$state_id[detected_cells])
  state_status[detected_states] <- "detected"
  
  state_status
}