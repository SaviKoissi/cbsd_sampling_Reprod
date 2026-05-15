# # surveillance_strategies.R
# 

strategy_baseline <- function(cells, n) {
  sample(cells, n)
}

strategy_host_density <- function(grid, cells, n) {
  probs <- grid$cassava_density[cells]
  probs <- probs / sum(probs, na.rm = TRUE)

  sample(cells, n, prob = probs)
}

update_states_after_detection <- function(state_status, detected_state) {
  state_status[detected_state] <- "detected"
  state_status
}

#========================================================
# SURVEILLANCE ENGINE (FIXED NA PROPAGATION)
#========================================================

run_surveillance <- function(sim, accessible_cells, grid_template,
                             strategy_fun,
                             n_surveys = 600,
                             params_field = list(n_plants = 100),
                             detection_accuracy = 1.0,
                             adjacency_list = NULL) {
  
  # REMOVE NA STATES EARLY (CRITICAL FIX)
  grid_template$state_id[is.na(grid_template$state_id)] <- "UNKNOWN"
  
  state_ids <- unique(grid_template$state_id)
  
  state_status <- setNames(rep("not_detected", length(state_ids)), state_ids)
  detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
  
  for (t in seq_along(sim)) {
    
    cat("Year:", t, "\n")
    grid <- sim[[t]]
    
    grid$state_id[is.na(grid$state_id)] <- "UNKNOWN"
    
    selected_cells <- strategy_fun(
      cells = accessible_cells,
      grid = grid,
      n = n_surveys
    )
    
    if (length(selected_cells) == 0) next
    
    detections <- sapply(selected_cells, function(i) {
      
      field <- simulate_field(
        time_since_infection = grid$infected_prop[i],
        params = params_field
      )
      
      detect_field(field, accuracy = detection_accuracy)
    })
    
    detected_cells <- selected_cells[detections]
    
    detected_states <- unique(grid$state_id[detected_cells])
    
    # SAFE UPDATE (NO NA INDEXING)
    detected_states <- detected_states[!is.na(detected_states)]
    
    state_status[detected_states] <- "detected"
    
    newly_detected <- setdiff(
      names(state_status)[state_status == "detected"],
      names(detection_year)[!is.na(detection_year)]
    )
    
    detection_year[newly_detected] <- t
    
    if (all(state_status == "detected")) break
  }
  
  detection_year
}