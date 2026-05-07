#run_surveillance.R

# # sampling estimation
# n_surveys = 2000
# run_surveillance <- function(sim, accessible_cells, grid_template,
#                              strategy_fun,
#                              n_surveys = 600,
#                              params_field = list(n_plants = 100),
#                              detection_accuracy = 1.0,
#                              adjacency_list = NULL) {
#   
#   state_ids <- unique(grid_template$state_id)
#   
#   state_status <- setNames(rep("not_detected", length(state_ids)), state_ids)
#   detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
#   
#   for (t in seq_along(sim)) {
#     
#     cat("Year:", t, "\n")
#     
#     grid <- sim[[t]]
#     
#     #-------------------------------
#     # Select cells depending on strategy
#     #-------------------------------
#     if (identical(strategy_fun, strategy_survey_adjacent)) {
#       selected_cells <- strategy_fun(
#         cells = accessible_cells,
#         grid = grid,
#         n = n_surveys,
#         state_status = state_status,
#         adjacency_list = adjacency_list
#       )
#       
#     } else if (identical(strategy_fun, strategy_single_detection) ||
#                identical(strategy_fun, strategy_host_density_single_detection)) {
#       
#       selected_cells <- strategy_fun(
#         cells = accessible_cells,
#         grid = grid,
#         n = n_surveys,
#         state_status = state_status
#       )
#       
#     } else if (identical(strategy_fun, strategy_host_density)) {
#       
#       selected_cells <- strategy_fun(
#         cells = accessible_cells,
#         grid = grid,
#         n = n_surveys
#       )
#       
#     } else {
#       
#       selected_cells <- strategy_fun(
#         cells = accessible_cells,
#         n = n_surveys
#       )
#     }
#     
#     if (length(selected_cells) == 0) next
#     
#     #-------------------------------
#     # Detection
#     #-------------------------------
#     detections <- sapply(selected_cells, function(i) {
#       
#       field <- simulate_field(
#         time_since_infection = grid$infected_prop[i],
#         params = params_field
#       )
#       
#       detect_field(field, accuracy = detection_accuracy)
#     })
#     
#     detected_cells <- selected_cells[detections == TRUE]
#     
#     #-------------------------------
#     # Update states
#     #-------------------------------
#     state_status <- update_state_detection(
#       state_status,
#       detected_cells,
#       grid
#     )
#     
#     newly_detected <- names(state_status)[
#       state_status == "detected" & is.na(detection_year)
#     ]
#     
#     detection_year[newly_detected] <- t
#     
#     if (all(state_status == "detected")) break
#   }
#   
#   return(detection_year)
# }

#========================================================
# RUN SURVEILLANCE SIMULATION
#========================================================

run_surveillance <- function(sim,
                             accessible_cells,
                             grid_template,
                             strategy_fun,
                             n_surveys = 600,
                             params_field = list(n_plants = 100),
                             detection_accuracy = 1.0,
                             adjacency_list = NULL) {
  
  #------------------------------------------------------
  # Initialize state tracking
  #------------------------------------------------------
  state_ids <- unique(grid_template$state_id)
  
  state_status <- setNames(rep("not_detected", length(state_ids)), state_ids)
  detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
  
  #------------------------------------------------------
  # Main yearly loop
  #------------------------------------------------------
  for (t in seq_along(sim)) {
    
    cat("Year:", t, "\n")
    
    grid <- sim[[t]]
    
    #--------------------------------------------------
    # Strategy selection (ALL assume row-index system)
    #--------------------------------------------------
    if (identical(strategy_fun, strategy_survey_adjacent)) {
      
      selected_cells <- strategy_fun(
        cells = accessible_cells,
        grid = grid,
        n = n_surveys,
        state_status = state_status,
        adjacency_list = adjacency_list
      )
      
    } else if (identical(strategy_fun, strategy_single_detection) ||
               identical(strategy_fun, strategy_host_density_single_detection)) {
      
      selected_cells <- strategy_fun(
        cells = accessible_cells,
        grid = grid,
        n = n_surveys,
        state_status = state_status
      )
      
    } else if (identical(strategy_fun, strategy_host_density)) {
      
      selected_cells <- strategy_fun(
        cells = accessible_cells,
        grid = grid,
        n = n_surveys
      )
      
    } else {
      
      selected_cells <- strategy_fun(
        cells = accessible_cells,
        n = n_surveys
      )
    }
    
    if (length(selected_cells) == 0) next
    
    #--------------------------------------------------
    # DETECTION STEP (FIXED)
    #--------------------------------------------------
    detections <- sapply(selected_cells, function(i) {
      
      # ❗ FIX:
      # infected_prop is NOT infection age.
      # We convert it into a proxy signal to prevent NA collapse.
      #
      # TEMPORARY assumption:
      # higher prevalence → longer infection duration
      
      #infection_proxy_time <- grid$infected_prop[i] * 12  # scale to months
      infection_proxy_time = sample(6:24, 1)
      
      field <- simulate_field(
        time_since_infection = infection_proxy_time,
        params = params_field
      )
      
      detect_field(field, accuracy = detection_accuracy)
    })
    
    detected_cells <- selected_cells[which(detections)]
    
    #--------------------------------------------------
    # UPDATE STATE DETECTIONS
    #--------------------------------------------------
    if (length(detected_cells) > 0) {
      
      state_status <- update_state_detection(
        state_status,
        detected_cells,
        grid
      )
      
      newly_detected <- names(state_status)[
        state_status == "detected" & is.na(detection_year)
      ]
      
      detection_year[newly_detected] <- t
    }
    
    #--------------------------------------------------
    # EARLY STOP CONDITION
    #--------------------------------------------------
    if (all(state_status == "detected")) break
  }
  
  return(detection_year)
}