## Complet section 2

# --- 1. Timing Start ---
cat("JOB START: ", as.character(Sys.time()), "\n")
start_tic <- proc.time() # Start the "tic"

# --- 1. Setup ---
library(terra)
library(sf)
library(tidyverse)

# Load your data (Assuming paths are correct)
sim <- readRDS("outputs/simulation.rds")


# --- Adjacency List for Nigeria (Simplified example of 5 key states) ---
# In a full run, use: poly2nb(state_shapefile)
nigeria_adj <- list(
  "Nasarawa" = c("Plateau", "Kogi", "Benue", "Kaduna", "Taraba"),
  "Plateau"  = c("Nasarawa", "Bauchi", "Kaduna", "Taraba"),
  "Kebbi"    = c("Sokoto", "Zamfara", "Niger"),
  "Ogun"     = c("Lagos", "Oyo", "Osun", "Ondo"),
  "Anambra"  = c("Delta", "Kogi", "Enugu", "Imo", "Abia")
  # Add other states as needed or load from a spatial object
)

safe_sample <- function(cells, n, weights = NULL) {
  if (length(cells) == 0) return(integer(0))
  n_to_pick <- min(n, length(cells))
  
  # Case 1: No weights provided, just random sample
  if (is.null(weights)) {
    return(sample(cells, n_to_pick))
  }
  
  # Handle NAs in weights
  weights[is.na(weights)] <- 0
  
  # Case 2: If we have some positive weights but fewer than 'n'
  # We pick all positive weight cells first, then fill with the rest
  pos_weight_indices <- which(weights > 0)
  
  if (length(pos_weight_indices) >= n_to_pick) {
    # We have enough cassava cells to satisfy the whole survey
    return(sample(cells, n_to_pick, prob = weights))
  } else {
    # Not enough cassava cells! Pick all of them, then fill the rest randomly
    picked_pos <- cells[pos_weight_indices]
    remaining_cells <- cells[-pos_weight_indices]
    
    if (length(remaining_cells) == 0) return(picked_pos)
    
    fill_n <- n_to_pick - length(picked_pos)
    picked_random <- sample(remaining_cells, min(fill_n, length(remaining_cells)))
    
    return(c(picked_pos, picked_random))
  }
}

run_surveillance_full <- function(sim, accessible_cells, grid_template, 
                                  strategy_name, n_surveys = 2000, 
                                  adj_list = nigeria_adj) {
  
  # 1. Initialize
  state_ids <- unique(grid_template$state_id)
  state_ids <- state_ids[state_ids != "UNKNOWN" & !is.na(state_ids)]
  
  state_status <- setNames(rep("not_detected", length(state_ids)), state_ids)
  detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
  
  # 2. Yearly Loop
  for (t in seq_along(sim)) {
    grid_year <- sim[[t]]
    
    # --- FILTER ACCESSIBLE POOL ---
    # Strategies that STOP surveying once a state is positive
    needs_filtering <- grepl("single_detection", strategy_name) || 
      strategy_name == "survey_adjacent" || 
      strategy_name == "host_single_detection"
    
    if (needs_filtering) {
      pos_states <- names(state_status)[state_status == "detected"]
      current_pool <- accessible_cells[!(grid_template$state_id[accessible_cells] %in% pos_states)]
    } else {
      current_pool <- accessible_cells
    }
    
    if (length(current_pool) == 0) break # All target states detected
    
    # --- STRATEGY LOGIC ---
    selected <- integer(0)
    
    if (strategy_name == "baseline" || strategy_name == "single_detection") {
      selected <- safe_sample(current_pool, n_surveys)
      
    } else if (strategy_name == "host_density" || strategy_name == "host_single_detection") {
      w <- grid_year$cassava_density[current_pool]
      selected <- safe_sample(current_pool, n_surveys, weights = w)
      
    } else if (strategy_name == "survey_adjacent") {
      # Custom logic for Adjacent
      det_states <- names(state_status)[state_status == "detected"]
      if (length(det_states) == 0) {
        selected <- safe_sample(current_pool, n_surveys)
      } else {
        neighbors <- unique(unlist(adj_list[det_states]))
        target_states <- setdiff(neighbors, det_states)
        adj_pool <- current_pool[grid_template$state_id[current_pool] %in% target_states]
        # If no neighbors available, the paper implies continuing baseline in negative states
        if(length(adj_pool) == 0) adj_pool <- current_pool 
        selected <- safe_sample(adj_pool, n_surveys)
      }
    }
    
    if (length(selected) == 0) next
    
    # 3. Detection Step
    detections <- sapply(selected, function(i) {
      p_val <- grid_year$infected_prop[i]
      if (is.na(p_val) || p_val <= 0.01) return(FALSE) 
      
      field <- simulate_field(inf_prop = p_val, params = list(n_plants = 100))
      detect_field(field, accuracy = 1.0) 
    })
    
    # 4. Update
    detected_cells <- selected[which(detections)]
    if (length(detected_cells) > 0) {
      new_found <- unique(grid_template$state_id[detected_cells])
      new_found <- new_found[new_found != "UNKNOWN" & is.na(detection_year[new_found])]
      if(length(new_found) > 0) {
        state_status[new_found] <- "detected"
        detection_year[new_found] <- t
      }
    }
  }
  return(detection_year)
}



# Define strategies
strategies <- c("baseline", "single_detection", "survey_adjacent", 
                "host_density", "host_single_detection")

all_results_list <- list()

for (s in strategies) {
  cat("Processing Strategy:", s, "...\n")
  all_results_list[[s]] <- run_surveillance_full(
    sim = sim, 
    accessible_cells = valid_accessible, 
    grid_template = grid_clean, 
    strategy_name = s,
    n_surveys = 20000
  )
}


# --- 3. Save Output ---
results_final <- as.data.frame(all_results_list)
results_final$State <- rownames(results_final)
results_final <- results_final %>% select(State, everything())

# Save as CSV for easy inspection and RDS for R-specific metadata
write.csv(results_final, "outputs/surveillance_results.csv", row.names = FALSE)
saveRDS(results_final, "outputs/surveillance_results.rds")

message("Pipeline complete. Results saved to outputs/surveillance_results.csv")

# --- 4. Timing End ---
total_toc <- proc.time() - start_tic
cat("\n--------------------------------------------------\n")
cat("JOB SUCCESSFULLY COMPLETED\n")
cat("Total time elapsed:", round(total_toc["elapsed"] / 60, 2), "minutes\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("--------------------------------------------------\n")