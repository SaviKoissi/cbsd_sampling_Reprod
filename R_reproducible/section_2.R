# within_field_model.R (REFINED)
simulate_field <- function(inf_prop, params) {
  n_plants <- params$n_plants
  
  # Map grid prevalence (0-1) to a field-level probability. 
  # We use a steeper scaling: if 5% of the cell is infected, 
  # the specific field sampled likely has a high internal prevalence.
  p_inf <- min(1, inf_prop * 15) 
  
  infected <- rbinom(n_plants, 1, prob = p_inf)
  
  # Folier symptoms appear in 80% of infected plants
  symptoms <- infected * rbinom(n_plants, 1, prob = 0.8)
  
  list(infected = infected, symptoms = symptoms)
}

detect_field <- function(field, accuracy = 1.0) {
  # Standard survey: inspect 30 plants
  n_sample <- min(30, length(field$symptoms))
  sampled_indices <- sample(seq_along(field$symptoms), n_sample)
  
  # Detection occurs if ANY sampled plant is symptomatic AND correctly identified
  detected <- any(field$symptoms[sampled_indices] == 1 & 
                    runif(n_sample) < accuracy)
  return(detected)
}

# surveillance_logic.R (REFINED)
update_state_detection <- function(state_status, detected_cells, grid) {
  if (length(detected_cells) == 0) return(state_status)
  
  # Ensure we only look at valid IDs
  valid_detections <- detected_cells[!is.na(detected_cells)]
  detected_states <- unique(grid$state_id[valid_detections])
  
  # Filter out "UNKNOWN" or NA states
  detected_states <- detected_states[!is.na(detected_states) & detected_states != "UNKNOWN"]
  
  if (length(detected_states) > 0) {
    state_status[detected_states] <- "detected"
  }
  return(state_status)
}





#======

run_surveillance <- function(sim, accessible_cells, grid_template, 
                             strategy_fun, n_surveys = 2000, 
                             params_field = list(n_plants = 100),
                             detection_accuracy = 1.0, 
                             adjacency_list = NULL) {
  
  # 1. Initialize
  state_ids <- unique(grid_template$state_id)
  state_ids <- state_ids[state_ids != "UNKNOWN" & !is.na(state_ids)]
  
  state_status <- setNames(rep("not_detected", length(state_ids)), state_ids)
  detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
  
  # 2. Yearly Loop
  for (t in seq_along(sim)) {
    grid_year <- sim[[t]]
    
    # Selection Strategy
    if (identical(strategy_fun, strategy_baseline)) {
      selected <- strategy_baseline(accessible_cells, n_surveys)
    } else {
      # Pass required args for complex strategies
      selected <- strategy_fun(cells = accessible_cells, grid = grid_year, n = n_surveys)
    }
    
    if (length(selected) == 0) next
    
    # 3. Detection Step
    detections <- sapply(selected, function(i) {
      p_val <- grid_year$infected_prop[i]
      
      # Safety Check: If NA or 0, no detection
      if (is.na(p_val) || p_val <= 0) return(FALSE)
      
      field <- simulate_field(inf_prop = p_val, params = params_field)
      detect_field(field, accuracy = detection_accuracy)
    })
    
    detected_cells <- selected[which(detections)]
    
    # 4. Update Status
    if (length(detected_cells) > 0) {
      state_status <- update_state_detection(state_status, detected_cells, grid_template)
      
      # Record year for newly detected states
      newly_found <- names(state_status)[state_status == "detected" & is.na(detection_year)]
      detection_year[newly_found] <- t
    }
    
    # Early Exit if all found
    if (all(state_status == "detected")) break
  }
  
  return(detection_year)
}

#===============

# --- 1. Setup ---
library(terra)
library(sf)

# Load your data (Assuming paths are correct)
sim <- readRDS("outputs/simulation.rds")
# Ensure grid is cleaned of UNKNOWNs to focus on target states
grid_clean <- grid[grid$state_id != "UNKNOWN", ]
valid_accessible <- accessible_cells[accessible_cells %in% grid_clean$row_id]

# --- 2. Define Baseline Strategy ---
strategy_baseline <- function(cells, n) {
  if (length(cells) <= n) return(cells)
  sample(cells, n)
}

# --- 3. Define Host Density Strategy ---
strategy_host_density <- function(cells, grid, n) {
  probs <- grid$cassava_density[cells]
  probs[is.na(probs)] <- 0
  if (sum(probs) == 0) return(sample(cells, n))
  sample(cells, n, prob = probs / sum(probs))
}

# --- 4. Run Pipeline ---

cat("Starting Surveillance Simulation...\n")

res_baseline <- run_surveillance(
  sim = sim,
  accessible_cells = valid_accessible,
  grid_template = grid_clean,
  strategy_fun = strategy_baseline,
  n_surveys = 2000
)

res_host <- run_surveillance(
  sim = sim,
  accessible_cells = valid_accessible,
  grid_template = grid_clean,
  strategy_fun = strategy_host_density,
  n_surveys = 2000
)

# --- 5. Compare Results ---
results_df <- data.frame(
  State = names(res_baseline),
  Baseline_Year = res_baseline,
  HostDensity_Year = res_host
)

print(results_df)

# Final Verification: How many states were NEVER detected?
missing_baseline <- sum(is.na(res_baseline))
cat("States not detected (Baseline):", missing_baseline, "out of", length(res_baseline), "\n")


infection_trend <- sapply(sim, function(df) sum(df$infected_prop > 0.01))
plot(infection_trend, type = "l", main = "Number of Infected Cells over 30 Years", 
     xlab = "Year (List Index)", ylab = "Count of Infected Cells")
