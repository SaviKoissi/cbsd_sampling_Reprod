## Complet section 2

# --- 1. Timing Start ---
cat("JOB START: ", as.character(Sys.time()), "\n")
start_tic <- proc.time() 

# --- 2. Setup Libraries ---
library(terra)
library(dplyr)

base_dir <- "/home/savi/project/SurveyCBSD/SurveyCBSD/R_reproducible"

#------------------------------------
# LOAD DATA (Matrix-backed recovery)
#------------------------------------
print("Loading raw numeric matrix inputs...")
cassava_mat  <- readRDS(file.path(base_dir, "data/cassavaMap/cassava_nigeria.rds"))
whitefly_mat <- readRDS(file.path(base_dir, "data/whitefly/whitefly.rds"))

# Reconstruct memory rasters and explicitly lock spatial reference frames
cassava <- terra::rast(cassava_mat)
terra::crs(cassava) <- "EPSG:4326"
terra::ext(cassava) <- c(2.6685, 14.6788, 4.2730, 13.8944)

whitefly <- terra::rast(whitefly_mat)
terra::crs(whitefly) <- "EPSG:4326"
terra::ext(whitefly) <- c(2.6685, 14.6788, 4.2730, 13.8944)

print("Loading vector spatial boundaries...")
states_vect <- terra::vect(file.path(base_dir, "data/state_boundaries/nga_admin1.shp"))
roads_vect  <- terra::vect(file.path(base_dir, "data/road_network/NGA_roads.shp"))

# Force coordinate structures on the vector objects in memory
terra::crs(states_vect) <- "EPSG:4326"
terra::crs(roads_vect)  <- "EPSG:4326"

#------------------------------------
# ALIGNMENT 
#------------------------------------
whitefly_aligned <- terra::resample(whitefly, cassava, method = "bilinear")

#-----------------------------------
# create_grid Function
#-----------------------------------
create_grid <- function(cassava, whitefly) {
  if (!terra::compareGeom(cassava, whitefly)) {
    stop("Rasters are not aligned")
  }
  
  cassava_vals  <- terra::values(cassava)
  whitefly_vals <- terra::values(whitefly)
  
  valid <- !is.na(cassava_vals) & !is.na(whitefly_vals) & (cassava_vals > 0)
  
  data.frame(
    row_id           = seq_len(sum(valid)),
    cell_id          = which(valid),
    cassava_density  = cassava_vals[valid],
    whitefly_density = whitefly_vals[valid],
    infected_prop    = 0
  )
}

grid <- create_grid(cassava, whitefly_aligned)

#-----------------------------------
# attach_states Function (Pure Terra Engine)
#-----------------------------------
attach_states <- function(grid, raster_ref, vec_states) {
  print("DEBUG: Executing zero-sf state association...")
  
  coords <- terra::xyFromCell(raster_ref, grid$cell_id)
  pts    <- terra::vect(coords, crs = terra::crs(vec_states))
  
  joined_data <- terra::extract(vec_states, pts)
  col_names   <- names(joined_data)
  
  state_col <- intersect(col_names, c("adm1_name", "ADM1_NAME", "name_1", "NAME_1", "NOM_1"))[1]
  if (is.na(state_col)) {
    stop("Could not find state name column. Columns found: ", paste(col_names, collapse=", "))
  }
  
  state_vec <- as.character(joined_data[[state_col]])
  state_vec[is.na(state_vec)] <- "UNKNOWN"
  
  grid$state_id <- state_vec
  return(grid)
}

grid <- attach_states(grid, cassava, states_vect)
stopifnot(!all(grid$state_id == "UNKNOWN"))

# Remove cells falling completely outside state boundaries
grid <- grid[grid$state_id != "UNKNOWN", ]

#------------------------------------
# ACCESSIBILITY (Zero-SF Distance Buffer)
#------------------------------------
print("Computing accessibility metrics via planar distances...")
coords <- terra::xyFromCell(cassava, grid$cell_id)
grid_pts <- terra::vect(coords, crs = "EPSG:4326")

# Using high-precision geodesic distance processing to bypass st_transform projection blocks
distances <- terra::distance(grid_pts, roads_vect)
accessible_indices <- which(distances <= 5000)

accessible_cells <- grid$cell_id[grid$cell_id %in% accessible_indices]

#------------------------------------
# DATA LOADING & SIMULATION CHECK
#------------------------------------
print("Loading epidemiological projection tensors...")
sim <- readRDS(file.path(base_dir, "outputs/simulation.rds"))

sim <- lapply(sim, function(df) {
  if("infected_prop" %in% names(df)) {
    df$infected_prop[is.na(df$infected_prop)] <- 0
  } else {
    df$infected_prop <- 0
  }
  return(df)
})

# --- Adjacency Matrix ---
nigeria_adj <- list(
  "Nasarawa" = c("Plateau", "Kogi", "Benue", "Kaduna", "Taraba"),
  "Plateau"  = c("Nasarawa", "Bauchi", "Kaduna", "Taraba"),
  "Kebbi"    = c("Sokoto", "Zamfara", "Niger"),
  "Ogun"     = c("Lagos", "Oyo", "Osun", "Ondo"),
  "Anambra"  = c("Delta", "Kogi", "Enugu", "Imo", "Abia")
)

# --- Sampling Function ---
safe_sample <- function(cells, n, weights = NULL) {
  if (length(cells) == 0) return(integer(0))
  n_to_pick <- min(n, length(cells))
  
  if (is.null(weights)) {
    return(sample(cells, n_to_pick))
  }
  
  weights[is.na(weights)] <- 0
  pos_weight_indices <- which(weights > 0)
  
  if (length(pos_weight_indices) >= n_to_pick) {
    return(sample(cells, n_to_pick, prob = weights))
  } else {
    picked_pos <- cells[pos_weight_indices]
    remaining_cells <- cells[-pos_weight_indices]
    if (length(remaining_cells) == 0) return(picked_pos)
    
    fill_n <- n_to_pick - length(picked_pos)
    picked_random <- sample(remaining_cells, min(fill_n, length(remaining_cells)))
    return(c(picked_pos, picked_random))
  }
}

simulate_field <- function(inf_prop, params) { return(inf_prop) }
detect_field   <- function(field, accuracy) { return(runif(1) < field) }

# --- Surveillance Loop Execution ---
run_surveillance_full <- function(sim, accessible_pool, grid_template, strategy_name, n_surveys = 2000, adj_list = nigeria_adj) {
  state_ids <- unique(grid_template$state_id)
  state_ids <- state_ids[state_ids != "UNKNOWN" & !is.na(state_ids)]
  
  state_status   <- setNames(rep("not_detected", length(state_ids)), state_ids)
  detection_year <- setNames(rep(NA, length(state_ids)), state_ids)
  
  for (t in seq_along(sim)) {
    grid_year <- sim[[t]]
    
    if (nrow(grid_year) != nrow(grid_template)) next
    
    needs_filtering <- grepl("single_detection", strategy_name) || 
                       strategy_name == "survey_adjacent" || 
                       strategy_name == "host_single_detection"
    
    if (needs_filtering) {
      pos_states   <- names(state_status)[state_status == "detected"]
      current_pool <- accessible_pool[!(grid_template$state_id[grid_template$cell_id %in% accessible_pool] %in% pos_states)]
    } else {
      current_pool <- accessible_pool
    }
    
    if (length(current_pool) == 0) break
    selected <- integer(0)
    
    if (strategy_name == "baseline" || strategy_name == "single_detection") {
      selected <- safe_sample(current_pool, n_surveys)
    } else if (strategy_name == "host_density" || strategy_name == "host_single_detection") {
      w <- grid_year$cassava_density[grid_template$cell_id %in% current_pool]
      selected <- safe_sample(current_pool, n_surveys, weights = w)
    } else if (strategy_name == "survey_adjacent") {
      det_states <- names(state_status)[state_status == "detected"]
      if (length(det_states) == 0) {
        selected <- safe_sample(current_pool, n_surveys)
      } else {
        neighbors     <- unique(unlist(adj_list[det_states]))
        target_states <- setdiff(neighbors, det_states)
        adj_pool      <- current_pool[grid_template$state_id[grid_template$cell_id %in% current_pool] %in% target_states]
        if(length(adj_pool) == 0) adj_pool <- current_pool 
        selected <- safe_sample(adj_pool, n_surveys)
      }
    }
    
    if (length(selected) == 0) next
    
    detections <- sapply(selected, function(cell) {
      idx   <- which(grid_template$cell_id == cell)
      p_val <- grid_year$infected_prop[idx]
      if (is.na(p_val) || p_val <= 0.01) return(FALSE) 
      field <- simulate_field(inf_prop = p_val, params = list(n_plants = 100))
      detect_field(field, accuracy = 1.0) 
    })
    
    detected_cells <- selected[which(detections)]
    if (length(detected_cells) > 0) {
      new_found <- unique(grid_template$state_id[grid_template$cell_id %in% detected_cells])
      new_found <- new_found[new_found != "UNKNOWN" & is.na(detection_year[new_found])]
      if(length(new_found) > 0) {
        state_status[new_found]   <- "detected"
        detection_year[new_found] <- t
      }
    }
  }
  return(detection_year)
}

# Execute strategies
strategies <- c("baseline", "single_detection", "survey_adjacent", "host_density", "host_single_detection")
all_results_list <- list()

for (s in strategies) {
  cat("Processing Strategy:", s, "...\n")
  all_results_list[[s]] <- run_surveillance_full(
    sim             = sim, 
    accessible_pool = accessible_cells, 
    grid_template   = grid, 
    strategy_name   = s,
    n_surveys       = 2000
  )
}

# --- 3. Save Output ---
results_final       <- as.data.frame(all_results_list)
results_final$State <- rownames(results_final)
results_final       <- results_final %>% select(State, everything())

write.csv(results_final, file.path(base_dir, "outputs/surveillance_results.csv"), row.names = FALSE)
saveRDS(results_final, file.path(base_dir, "outputs/surveillance_results.rds"))

message("Pipeline complete. Results saved to outputs/surveillance_results.csv")

# --- 4. Timing End ---
total_toc <- proc.time() - start_tic
cat("\n--------------------------------------------------\n")
cat("JOB SUCCESSFULLY COMPLETED\n")
cat("Total time elapsed:", round(total_toc["elapsed"] / 60, 2), "minutes\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("--------------------------------------------------\n")
