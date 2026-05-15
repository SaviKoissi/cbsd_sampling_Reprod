# pipeline.R

run_surveillance_pipeline <- function(
    host_raster,          # terra SpatRaster (e.g. Cassava)
    vector_raster,        # terra SpatRaster (e.g. Whitefly)
    admin_boundaries,     # sf object (States/Provinces)
    init_params = list(alpha=2, beta=1, max_dist=10000, years=20, n_surveys=2000),
    strategy = "baseline"
) {
  
  # 1. Align Environments
  template <- host_raster
  vector_res <- resample(vector_raster, template)
  
  # 2. Create Grid
  v_host <- values(host_raster)
  v_vect <- values(vector_res)
  valid <- !is.na(v_host) & !is.na(v_vect)
  
  grid <- data.frame(
    cell_id = which(valid),
    host_density = v_host[valid],
    vector_density = v_vect[valid],
    infected_prop = 0,
    admin_id = terra::extract(admin_boundaries, crs(template), xy=TRUE) # Dynamic Admin Assignment
  )
  
  # 3. Initialize Epidemic (Random start based on host density)
  start_node <- run_sampler(1:nrow(grid), 1, weights = grid$host_density)
  grid$infected_prop[start_node] <- 0.1
  
  coords <- crds(template)[grid$cell_id, ]
  results <- list()
  admin_status <- setNames(rep(NA, nrow(admin_boundaries)), admin_boundaries$name) # Flexible name col
  
  # 4. Simulation + Surveillance Loop
  for (t in 1:init_params$years) {
    # A. Pathogen Spread
    pressure <- compute_pressure(grid, coords, init_params$alpha, init_params$beta, init_params$max_dist)
    # Normalize and update
    pressure <- (pressure / max(pressure, 1e-6)) * (1 - grid$infected_prop)
    grid$infected_prop <- pmin(grid$infected_prop + (rbinom(nrow(grid), 1, pressure) * 0.1), 1)
    
    # B. Surveillance Selection
    # (Apply logic based on 'strategy' string)
    # [Insert Strategy Switch Logic Here using run_sampler]
    
    # C. Record detections
    # ... Detection logic ...
  }
  
  return(results)
}