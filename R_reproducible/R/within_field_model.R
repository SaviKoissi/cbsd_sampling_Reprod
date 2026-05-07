#within_field_model.R

#================
simulate_field <- function(time_since_infection, params) {
  
  n_plants <- params$n_plants
  
  # stronger amplification + floor
  p_inf <- min(1, max(0.02, time_since_infection * 20))
  
  infected <- rbinom(n_plants, 1, prob = p_inf)
  symptoms <- infected * rbinom(n_plants, 1, prob = 0.8)
  
  list(infected = infected, symptoms = symptoms)
}

#==================


detect_field <- function(field, accuracy = 1.0) {
  
  n_sample <- min(30, length(field$symptoms))
  sampled <- sample(seq_along(field$symptoms), n_sample)
  
  any(
    field$symptoms[sampled] == 1 &
      runif(n_sample) < accuracy
  )
}