#within_field_model.R

#================
# simulate_field <- function(time_since_infection, params) {
#   
#   n_plants <- params$n_plants
#   
#   # stronger amplification + floor
#   p_inf <- min(1, max(0.02, time_since_infection * 20))
#   
#   infected <- rbinom(n_plants, 1, prob = p_inf)
#   symptoms <- infected * rbinom(n_plants, 1, prob = 0.8)
#   
#   list(infected = infected, symptoms = symptoms)
# }

# simulate_field <- function(time_since_infection, params) {
#   
#   n_plants <- params$n_plants
#   
#   # FIX: bounded infection probability
#   p_inf <- plogis(time_since_infection * 10 - 2)
#   
#   infected <- rbinom(n_plants, 1, prob = p_inf)
#   
#   # symptoms depend on infection (NOT independent)
#   symptom_prob <- 0.8
#   symptoms <- rbinom(n_plants, 1, prob = symptom_prob) * infected
#   
#   list(
#     infected = infected,
#     symptoms = symptoms
#   )
# }
# within_field_model.R

simulate_field <- function(inf_prop, params) {
  n_plants <- params$n_plants
  # Boost p_inf so that if a cell is 'infected', the field is clearly infected
  p_inf <- min(1, inf_prop * 20) 
  
  infected <- rbinom(n_plants, 1, prob = p_inf)
  # High symptom visibility for testing
  symptoms <- infected * rbinom(n_plants, 1, prob = 0.9)
  
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