# surveillance_strategies.R

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