# metric.R

compute_delay <- function(infection_time, detection_time) {
  detection_time - infection_time
}

compute_epidemic_size <- function(grid) {
  sum(grid$infected_prop > 0)
}

compute_detection_threshold <- function(results) {
  # Placeholder: requires definition from paper
  NA
}