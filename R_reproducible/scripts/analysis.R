# analys.R

source("R/metrics.R")

sim <- readRDS("outputs/simulation.rds")

sizes <- sapply(sim, compute_epidemic_size)

plot(sizes, type = "l", main = "Epidemic size over time")
