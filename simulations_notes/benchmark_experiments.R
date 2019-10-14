library(microbenchmark)
source("simulations_notes/old_version.R")
source("simulations_notes/compromise_version.R")
source("simulations_notes/simulation.R")

metrics <- list("mse" = mse, "bias" = bias, "new" = new_estimator)

benchmark_results <- microbenchmark(
    old_version = behaviour(),
    compromised_version = posrednia(100),
    modified_version2 = simulation2(100, 100, 0.1, 0.5, 5),
    modified_verision1 = simulation(100, 100, 0.1, 0.5, 5, metrics),
    modified_verision = simulation3(100, 100, 0.1, 0.5, 5, metrics),
    times = 500
  )

benchmark_results


