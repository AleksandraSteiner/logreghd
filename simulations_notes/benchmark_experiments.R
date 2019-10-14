library(microbenchmark)
source("simulations_notes/old_version.R")
source("simulations_notes/compromise_version.R")
source("simulations_notes/simulation.R")


benchmark_results <- microbenchmark(
    old_version = behaviour(),
    compromised_version = posrednia(10),
    modified_version = simulation2(10, 100, 0.1, 0.5, 5),
    times = 10
  )

benchmark_results
