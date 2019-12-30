library(testthat) 

source("locations.R")
source("tests/test_locations.R")
source("bayesian_decoder.R")
source("tests/test_bayesian_decoder.R")

test_results <- test_dir("tests", reporter="summary")