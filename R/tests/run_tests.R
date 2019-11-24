library(testthat) 

source("tests/test_locations.R")
source("locations.R")

test_results <- test_dir("tests", reporter="summary")