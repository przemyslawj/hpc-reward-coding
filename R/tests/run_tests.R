library(testthat) 

source("locations.R")
source("tests/test_locations.R")

source("place_field_utils.R")
source("tests/test_place_field_utils.R")

test_results <- test_dir("tests", reporter="summary")