library(testthat) 

source("locations.R")
source("tests/test_locations.R")
source("bayesian_decoder.R")
Rcpp::sourceCpp("bayesian_decoder.cpp")
source("tests/test_bayesian_decoder.R")

source("trace_utils.R")
source("tests/test_trace_utils.R")

Rcpp::sourceCpp("palce_field.cpp")
source("tests/test_place_field.R")

test_results <- test_dir("tests", reporter="summary")