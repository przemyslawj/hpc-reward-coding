library(testthat)

context('Tests for min distance to reward')

test_that('Returns NAs when no rewards', {
  locs = data.frame()
  field.x = c(1, 50)
  field.y = c(1, 50)
  rew.distances = calc.min.rew.dist(locs, field.x, field.y)
  expect_equal(rew.distances$location.ordinal, c(NA, NA))
  expect_equal(rew.distances$rew.dist, c(NA, NA))
  expect_equal(rew.distances$rew.angle, c(NA, NA))
})

test_that('Returns distances to closer of two rewards', {
  locs = data.frame(location_ordinal = c(1, 2),
                    trans_x = c(1, 46),
                    trans_y = c(2, 47))
  field.x = c(1, 50, 1)
  field.y = c(1, 50, 2)
  rew.distances = calc.min.rew.dist(locs, field.x, field.y)
  expect_equal(rew.distances$location.ordinal, c(1, 2, 1))
  expect_equal(rew.distances$rew.dist, c(1, 5, 0))
  expect_equal(rew.distances$rew.angle, c(90.0, 216.8699, 0))
})