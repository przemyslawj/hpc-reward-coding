library(testthat)

#############################################
context('Tests for mobility bouts detection')
#############################################

test_that('Returns two mobility bouts', {
  is_running = c(F,F,F, rep(TRUE, 4), F,F,F,F, rep(TRUE, 5), F)
  velocity = is_running %>% as.double()
  dist_reward = rep(0, length(velocity))
  atReward = c(rep(0, 8), 1, rep(0, 7), 2)
  res = mobility.bouts.tibble(11:(10+length(is_running)), is_running, atReward, velocity, dist_reward, window.len=3)
  expect_equal(c(4,12), res$index_start)
  expect_equal(c(14,22), res$timestamp_start)
  expect_equal(c(7, 16), res$index_end)
  expect_equal(c(17, 26), res$timestamp_end)
  expect_equal(c(1,2), res$rew_at_end)
})

test_that('Returns mobility bouts when NAs present', {
  is_running = c(F,F, NA, rep(TRUE, 4))
  atReward = c(rep(0, 7), 1)
  velocity = is_running %>% as.double()
  dist_reward = rep(0, length(velocity))
  res = mobility.bouts.tibble(11:(10+length(is_running)), is_running,atReward, velocity, dist_reward, window.len=3)
  expect_equal(c(4), res$index_start)
  expect_equal(c(14), res$timestamp_start)
  expect_equal(c(6), res$index_end)
  expect_equal(c(16), res$timestamp_end)
  expect_equal(c(1), res$rew_at_end)
})


#################################################
context('Tests for reward approaches detection')
#################################################
test_that('Returns two approaches to reward', {
  dist.reward0 = c(5, 5, 5, 4:1, 2:6, 6:3, 4)
  dist.reward1 = rep(10, length(dist.reward0))
  running.vals = c(F,F,F, rep(T, length(dist.reward0) - 3))
  velocity = running.vals %>% as.double()
  res = reward.approaches.tibble(1:length(running.vals), dist.reward0, dist.reward1, running.vals, velocity, 
                                 window.len=1, min.dist.run=1, min.approach.dist=3)
  expect_equal(c(4,8), res$index_start)
  expect_equal(c(7,16), res$index_end)
  expect_equal(c(4, 8), res$timestamp_start)
  expect_equal(c(7, 16), res$timestamp_end)
})
