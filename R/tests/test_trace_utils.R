library(testthat)
library(dplyr)

context('Tests for the filter.running')

test_that('filter.running succeeds for one trial', {
  running.vel.thr = 2
  epoch.dur.ms = 100
  fr = 20
  running.pos = seq(0, by=running.vel.thr*8 / fr, length.out=4)
  df = data.frame(
    id=1:20,
    trial_id = '1',
    cell_id = 1,
    exp_title = 'trial',
    timestamp = seq(0, by=1000/fr, length.out=20),
    smooth_trans_x = c(running.pos, rep(4.0, 6), running.pos + 2.0, rep(10.0, 6)),
    smooth_trans_y = 1,
    velocity = c(rep(running.vel.thr * 4, 5), rep(0, 5), rep(running.vel.thr * 4, 5), rep(0, 5))
  )
  running.index = isRunning(df, 2, running.vel.thr, epoch.dur.ms)
  #df.filtered = filter.running(df, min.run.velocity=2, mean.run.velocity=running.vel.thr, window.dur.ms=epoch.dur.ms)
  df.filtered = df[which(running.index),]
  expect_equal(nrow(df.filtered), 10)
  expect_equal(df.filtered$id, c(1:5, 11:15))
}) 