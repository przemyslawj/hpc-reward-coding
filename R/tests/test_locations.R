library(testthat)

context('Tests for setting up in reward locations')

test_that('Returns empty df when no rewards', {
  loc.df = add_location_set(data.frame())
  expect_equal(nrow(loc.df), 0)
})

test_that('Sets location set when both rewards changed on the same day', {
  d1=as.Date('2019-01-01')
  loc.df = data.frame(Animal=rep('A1', 6),
                      Well_row=c(1,2,1,2,3,4),
                      Well_col=c(1,2,1,2,3,4),
                      date=c(d1,d1,d1+1,d1+1,d1+2,d1+2),
                      is_test=rep(FALSE,6),
                      Valence=rep('Positive',6))
  loc.df$date = format(loc.df$date)
  rew.df = add_location_set(loc.df)
  expect_equal(rew.df$location_set, c(1,1,1,1,2,2))
  expect_equal(rew.df$location_ordinal, c(1,2,1,2,3,4))
})

test_that('Sets location set when the test trial is present', {
  d1=as.Date('2019-01-01')
  loc.df = data.frame(Animal=rep('A1', 6),
                      Well_row=c(1,2,1,2,3,4),
                      Well_col=c(1,2,1,2,3,4),
                      date=c(d1,d1,d1+1,d1+1,d1+1,d1+1),
                      Valence=rep('Positive',6))
  loc.df$is_test = FALSE
  loc.df$is_test[3:4] = rep(TRUE,2)
  loc.df$date = format(loc.df$date)
  rew.df = add_location_set(loc.df)
  expect_equal(rew.df$location_set, c(1,1,1,1,2,2))
})

test_that('Sets location set when only one reward changes on a day', {
  d1=as.Date('2019-01-01')
  loc.df = data.frame(Animal=rep('A1', 8),
                      Well_row=c(1,2,1,3,4,3,4,3),
                      Well_col=c(1,2,1,3,4,3,4,3),
                      date=c(d1,d1,d1+1,d1+1,d1+2,d1+2,d1+3,d1+3),
                      Valence=rep('Positive',8))
  loc.df$is_test = FALSE
  loc.df$date = format(loc.df$date)
  rew.df = add_location_set(loc.df)
  expect_equal(rew.df$location_set, c(1,1,2,2,3,3,3,3))
  expect_equal(rew.df$location_ordinal, c(1,2,1,3,4,3,4,3))
})