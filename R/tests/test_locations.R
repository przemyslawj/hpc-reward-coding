library(testthat)

context('Tests for reward locations')

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
                      exp_title='trial',
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
                      Valence=rep('Positive',6),
                      exp_title=c(rep('trial', 2), rep('beforetest',2), rep('trial',2)))
  loc.df$is_test = FALSE
  loc.df$is_test[3:4] = rep(TRUE,2)
  loc.df$date = format(loc.df$date)
  rew.df = add_location_set(loc.df)
  expect_equal(rew.df$location_set, c(1,1,1,1,2,2))
})

test_that('Sets location set when only one reward changes on a day', {
  d1 = as.Date('2019-01-01')
  loc.df = data.frame(Animal='A1',
                      Well_row=c(1,2,1,2,1,3,1,3,4,3,4,3),
                      Well_col=c(1,2,1,2,1,3,1,3,4,3,4,3),
                      date=c(d1,d1,d1+1,d1+1,d1+1,d1+1,d1+1,d1+1,d1+2,d1+2,d1+3,d1+3),
                      exp_title=c(rep('trial',2),rep('beforetest', 2),rep('trial',2), rep('aftertest', 2), rep('trial',4)),
                      is_test=FALSE,
                      Valence='Positive')
  loc.df$is_test[3:4] = TRUE
  loc.df$is_test[7:8] = TRUE
  loc.df$date = format(loc.df$date)
  rew.df = add_location_set(loc.df)
  expect_equal(rew.df$location_set, c(1,1,1,1,2,2,2,2,3,3,3,3))
  expect_equal(rew.df$location_ordinal, c(1,2,1,2,1,3,1,3,4,3,4,3))
})


context('Tests for previous reward locations')

test_that('Adds correct previous rewards when both rewards moved', {
  d1 = as.Date('2019-01-01')
  loc.df = data.frame(animal='A1',
                      Well_row=c(1,2,3,4),
                      Well_col=c(1,2,3,4),
                      date=c(d1,d1,d1+1,d1+1),
                      is_test=FALSE,
                      exp_title='trial',
                      Valence='Positive',
                      location_set=c(1,1,2,2),
                      location_ordinal=c(1,2,3,4))
  prev.loc.df = add_prev_locations(loc.df)
  expect_equal(4, nrow(subset(prev.loc.df, current_loc)))
  expect_equal(2, nrow(subset(prev.loc.df, !current_loc)))
  expect_equal(prev.loc.df$location_set, c(1, 1, 2, 2, 2, 2))
  expect_equal(prev.loc.df$location_ordinal, c(1, 2, 3, 4, 1, 2))
})

test_that('Adds correct previous rewards when one reward moved', {
  d1 = as.Date('2019-01-01')
  loc.df = data.frame(animal='A1',
                      Well_row=c(1,2,3,2),
                      Well_col=c(1,2,3,2),
                      date=c(d1,d1,d1+1,d1+1),
                      is_test=FALSE,
                      exp_title='trial',
                      Valence='Positive',
                      location_set=c(1,1,2,2),
                      location_ordinal=c(1,2,3,2))
  prev.loc.df = add_prev_locations(loc.df)
  expect_equal(4, nrow(subset(prev.loc.df, current_loc)))
  expect_equal(1, nrow(subset(prev.loc.df, !current_loc)))
  expect_equal(prev.loc.df$location_set, c(1, 1, 2, 2, 2, 2))
  expect_equal(prev.loc.df$location_ordinal, c(1, 2, 3, 4, 1, 2))
})

test_that('Adds correct two previous rewards after moving rewards twice', {
  d1 = as.Date('2019-01-01')
  loc.df = data.frame(animal='A1',
                      Well_row=c(1,2,3,2,3,4),
                      Well_col=c(1,2,3,2,3,4),
                      date=c(d1,d1,d1+1,d1+1,d1+2,d1+2),
                      is_test=FALSE,
                      exp_title='trial',
                      Valence='Positive',
                      location_set=c(1,1,2,2,3,3),
                      location_ordinal=c(1,2,3,2,3,4))
  prev.loc.df = add_prev_locations(loc.df, prev.loc.set.diff=1) %>%
    add_prev_locations(prev.loc.set.diff=2)
  expect_equal(6, nrow(subset(prev.loc.df, current_loc)))
  expect_equal(3, nrow(subset(prev.loc.df, !current_loc)))
  expect_equal(prev.loc.df$location_set, c(1, 1, 2, 2, 3, 3, 2, 3, 3))
  expect_equal(prev.loc.df$location_ordinal, c(1, 2, 3, 2, 3, 4, 1, 2, 1))
})
