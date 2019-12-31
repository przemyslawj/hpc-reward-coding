library(testthat)
library(dplyr)

context('Tests for the bayesian model')


test_that('to_1dim succeeds', {
  expect_equal(to_1dim(0, 0), 1)
  expect_equal(to_1dim(0, 1), 2)
  expect_equal(to_1dim(1, 0), 21)
}) 

test_that('from1dim succeeds', {
  expect_equal(from_1dim(1), list(x=0, y=0))
  expect_equal(from_1dim(2), list(x=0, y=1))
  expect_equal(from_1dim(21), list(x=1, y=0))
}) 

test_that('timebin.traces success for one trial many cells', {
  df = data.frame(
    animal='A',
    date='2019-01-01',
    cell_id=c(1,1,1,2,2,2),
    trial_id='1',
    trial=1,
    timestamp=c(1,2,11,1,2,12),
    trace=c(0,2,2,1,3,4),
    smooth_trans_x=1,
    smooth_trans_y=2
  )
  actual.binned = timebin.traces(df, timebin.dur.msec=10)
  actual.binned = actual.binned %>% arrange(cell_id, time_bin)
  expect_equal(nrow(actual.binned), 4)
  expect_equal(actual.binned$cell_id, c(1,1,2,2))
  expect_equal(actual.binned$time_bin, c(0,1,0,1))
  expect_equal(actual.binned$bin.trace, c(1,2,2,4))
  expect_equal(actual.binned$mean.x, rep(1,4))
  expect_equal(actual.binned$mean.y, rep(2,4))
  
})

test_that('timebin.traces success for two trials', {
  df = data.frame(
    animal='A',
    date='2019-01-01',
    cell_id=rep(1,6),
    trial_id=rep(c('trial1','trial2'), each=3),
    trial=rep(c(1,2), each=3),
    timestamp=c(1,2,21,1,2,12),
    trace=c(0,2,2,1,3,4),
    smooth_trans_x=1,
    smooth_trans_y=2
  )
  actual.binned = timebin.traces(df, timebin.dur.msec=10)
  actual.binned = actual.binned %>% arrange(cell_id, time_bin)
  expect_equal(nrow(actual.binned), 4)
  expect_equal(actual.binned$cell_id, rep(1,4))
  expect_equal(actual.binned$time_bin, c(0,2,3,4))
  expect_equal(actual.binned$mean.trace, c(1,2,2,4))
  expect_equal(actual.binned$mean.x, rep(1,4))
  expect_equal(actual.binned$mean.y, rep(2,4))
})

test_that("bin.responses success for two cells", {
  df = data.frame(
    animal='A',
    date='2019-01-01',
    cell_id=rep(c(1,2),each=4),
    trial_id='trial1',
    time_bin=rep(1:4),
    mean.trace=c(1:4, 9:6),
    mean.x=1:8,
    mean.y=1
  )
  df$cell_id = as.factor(df$cell_id)
  quantile.fractions=c(0.5, 1.0)
  actual.binned = bin.responses(df, quantile.fractions)
  actual.binned = actual.binned %>% arrange(cell_id, time_bin)
  expect_equal(nrow(actual.binned), 8)
  expect_equal(actual.binned$response_bin, c(1, 1, 2, 2, 2, 2, 1, 1))
})

test_that('create.model success for one cell', {
  df = data.frame(
    animal='A',
    date='2019-01-01',
    cell_id=1,
    time_bin=1:6,
    response_bin=c(2,2,1,1,1,2),
    bin.x=      c(1,1,2,2,4,4),
    bin.y=1
  )
  
  model = create.model(df, nresponse.bins=2)
  
  likelihood = model$likelihood %>% arrange(bin.x, bin.y, bin.response)
  expect_equal(nrow(likelihood), 6)
  expect_equal(likelihood$prob.r, c(0.25, 0.75, 0.75, 0.25, 0.5, 0.5))
  
  expect_equal(model$prior$prob, rep(1/3,3))
})

test_that('create.model2 success for one cell', {
  stimulus = 1:4
  response = matrix(c(2,2,1,1,1,1,2,2), nrow = 2, byrow=TRUE) 
  rownames(response) = c('10', '12') # cell names
  
  model = create.model2(response, stimulus, 2, nstim.bins = 20)
  
  likelihood = model$likelihood 
  expect_equal(nrow(likelihood), 20)
  expect_equal(likelihood[1:4,1,1], c(1/3, 1/3, 2/3, 2/3))
  expect_equal(likelihood[1:4,1,2], c(2/3, 2/3, 1/3, 1/3))
  expect_equal(colnames(likelihood), rownames(response))
  
  expect_equal(model$prior, c(rep(0.25,4), rep(0, 16)))
})

test_that('bayes.max.s2 success with equal prior and same cells', {
  likelihoodM = array(0.1, dim=c(4,2,2))
  colnames(likelihoodM) = c('1', '2')
  likelihoodM[3,'1',1] = 0.9
  likelihoodM[3,'2',2] = 0.9
  
  model.bayes = list(prior=rep(0.25, 4),
                     likelihood=likelihoodM)
  pv = c(`1`=1, `2`=2)
  #res = bayes.max.s2(model.bayes, pv)
  res = bayesmax(model.bayes$prior, model.bayes$likelihood, pv)
  expect_equal(3, res$s)
})

test_that('bayes.max.s2 success with equal prior and different cells', {
  likelihoodM = array(0.1, dim=c(4,2,2))
  colnames(likelihoodM) = c('1', '2')
  likelihoodM[1,'1',1] = 0.9
  likelihoodM[4,'2',2] = 0.7
  likelihoodM[4,'2',1] = 0.3
  
  model.bayes = list(prior=rep(0.25, 4),
                     likelihood=likelihoodM)
  pv = c(`2`=1, `3`=2)
  #res = bayes.max.s2(model.bayes, pv)
  res = bayesmax(model.bayes$prior, model.bayes$likelihood, pv)
  expect_equal(4, res$s)
})

test_that('bayes.max.s2 success with unequal prior', {
  likelihoodM = array(0.25, dim=c(4,2,2))
  colnames(likelihoodM) = c('1', '2')
  likelihoodM[3,'1',1] = 0.7
  likelihoodM[3,'1',2] = 0.3
  
  model.bayes = list(prior=c(0.7, rep(0.1, 3)),
                     likelihood=likelihoodM)
  pv = c(`1`=1)
  res = bayesmax(model.bayes$prior, model.bayes$likelihood, pv)
  #res = bayes.max.s2(model.bayes, pv)
  expect_equal(1, res$s)
})

