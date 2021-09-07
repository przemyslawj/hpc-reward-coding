library(testthat)

###################################
context('Tests for find peaks')
###################################

test_that('Returns a single point peak', {
  field = matrix(0, nrow=20, ncol=20)
  field[5, 5] = 1
  peaks.rowcol = fieldPeaks(field, minpeaksize=1, min.peakheight=0.5, sigma=0.0)
  expect_equal(dim(peaks.rowcol), c(1,2))
  expect_equal(c(5,5), as.vector(peaks.rowcol[1,]))
})

test_that('Returns no peaks when peak size too small', {
  field = matrix(0, nrow=20, ncol=20)
  field[5, 5] = 1
  peaks.rowcol = fieldPeaks(field, minpeaksize=2, min.peakheight=0.5, sigma=0.0)
  expect_equal(dim(peaks.rowcol), c(0,2))
})

test_that('Returns two peaks', {
  field = matrix(0, nrow=20, ncol=20)
  # Peak 1
  field[4, 1] = 0.5
  field[5, 1] = 1
  field[6, 1] = 0.5

  # Peak 2
  field[7, 4] = 0.5
  field[7, 5] = 1
  field[7, 6] = 0.5
  peaks.rowcol = fieldPeaks(field, minpeaksize=2, min.peakheight=0.4, sigma=0.0)
  expect_equal(dim(peaks.rowcol), c(2,2))
  expect_equal(c(5,1), as.vector(peaks.rowcol[1,]))
  expect_equal(c(7,5), as.vector(peaks.rowcol[2,]))
})

test_that('Returns peak COM', {
  field = matrix(0, nrow=20, ncol=20)
  # Peak 1
  field[5, 4] = 1
  field[5, 5] = 0.5
  field[5, 6] = 0.5

  # Peak 2
  field[1, 10] = 0.5
  field[2, 10] = 1
  field[3, 10] = 0.25

  peaks.rowcol = fieldPeaks(field, minpeaksize=2, min.peakheight=0.5, sigma=0.0)
  expect_equal(dim(peaks.rowcol), c(2,2))
  expect_equal(c(5/3, 10), as.vector(peaks.rowcol[1,]))
  expect_equal(c(5, 4.75), as.vector(peaks.rowcol[2,]))
})

###########################################
context('Tests for attraction strength')
###########################################

test_that('Attraction strength correct', {
  res = calc.attraction.strength(sqrt(2), 45, sqrt(2), 45)
  expect_equal(res, 0)
  
  res = calc.attraction.strength(1, 45, 0, 45)
  expect_equal(res, 1)
  
  res = calc.attraction.strength(1, 45, 1, 0)
   expect_equal(res, 0)
  
  res = calc.attraction.strength(1, 45, 2, 45)
  expect_equal(res, 1/3)
  
  # Right Triangle with dist1=5, dist2=4, x=3
  res = calc.attraction.strength(5, acos(4/5) * 180/pi, 4, 0)
  expect_equal(res, 3/5 * (5-4) / (5+4))
  
  # Test parallel processing
  res = calc.attraction.strength(c(1, 1, 1, 1), 
                                 c(45, 45, 45, 45), 
                                 c(1, 0, 1, 2), 
                                 c(45, 45, 0, 45))
  expect_equal(res, c(0, 1, 0, 1/3))
})