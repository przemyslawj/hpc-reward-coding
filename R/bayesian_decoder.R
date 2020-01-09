library(data.table)
library(dplyr)
library(purrr)
library(Rcpp)

Rcpp::sourceCpp('bayesian_decoder.cpp')

xybins = 20

to_1dim = function(x, y) {
  # Vals from 1 to xybins * xybins
  xybins * x + y + 1
}

from_1dim = function(z) {
  z = z - 1
  list(x=floor(z/xybins), y=z%%xybins)
}

get.response.bin = function(vals, quantile.fractions) {
  trace.quantiles = quantile(vals, quantile.fractions) + 0.001
  map_int(vals, ~ dplyr::first(which(.x <= trace.quantiles)))
}


bin.responses = function(df, quantile.fractions) {
  binned.df = df[, response_bin := get.response.bin(.SD$mean.trace, quantile.fractions), by=c('animal', 'date', 'cell_id')]
  binned.df
}

nevents.bin.responses = function(df, nevents.thr=0.5) {
  binned.df = df[, response_bin := ifelse(nevents >= nevents.thr, 2, 1)]
  binned.df
}

# Calculates prior and likelihood used by Bayesian decoder
# response: matrix of trace values cell x time
# stimulus: vector of stimulus values length equal time samples
# nresponse.bins: count of bins for trace values
#
# Returns a list with fields
# prior: prior probability of each stimulus, length equal nstim.bins
# likelihood: 3D array: nstim.bins x ncells x nresponse.bins with a probability
# of a cell giving particular response for a given stimulus bin
create.model2 = function(training.df, nresponse.bins, nstim.bins=20*20, value.var='response_bin') {
  # to matrix representation
  response = reshape2::acast(training.df, cell_id ~ time_bin, value.var=value.var)
  setorder(training.df, time_bin)
  stimulus = training.df %>%
    dplyr::select(time_bin, bin.xy) %>%
    dplyr::arrange(time_bin) %>%
    dplyr::distinct()
  stimulus = stimulus$bin.xy
  
  occupancy.hist = hist(stimulus, breaks=0:nstim.bins, plot=FALSE)
  Moccupancy = occupancy.hist$counts
  visited = which(Moccupancy > 0)
  
  Moccupancy[visited] = Moccupancy[visited] + nresponse.bins
  ncells = nrow(response)
  
  M = array(0, dim=c(nstim.bins, ncells, nresponse.bins))
  # Make unseen bins probability > 0 by adding 1 to each response_bin count
  M[visited,,] = 1
  for (s_i in 1:length(stimulus)) {
    for (cell_i in 1:ncells) {
      response_bin = response[cell_i, s_i];
      M[stimulus[s_i], cell_i, response_bin] = M[stimulus[s_i], cell_i, response_bin]  + 1
    }
  }
  
  Mlikelihood = M / Moccupancy
  colnames(Mlikelihood) = rownames(response)
  return(list(prior=occupancy.hist$density,
              likelihood=Mlikelihood))
  
}

create.mfr.model = function(training.df, nstim.bins=20*20) {
  occupancy.hist = hist(training.df$bin.xy, breaks=0:nstim.bins, plot=FALSE)
  cell.ids = unique(training.df$cell_id)
  full.xy.df = data.frame(bin.xy = rep(1:nstim.bins, length(cell.ids)), 
                          nevents=NA, 
                          cell_id = rep(cell.ids, each=nstim.bins))
  # Add one event to each visited bin to avoid 0 probability of the bin
  one.event.df = select(training.df, cell_id, bin.xy) %>% 
    distinct() %>%
    mutate(nevents=1)
  
  model.df = bind_rows(select(training.df, cell_id, bin.xy, nevents), full.xy.df, one.event.df)
  mfr.matrix = reshape2::acast(model.df, cell_id ~ bin.xy, value.var='nevents', fun.aggregate=function(x) {mean(x, na.rm=TRUE)})
  sd.matrix = reshape2::acast(model.df, cell_id ~ bin.xy, value.var='nevents', fun.aggregate=function(x) {sd(x, na.rm=TRUE)})
  M = abind::abind(mfr.matrix, sd.matrix, along=3)
  return(list(prior=occupancy.hist$density,
              likelihood=M))
}


smooth.likelihoods = function(model.bayes, sigma=0.8) {
  nstim = dim(model.bayes$likelihood)[1]
  ncells = dim(model.bayes$likelihood)[2]
  nresponses = dim(model.bayes$likelihood)[3]
  
  likelihood.smooth = array(0, dim=dim(model.bayes$likelihood))
  
  for (cell_i in 1:ncells) {
    for (r in 1:nresponses) {
      M = matrix(model.bayes$likelihood[, cell_i, r], ncol=xybins, nrow=xybins, byrow = TRUE)
      M1 = gauss2dsmooth(M, lambda=sigma, nx=8, ny=8)
      M1[which(is.na(M))] = NA
          
      likelihood.smooth[, cell_i, r] = as.vector(t(M1))
    }
    # normalize to prob values
    cell.norm = rowSums(likelihood.smooth[, cell_i, ])
    likelihood.smooth[, cell_i, ] = likelihood.smooth[, cell_i,] / cell.norm
  }
  dimnames(likelihood.smooth) = dimnames(model.bayes$likelihood)
  model.bayes$likelihood = likelihood.smooth
  
  return(model.bayes)
}

# # Bayes rule (will assume P(s) uniform):
# # max_s [P(s|r)] ~ max_s P(r|s) * P(s)
# #
# # Assuming independent activity of the cells:
# # P(r|s) = II_i P(r_i|s)
# # pv - data frame with activation
# bayes.max.s = function(model.bayes, pv) {
#   model.df = model.bayes$likelihood
#   
#   if (length(pv) == 0) {
#     print('Error, pv empty')
#   }
#   #setkey(pv, cell_id, bin.response)
#   max.s = -1
#   max.s.prob = -1.0
#   for (s in model.df[,unique(bin.xy)]) {
#     #cells.probs = model.df[bin.xy==s,]
#     #if (nrow(cells.probs) > 0) {
#       # if a cell not present on the day, it will be excluded from probability calculation
#     #probs = pv[model.df[bin.xy==s,], prob.r, on=c('cell_id', 'bin.response')]
#     probs = model.df[bin.xy==s,][pv, prob.r, on=c('cell_id', 'bin.response')]
#     
#     if (length(probs) > 0) {
#       # P(r|s) = II_i P(r_i|s)
#       s.prob = cumprod(probs) %>% last
#       s.prob = s.prob * model.bayes$prior[bin.xy==s, prob]
#       
#       if (s.prob > max.s.prob) {
#         max.s.prob = s.prob
#         max.s = s
#       }
#     }
#   }
#   return(list(s=max.s, prob=max.s.prob))
# }
# 
# 
# eval.testdata = function(test.df, model.bayes) {
#   test.timestamps = test.df[,unique(time_bin)]
#   ntimestamps = length(test.timestamps)
#   error.norms = rep(0, ntimestamps)
#   expected.s = data.frame(bin.x=error.norms, bin.y=error.norms)
#   actual.s = data.frame(x=error.norms, y=error.norms)
# 
#   for (i in 1:length(test.timestamps)) {
#     test.timestamp = test.timestamps[i]
#     pv = test.df[time_bin==test.timestamp, .(cell_id, bin.response=response_bin)]
#     expected.s.row = test.df[time_bin==test.timestamp, .(bin.x, bin.y)] %>% first
#     expected.s[i,] = expected.s.row
#     actual.s.res = bayes.max.s(model.bayes, pv)
#     actual.s.bin = actual.s.res$s %>% from_1dim()
#     actual.s[i, ] = actual.s.bin
#     
#     error.norms[i] = norm2(actual.s.bin$x - expected.s.row$bin.x, actual.s.bin$y - expected.s.row$bin.y)
#   }
#   
#   results.df = cbind(expected.s, actual.s) %>%
#     dplyr::rename(actual.x=x, actual.y=y)
#   results.df$error = error.norms
#   return(results.df)
# }

eval.testdata2 = function(test.df, model.bayes, classifier.fun, value.var='response_bin') {
  test.response = reshape2::acast(test.df, cell_id ~ time_bin, value.var=value.var)
  expected.s = test.df %>%
    mutate(prior = model.bayes$prior[bin.xy]) %>%
    select(time_bin, bin.x, bin.y, bin.xy, prior) %>%
    arrange(time_bin) %>%
    distinct()
  
  error.norms = rep(0, ncol(test.response))
  actual.s = data.frame(x=error.norms, y=error.norms)
  density.sum = matrix(0, nrow=xybins^2,ncol=xybins^2)
  stimulus.counts = rep(0, xybins^2)
  for (i in 1:ncol(test.response)) {
    pv = test.response[, i]
    test.timestamp = as.integer(colnames(test.response)[i])
    expected.s.row = expected.s[i,]
    actual.s.res = classifier.fun(model.bayes$prior, model.bayes$likelihood, pv)
    actual.s.bin = actual.s.res$s %>% from_1dim()
    actual.s[i, ] = actual.s.bin
    
    density = actual.s.res$density
    density[is.na(density)] = 0
    density.prob = density / sum(density)
    density.sum[expected.s.row$bin.xy, ] = density.sum[expected.s.row$bin.xy, ] + density
    stimulus.counts[expected.s.row$bin.xy] = stimulus.counts[expected.s.row$bin.xy] + 1
    
    error.norms[i] = norm2(actual.s.bin$x - expected.s.row$bin.x, actual.s.bin$y - expected.s.row$bin.y)
  }
  
  results.df = cbind(expected.s, actual.s) %>%
    dplyr::rename(actual.x=x, actual.y=y)
  results.df$error = error.norms
  stimulus.counts[which(stimulus.counts == 0)] = 1
  return(list(df=results.df, density=density.sum / stimulus.counts))
}

random.prior.classifier = function(prior, likelihood, pv) {
  r = runif(1)
  prior[which(is.na(prior))] = 0
  s = which(cumsum(prior) > r) %>% first
  list(s=s, prob=prior[s], density=rep(1 / length(prior), length(prior)))
}


find.first.timestamp = function(timestamps, ind, inc=1) {
  while ((ind + inc > 0) &&
         (ind + inc <= length(timestamps)) &&
         (timestamps[ind] == timestamps[ind+inc])) {
    ind = ind + inc
  }
  
  return(ind)
}

filter.sampled = function(training.df, test.df, min.samples=10) {
  bin.samples = training.df[, .(nsamples=length(unique(time_bin))), by=bin.xy]
  train.filtered = bin.samples[training.df, on='bin.xy'][nsamples >= min.samples,]
  test.filtered = bin.samples[test.df, on='bin.xy'][nsamples >= min.samples,]
  
  list(train=training.df,
       test=test.filtered)
}

eval.decoder = function(binned.traces, nresponse.bins, training.split.fraction=0.8, cv=TRUE) {
  ncv = ifelse(cv, (1.0 / (1.0 - training.split.fraction)) %>% floor %>% as.integer, 1)
  
  nrows.training = floor(training.split.fraction * nrow(binned.traces))
  nrows.testing = nrow(binned.traces) - nrows.training
  result.df = data.frame()
  
  binned.traces = data.table(binned.traces)
  
  test.index.start = nrow(binned.traces)
  for (i in 1:ncv) {
    test.index.end = find.first.timestamp(binned.traces$time_bin, test.index.start, 1)
    test.index.start = find.first.timestamp(binned.traces$time_bin, max(1, test.index.start - nrows.testing), -1)
    test.indecies = test.index.start:test.index.end
    train.indecies = setdiff(1:nrow(binned.traces), test.indecies)
    filtered.dfs = filter.sampled(binned.traces[train.indecies,], binned.traces[test.indecies,])
    training.df = filtered.dfs$train
    test.df = filtered.dfs$test

    model.bayes = create.model2(training.df, nresponse.bins, nstim.bins=xybins*xybins, value.var='response_bin')
    #mfr.model = create.mfr.model(training.df, nstim.bins=xybins*xybins)
    #smooth.model.bayes = smooth.likelihoods(model.bayes)
    if (nrow(test.df) == 0) {
      warning('0 rows for testing after filtering with min samples')
    } else {
      system.time(eval.res <- eval.testdata2(test.df, model.bayes, bayesmax))
      #system.time(eval.res <- eval.testdata2(test.df, mfr.model, bayesmax_mfr, value.var='nevents'))
      partial.df = eval.res$df
      
      random.classifier.res <- eval.testdata2(test.df, model.bayes, random.prior.classifier)
      partial.df$random_error = random.classifier.res$df$error
      partial.df$cv=i
      result.df = bind_rows(result.df, partial.df)
    }
  }
  
  return(result.df)
}
