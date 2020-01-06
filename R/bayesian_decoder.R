library(data.table)
library(dplyr)
library(purrr)
library(Rcpp)

Rcpp::sourceCpp('bayesian_decoder.cpp')

xybins = 20
run.threshold  = 2

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


# Calculates prior and likelihood used by Bayesian decoder
# response: matrix of trace values cell x time
# stimulus: vector of stimulus values length equal time samples
# nresponse.bins: count of bins for trace values
#
# Returns a list with fields
# prior: prior probability of each stimulus, length equal nstim.bins
# likelihood: 3D array: nstim.bins x ncells x nresponse.bins with a probability
# of a cell giving particular response for a given stimulus bin
create.model2 = function(response, stimulus, nresponse.bins, nstim.bins=20) {
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

# # Returns DF with probability of a binned response for a given location
# create.model = function(training.df, nresponse.bins) {
#   training.df = data.table(training.df)
#   occupancy.df = training.df %>%
#     group_by(time_bin) %>%
#     slice(1) # One row per cell, take one
#   
#   Moccupancy = reshape2::acast(occupancy.df, bin.x ~ bin.y, fun.aggregate=length, value.var='time_bin')
#   occupancy.probs.df = melt(Moccupancy)
#   total.occupancy = sum(occupancy.probs.df$value)
#   occupancy.probs.df = mutate(occupancy.probs.df, 
#                               prob=value/total.occupancy,
#                               bin.xy=to_1dim(Var1, Var2)) %>%
#     dplyr::rename(bin.x=Var1, bin.y=Var2) %>%
#     select(-value) %>%
#     data.table()
#   setkey(occupancy.probs.df, bin.xy)
#   
#   visited = which(Moccupancy > 0)
#   # Make unseen bins probability > 0 by adding 1 to each response_bin count
#   Moccupancy[visited] = Moccupancy[visited] + nresponse.bins
#   M.add = matrix(0, nrow=dim(Moccupancy)[1], ncol=dim(Moccupancy)[2])
#   M.add[visited] = M.add[visited] + 1
#   rownames(M.add) = rownames(Moccupancy)
#   colnames(M.add) = colnames(Moccupancy)
#   
#   all.field.df = data.frame()
#   for (cell_name in training.df$cell_id %>% unique) {
#     cell.df = training.df[cell_id == cell_name, .(bin.x, bin.y, response_bin)]
# 
#     field.df = data.frame()
#     for (i in 1:nresponse.bins) {
#       bin_name = paste('bin', i, sep='_')
#       cell.df[, bin_name] = i == cell.df$response_bin
#       M = reshape2::acast(cell.df, bin.x ~ bin.y, fun.aggregate=sum, value.var=bin_name)
#       Moccupancy.sub = Moccupancy[rownames(M),colnames(M)] %>% as.matrix
#       M.add.sub = M.add[rownames(M),colnames(M)] %>% as.matrix
#       
#       M = M + M.add.sub
#       M = M / Moccupancy.sub
#       bin.field.df = melt(M) %>% 
#         dplyr::filter(!is.na(value)) %>%
#         dplyr::mutate(bin.xy = to_1dim(Var1, Var2)) %>%
#         dplyr::rename(prob.r=value) %>% 
#         dplyr::filter(is.finite(prob.r), !is.na(prob.r)) %>% # Skip if occupancy was 0 (x/0 = Inf)
#         dplyr::rename(bin.x=Var1, bin.y=Var2) %>%
#         dplyr::select(bin.xy, bin.x, bin.y, prob.r)
#       bin.field.df$bin.response = i
#       
#       field.df = bind_rows(field.df, bin.field.df)
#     }
#     
#     field.df$cell_id = cell_name
#     field.df$animal = training.df$animal[1]
#     
#     all.field.df = bind_rows(all.field.df, field.df)
#   }
#   
#   all.field.df = data.table(all.field.df)
#   setkey(all.field.df, bin.xy, cell_id, bin.response)
#   
#   return(list(prior=occupancy.probs.df, 
#               likelihood=all.field.df))
# }
# 
# 
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

eval.testdata2 = function(test.df, model.bayes, classifier.fun) {
  test.response = reshape2::acast(test.df, cell_id ~ time_bin, value.var='response_bin')
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

filter.sampled = function(training.df, test.df, min.samples=6) {
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
  
  binned.traces = binned.traces %>%
    mutate(bin.xy = to_1dim(bin.x, bin.y))
  binned.traces = data.table(binned.traces)
  setorder(binned.traces, time_bin, cell_id)
  setkey(binned.traces, time_bin, cell_id)
  
  test.index.start = nrow(binned.traces)
  for (i in 1:ncv) {
    test.index.end = find.first.timestamp(binned.traces$time_bin, test.index.start, 1)
    test.index.start = find.first.timestamp(binned.traces$time_bin, max(1, test.index.start - nrows.testing), -1)
    test.indecies = test.index.start:test.index.end
    train.indecies = setdiff(1:nrow(binned.traces), test.indecies)
    #training.df = binned.traces[train.indecies,]
    #test.df = binned.traces[test.indecies,]
    filtered.dfs = filter.sampled(binned.traces[train.indecies,], binned.traces[test.indecies,])
    training.df = filtered.dfs$train
    test.df = filtered.dfs$test

    # to matrix representation
    training.response = reshape2::acast(training.df, cell_id ~ time_bin, value.var='response_bin')
    stimulus = training.df %>%
      select(time_bin, bin.xy) %>%
      arrange(time_bin) %>%
      distinct()
    model.bayes = create.model2(training.response, stimulus$bin.xy, nresponse.bins, nstim.bins=xybins*xybins)
    #smooth.model.bayes = smooth.likelihoods(model.bayes)
    if (nrow(test.df) == 0) {
      warning('0 rows for testing after filtering with min samples')
    } else {
      system.time(eval.res <- eval.testdata2(test.df, model.bayes, bayesmax))
      partial.df = eval.res$df
      
      random.classifier.res <- eval.testdata2(test.df, model.bayes, random.prior.classifier)
      partial.df$random_error = random.classifier.res$df$error
      partial.df$cv=i
      result.df = bind_rows(result.df, partial.df)
    }
  }
  
  return(result.df)
}
