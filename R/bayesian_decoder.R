library(data.table)
library(dplyr)
library(purrr)
library(Rcpp)

Rcpp::sourceCpp('bayesian_decoder.cpp')

xybins = 20

# Calculates prior and likelihood used by Bayesian decoder
# response: matrix of trace values cell x time
# stimulus: vector of stimulus values length equal time samples
# nresponse.bins: count of bins for trace values
#
# Returns a list with fields
# prior: prior probability of each stimulus, length equal nstim.bins
# likelihood: 3D array: nstim.bins x ncells x nresponse.bins with a probability
# of a cell giving particular response for a given stimulus bin
create.discrete.bayes = function(training.df, nstim.bins=20*20, value.var=response_bin, stim.var=bin.xy) {
  stim.var = enquo(stim.var)
  value.var.name = quo_name(enquo(value.var))
  
  # to matrix representation
  response = reshape2::acast(training.df, cell_id ~ time_bin, value.var=value.var.name)
  setorder(training.df, time_bin)
  stimulus = training.df %>%
    dplyr::select(time_bin, !! stim.var) %>%
    dplyr::arrange(time_bin) %>%
    dplyr::distinct()
  stimulus = stimulus[[2]]
  
  occupancy.hist = hist(stimulus, breaks=0:nstim.bins, plot=FALSE)
  Moccupancy = occupancy.hist$counts
  visited = which(Moccupancy > 0)
  
  nresponse.bins = max(training.df[[value.var.name]])
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

create.mfr.bayes = function(training.df, nstim.bins=20*20, value.var=nevents, stim.var=bin.xy) {
  stim.var = enquo(stim.var)
  stim.var.name = quo_name(stim.var)
  value.var = enquo(value.var)
  value.var.name = quo_name(value.var)
  training.df = dplyr::rename(training.df, s = !!stim.var)
  
  sample.cell.df = training.df[, head(.SD, 1), by = time_bin]
  occupancy.hist = hist(sample.cell.df$s, breaks=0:nstim.bins, plot=FALSE)
  cell.ids = unique(training.df$cell_id)
  full.xy.df = data.frame(s = rep(1:nstim.bins, length(cell.ids)), 
                          nevents=NA, 
                          cell_id = rep(cell.ids, each=nstim.bins))
  # Add one event to each visited bin to avoid 0 probability of the bin
  one.event.df = select(training.df, cell_id, s) %>% 
    distinct() %>%
    mutate(nevents=1)
  
  model.df = bind_rows(select(training.df, cell_id, s, !!value.var), full.xy.df, one.event.df)
  mfr.matrix = reshape2::acast(model.df, cell_id ~ s, value.var=value.var.name, fun.aggregate=function(x) {mean(x, na.rm=TRUE)})
  sd.matrix = reshape2::acast(model.df, cell_id ~ s, value.var=value.var.name, fun.aggregate=function(x) {sd(x, na.rm=TRUE)})
  M = abind::abind(mfr.matrix, sd.matrix, along=3)
  return(list(occupancy.counts=occupancy.hist$counts,
              prior=occupancy.hist$density,
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

bin.distance.error = function(expected.xy, actual.xy) {
  actual.s.bin = from_1dim(actual.xy)
  expected.s.bin = from_1dim(expected.xy)
  dist = norm2(actual.s.bin$x - expected.s.bin$x, actual.s.bin$y - expected.s.bin$y)
  return(dist)
}

eval.testdata2 = function(test.df, 
                          model.bayes, 
                          predict.fun=bayesmax, 
                          error.fun=bin.distance.error,
                          stim.var=bin.xy,
                          value.var=response_bin, 
                          nstim.bins=xybins^2) {
  stim.var = enquo(stim.var)
  stim.var.name = quo_name(stim.var)
  value.var.name = quo_name(enquo(value.var))
  actual.stim.var.name = paste0('actual.', stim.var.name)
  
  test.response = reshape2::acast(test.df, cell_id ~ time_bin, value.var=value.var.name)
  results.df = test.df %>%
    mutate(prior = model.bayes$prior[test.df[[stim.var.name]]]) %>%
    select(time_bin, !!stim.var, prior) %>%
    arrange(time_bin) %>%
    distinct()
  
  results.df$error = 0
  results.df[[actual.stim.var.name]] = 0
  density.sum = matrix(0, nrow=nstim.bins,ncol=nstim.bins)
  stimulus.counts = rep(0, nstim.bins)
  for (i in 1:ncol(test.response)) {
    pv = test.response[, i]
    test.timestamp = as.integer(colnames(test.response)[i])
    expected.s = results.df[[stim.var.name]][i]
    actual.s.res = predict.fun(model.bayes$prior, model.bayes$likelihood, pv)
    results.df[[actual.stim.var.name]][i] = actual.s.res$s
    
    density = actual.s.res$density
    density[is.na(density)] = 0
    density.prob = density / sum(density)

    density.sum[expected.s, ] = density.sum[expected.s, ] + density
    stimulus.counts[expected.s] = stimulus.counts[expected.s] + 1
    
    results.df$error[i] = error.fun(expected.s, actual.s.res$s) 
  }
  

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

filter.sampled = function(training.df, test.df, min.samples=20) {
  bin.samples = training.df[, .(nsamples=length(unique(time_bin))), by=bin.xy]
  train.filtered = bin.samples[training.df, on='bin.xy'][nsamples >= min.samples,]
  test.filtered = bin.samples[test.df, on='bin.xy'][nsamples >= min.samples,]
  
  list(train=training.df,
       test=test.filtered)
}

eval.decoder = function(binned.traces, 
                        train.fun=create.discrete.bayes,
                        predict.fun=bayesmax,
                        stim.var=bin.xy,
                        value.var=response_bin,
                        nstim.bins=xybins * xybins,
                        training.split.fraction=0.8, 
                        cv=TRUE,
                        error.fun=bin.distance.error) {
  stim.var = enquo(stim.var)
  value.var = enquo(value.var)
  
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

    model.bayes = train.fun(training.df, 
                            nstim.bins=nstim.bins, 
                            value.var=!!value.var,
                            stim.var=!!stim.var)
    
    #smooth.model.bayes = smooth.likelihoods(model.bayes)
    if (nrow(test.df) == 0) {
      warning('0 rows for testing after filtering with min samples')
    } else {
      system.time(eval.res <- eval.testdata2(test.df, 
                                             model.bayes, 
                                             predict.fun=predict.fun,
                                             stim.var=!!stim.var,
                                             value.var=!!value.var,
                                             nstim.bins=nstim.bins,
                                             error.fun=error.fun))
      partial.df = eval.res$df
      
      random.classifier.res <- eval.testdata2(test.df, 
                                              model.bayes, 
                                              predict.fun=random.prior.classifier,
                                              stim.var=!!stim.var,
                                              value.var=!!value.var,
                                              nstim.bins=nstim.bins,
                                              error.fun=error.fun)
      partial.df$random_error = random.classifier.res$df$error
      partial.df$cv=i
      result.df = bind_rows(result.df, partial.df)
    }
  }
  
  return(result.df)
}
