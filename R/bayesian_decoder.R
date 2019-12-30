library(data.table)
library(dplyr)
library(purrr)

#TODO: bayes.max.s is slow

xybins = 20
bin.width = 100/xybins
run.threshold  = 2

to_1dim = function(x, y) {
  xybins * x + y
}

from_1dim = function(z) {
  list(x=floor(z/xybins), y=z%%xybins)
}

timebin.traces = function(data.traces, timebin.dur.msec=100) {
  max.timestamp = max(data.traces$timestamp)
  timebinned.traces = data.traces %>%
    mutate(abs_timestamp = (trial-1)*(max.timestamp+timebin.dur.msec) + timestamp,
           time_bin = floor(abs_timestamp/timebin.dur.msec) %>% as.integer) %>%
    group_by(animal, date, trial_id, cell_id, time_bin) %>%
    dplyr::summarise(mean.trace = mean(trace),
                     mean.x = mean(smooth_trans_x),
                     mean.y = mean(smooth_trans_y)) %>%
    mutate(bin.x = floor(mean.x / bin.width),
           bin.y = floor(mean.y / bin.width))
  
  timebinned.traces = data.table(timebinned.traces) %>%
    setorder(time_bin, cell_id)
  timebinned.traces
}


get.response.bin = function(vals, animal_name, date_str, cell_name, cell.thresholds, quantile.fractions) {
  trace.quantiles.row = cell.thresholds[animal==animal_name[1] & date==date_str[1] & cell_id==cell_name[1],]
  trace.quantiles = trace.quantiles.row[1, 4: (3 + length(quantile.fractions))] %>% as.matrix
  map_int(vals, ~ dplyr::first(which(.x <= trace.quantiles)))
}


bin.responses = function(df, quantile.fractions) {
  cell.thresholds = df %>% 
    ddply(.(animal, date, cell_id), function(cell.df) {
      quantile(cell.df$mean.trace, quantile.fractions) + 0.001
    })
  cell.thresholds = data.table(cell.thresholds)
  setkey(cell.thresholds, animal, date, cell_id)
  
  binned.df = df %>%
    dplyr::group_by(animal, date, cell_id) %>%
    dplyr::mutate(response_bin = get.response.bin(mean.trace, animal, date, cell_id, cell.thresholds, quantile.fractions)) %>%
    dplyr::ungroup()
          
  binned.df
}

# Returns DF with probability of a binned response for a given location
create.model = function(training.df, nresponse.bins) {
  training.df = data.table(training.df)
  occupancy.df = training.df %>%
    group_by(time_bin) %>%
    slice(1) # One row per cell, take one
  
  Moccupancy = reshape2::acast(occupancy.df, bin.x ~ bin.y, fun.aggregate=length, value.var='time_bin')
  occupancy.probs.df = melt(Moccupancy)
  total.occupancy = sum(occupancy.probs.df$value)
  occupancy.probs.df = mutate(occupancy.probs.df, 
                              prob=value/total.occupancy,
                              bin.xy=to_1dim(Var1, Var2)) %>%
    dplyr::rename(bin.x=Var1, bin.y=Var2) %>%
    select(-value) %>%
    data.table()
  setkey(occupancy.probs.df, bin.xy)
  
  visited = which(Moccupancy > 0)
  # Make unseen bins probability > 0 by adding 1 to each response_bin count
  Moccupancy[visited] = Moccupancy[visited] + nresponse.bins
  M.add = matrix(0, nrow=dim(Moccupancy)[1], ncol=dim(Moccupancy)[2])
  M.add[visited] = M.add[visited] + 1
  rownames(M.add) = rownames(Moccupancy)
  colnames(M.add) = colnames(Moccupancy)
  
  all.field.df = data.frame()
  for (cell_name in training.df$cell_id %>% unique) {
    cell.df = training.df[cell_id == cell_name, .(bin.x, bin.y, response_bin)]

    field.df = data.frame()
    for (i in 1:nresponse.bins) {
      bin_name = paste('bin', i, sep='_')
      cell.df[, bin_name] = i == cell.df$response_bin
      M = reshape2::acast(cell.df, bin.x ~ bin.y, fun.aggregate=sum, value.var=bin_name)
      Moccupancy.sub = Moccupancy[rownames(M),colnames(M)] %>% as.matrix
      M.add.sub = M.add[rownames(M),colnames(M)] %>% as.matrix
      
      M = M + M.add.sub
      M = M / Moccupancy.sub
      bin.field.df = melt(M) %>% 
        dplyr::filter(!is.na(value)) %>%
        dplyr::mutate(bin.xy = to_1dim(Var1, Var2)) %>%
        dplyr::rename(prob.r=value) %>% 
        dplyr::filter(is.finite(prob.r), !is.na(prob.r)) %>% # Skip if occupancy was 0 (x/0 = Inf)
        dplyr::rename(bin.x=Var1, bin.y=Var2) %>%
        dplyr::select(bin.xy, bin.x, bin.y, prob.r)
      bin.field.df$bin.response = i
      
      field.df = bind_rows(field.df, bin.field.df)
    }
    
    field.df$cell_id = cell_name
    field.df$animal = training.df$animal[1]
    
    all.field.df = bind_rows(all.field.df, field.df)
  }
  
  all.field.df = data.table(all.field.df)
  setkey(all.field.df, bin.xy, cell_id, bin.response)
  
  return(list(prior=occupancy.probs.df, 
              likelihood=all.field.df))
}


# Bayes rule (will assume P(s) uniform):
# max_s [P(s|r)] ~ max_s P(r|s) * P(s)
#
# Assuming independent activity of the cells:
# P(r|s) = II_i P(r_i|s)
# pv - data frame with activation
bayes.max.s = function(model.bayes, pv) {
  model.df = model.bayes$likelihood
  
  if (length(pv) == 0) {
    print('Error, pv empty')
  }
  #setkey(pv, cell_id, bin.response)
  max.s = -1
  max.s.prob = -1.0
  for (s in model.df[,unique(bin.xy)]) {
    #cells.probs = model.df[bin.xy==s,]
    #if (nrow(cells.probs) > 0) {
      # if a cell not present on the day, it will be excluded from probability calculation
    #probs = pv[model.df[bin.xy==s,], prob.r, on=c('cell_id', 'bin.response')]
    probs = model.df[bin.xy==s,][pv, prob.r, on=c('cell_id', 'bin.response')]
    
    if (length(probs) > 0) {
      # P(r|s) = II_i P(r_i|s)
      s.prob = cumprod(probs) %>% last
      s.prob = s.prob * model.bayes$prior[bin.xy==s, prob]
      
      if (s.prob > max.s.prob) {
        max.s.prob = s.prob
        max.s = s
      }
    }
  }
  return(list(s=max.s, prob=max.s.prob))
}


eval.testdata = function(test.df, model.bayes) {
  test.timestamps = test.df[,unique(time_bin)]
  expected.s = data.frame()
  actual.s = data.frame()
  #setkey(test.df, time_bin)
  error.norms = rep(0, length(test.timestamps))
  expected.s = data.frame(bin.x=error.norms, bin.y=error.norms)
  actual.s = data.frame(x=error.norms, y=error.norms)
  for (i in 1:length(test.timestamps)) {
    test.timestamp = test.timestamps[i]
    pv = test.df[time_bin==test.timestamp, .(cell_id, bin.response=response_bin)]
    expected.s.row = test.df[time_bin==test.timestamp, .(bin.x, bin.y)] %>% first
    #expected.s = bind_rows(expected.s, expected.s.row)
    expected.s[i,] = expected.s.row
    actual.s.res = bayes.max.s(model.bayes, pv)
    actual.s.bin = actual.s.res$s %>% from_1dim()
    #actual.s = bind_rows(actual.s, actual.s.bin)
    actual.s[i, ] = actual.s.bin
    
    error.norms[i] = norm2(actual.s.bin$x - expected.s.row$bin.x, actual.s.bin$y - expected.s.row$bin.y)
  }
  
  results.df = cbind(expected.s, actual.s) %>%
    dplyr::rename(actual.x=x, actual.y=y)
  results.df$error = error.norms
  return(results.df)
}

eval.decoder = function(binned.traces, nresponse.bins, training.split.fraction=0.8, cv=TRUE) {
  ncv = ifelse(cv, (1.0 / (1.0 - training.split.fraction)) %>% floor %>% as.integer, 1)
  
  nrows.training = floor(training.split.fraction * nrow(binned.traces))
  nrows.testing = nrow(binned.traces) - nrows.training
  result.df = data.frame()
  
  binned.traces = data.table(binned.traces)
  setorder(binned.traces, time_bin, cell_id)
  setkey(binned.traces, time_bin, cell_id)
  
  train.index.start = 1
  for (i in 1:ncv) {
    train.indecies = train.index.start:(train.index.start + nrows.training)
    training.df = binned.traces[train.indecies,]
    test.indecies = setdiff(1:nrow(binned.traces), train.indecies)
    model.bayes = create.model(training.df, nresponse.bins)
    test.df = binned.traces[test.indecies, ]
    
    partial.df = eval.testdata(test.df, model.bayes)
    partial.df$cv=i
    result.df = bind_rows(result.df, partial.df)
    
    train.index.start = train.index.start + nrows.testing
  }
  
  return(result.df)
}
