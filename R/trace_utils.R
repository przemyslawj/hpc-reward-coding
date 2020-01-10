library(data.table)
library(purrr)

source('utils.R')

zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
}

sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))


read.data.trace = function(caimg_result_dir, filter_exp_title = NA) {
  print(paste('Processing dir: ', caimg_result_dir))
  cellmapping_file = file.path(caimg_result_dir, 'cell_mapping.csv')
  cellmapping.df = read.csv(cellmapping_file)
  traces_file = file.path(caimg_result_dir, 'traces_and_positions.csv')
  data = fread(traces_file)
  
  if (!is.na(filter_exp_title)) {
    data = data[exp_title == filter_exp_title,]  
  }
  
  data.traces = melt.traces(data)
  data.traces = left_join(data.traces, cellmapping.df, by=c('cell'='cell_no'))
  dir_parts = str_split(caimg_result_dir, '/')
  animal = dir_parts[[1]][length(dir_parts[[1]]) - 2]
  data.traces$animal = animal
  data.traces = data.table(data.traces)
  data.traces = data.traces[smooth_trans_x >= 0 & smooth_trans_y >= 0, ]
  
  setnames(data.traces, c('smooth_trans_x', 'smooth_trans_y', 'smooth_heading_angle'), c('x', 'y', 'angle'))
  return(data.traces)
}

calc.event.vec = function(deconv_trace, deconv.threshold=0.1) {
  deconv.thr.val = deconv.threshold * max(deconv_trace)
  res = map_lgl(deconv_trace, ~ .x >= deconv.thr.val)
  res
}

detect.events = function(data.traces, deconv.threshold=0.1) {
  data.traces[, is.event := calc.event.vec(.SD$deconv_trace, deconv.threshold) , by=c('animal', 'date', 'cell_id')]
  return(data.traces)
}

timebin.traces = function(data.traces, timebin.dur.msec=200) {
  if (nrow(data.traces) == 0) {
    return(data.frame())
  }
  
  max.timestamp = max(data.traces$timestamp)
  
  data.traces[, abs_timestamp := (trial-1)*(max.timestamp+timebin.dur.msec) + timestamp]
  data.traces[, time_bin := floor(abs_timestamp/timebin.dur.msec) %>% as.integer]
  timebinned.traces = data.traces[, 
                                  lapply(.SD, mean), 
                                  by=.(animal, date, trial_id, trial, exp_title, cell_id, time_bin),
                                  .SDcols = !c('is.event', 'dist')]
  setorder(timebinned.traces, time_bin, cell_id)
  nevents.traces = data.traces[, .(nevents=sum(is.event)), 
                                 by=.(animal, date, trial_id, trial, exp_title, cell_id, time_bin)]
  setorder(nevents.traces, time_bin, cell_id)
  timebinned.traces$nevents = nevents.traces$nevents
  timebinned.traces
}

stimbin.traces = function(data.traces, stim.var, bin.width=100/xybins) {
  stim.var = enquo(stim.var)
  
  bin.var.name = paste0('bin.', quo_name(stim.var))
  data.traces %>%
    dplyr::mutate(!!bin.var.name := as.integer(floor(!!stim.var / bin.width)))
}


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


bin.responses = function(df, quantile.fractions, binned.var='trace') {
  df = data.table(df)
  binned.df = df[, response_bin := get.response.bin(.SD[[binned.var]], quantile.fractions), 
                 by=c('animal', 'date', 'cell_id')]
  binned.df
}

nevents.bin.responses = function(df, nevents.thr=0.5) {
  df = data.table(df)
  binned.df = df[, response_bin := ifelse(nevents >= nevents.thr, 2, 1)]
  binned.df
}

bin.time.space = function(data.traces, bin.quantile.fractions=NULL, binned.var='trace', timebin.dur.msec=200, bin.width=100/xybins) {
  timebinned.traces = timebin.traces(data.traces[x >= 0 & y >= 0, ],
                                     timebin.dur.msec = timebin.dur.msec) %>%
    stimbin.traces(x, bin.width) %>%
    stimbin.traces(y, bin.width) %>%
    data.table()
  if (!is.null(bin.quantile.fractions)) {
    binned.traces = bin.responses(timebinned.traces, bin.quantile.fractions, binned.var=binned.var)
  }
  binned.traces[, bin.xy := to_1dim(bin.x, bin.y)]
  
  binned.traces
}

zscore.traces = function(data) {
  ddply(data, .(animal, date, trial_id, cell), plyr::mutate,
        timestamp = timestamp,
        smooth_trans_x = smooth_trans_x,
        smooth_trans_y = smooth_trans_y,
        dist = dist,
        velocity = velocity,
        inside_roi = inside_roi,
        nevents = nevents,
        trace = trace,
        ztrace=zscore(trace))
}

# Fast implementation of melting with data.tables
melt.traces = function(data) {
  data[, grep("^events_", colnames(data)):=NULL]
  trace.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'trace_')]
  #events.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'events_')]
  deconv.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'deconvTrace_')]
  melted.df = melt(data, 
       measure = list(trace.measure.vars, 
                      #events.measure.vars, 
                      deconv.measure.vars), 
       value.name = c('trace', 
                      #'nevents', 
                      'deconv_trace'))
  
  melted.df = melted.df[, ('cell') := variable %>% trimws %>% as.integer][order(date,trial_id,cell,timestamp), !'variable']
  return(melted.df)
}

gather.traces = function(data) {
  data.traces = data %>%
    select(-starts_with('events_'), -starts_with('deconvTrace_')) %>%
    gather('cell', 'trace', starts_with('trace_')) %>%
    mutate(cell = str_replace(cell, 'trace_[.]?[.]?','')) %>%
    arrange(date, trial_id, cell, timestamp)
  
  data.deconv_traces = data %>%
    select('date', 'trial_id', 'timestamp', starts_with('deconvTrace_')) %>%
    gather('cell', 'deconv_trace', starts_with('deconvTrace_')) %>%
    mutate(cell = str_replace(cell, 'deconvTrace_[.]?[.]?','')) %>%
    arrange(date, trial_id, cell, timestamp)
  
  #data.events = data %>%
  #  select('date', 'trial_id', 'timestamp', starts_with('events_')) %>%
  #  gather('cell', 'nevents', starts_with('events_')) %>%
  #  mutate(cell = str_replace(cell, 'events_[.]?[.]?','')) %>%
  #  arrange(date, trial_id, cell, timestamp)
  
  data = cbind(data.traces, 
               data.deconv_traces$deconv_trace)
               #data.events$nevents)
  #colnames(data)[(ncol(data)-1):ncol(data)] = c('deconv_trace', 'nevents')
  colnames(data)[(ncol(data))] = 'deconv_trace'
  data$cell = as.integer(data$cell)
  
  return(data)
}


# Filter based on velocity calculated for the time window. Filtering is done on epochs
# Slow code, use isRunning for better performance.
filter.running = function(df, min.run.velocity=2, mean.run.velocity=3, window.dur.ms=500) {
  df = data.table(df)
  setorder(df, trial_id, exp_title, cell_id, timestamp)
  is.running=rep(FALSE, nrow(df))
  i = 1
  while (i < nrow(df)) {
    j = i + 1
    while(j < nrow(df) && 
          df$velocity[i] >= min.run.velocity && 
          df$trial_id[i] == df$trial_id[j] &&
          df$velocity[j] >= min.run.velocity) {
      j = j + 1
    }
    j = j - 1
    
    dist = norm2(df$smooth_trans_x[j] - df$smooth_trans_x[i],
                 df$smooth_trans_y[j] - df$smooth_trans_y[i])
    dur.ms = max(df$timestamp[j] - df$timestamp[i], 1)
    vel = dist / dur.ms * 1000
    is.running[i:j] = (vel >= mean.run.velocity) && (dur.ms >= window.dur.ms)
    i = j + 1
  }
  
  df[which(is.running),]
}

# Map traces to principal components
pca.binned.traces = function(binned.traces) {
  response.matrix = reshape2::acast(binned.traces, time_bin ~ cell_id, value.var='response_bin') 
  pca.res = prcomp(response.matrix, center = TRUE, scale. = TRUE)
  cumvar = cumsum(pca.res$sdev*pca.res$sdev)
  cumvar.portion = cumvar / cumvar[length(cumvar)]
  npc = which(cumvar.portion >= 0.9) %>% first
  response.df = as.data.frame(pca.res$x[,1:npc])
  response.df$time_bin = as.integer(rownames(response.matrix))
  pc.measure.vars = colnames(response.df)[stringr::str_starts(colnames(response.df), 'PC')]
  melt.pc.df = melt(response.df, measure=pc.measure.vars, value.name='mean.trace') %>%
    data.table()
  melt.pc.df[,cell_id:= as.integer(stringr::str_replace(variable, 'PC', ''))]
  metadata.df = dplyr::select(binned.traces, animal, date, trial_id, trial, time_bin, bin.xy, bin.x, bin.y) %>%
    distinct() %>% data.table()
  melt.pc.df = melt.pc.df[metadata.df, on='time_bin']
  
  quantile.fractions = c(0.2, 0.4, 0.6, 0.8, 1.0)
  binned.pc.df = bin.responses(melt.pc.df, quantile.fractions)
  setkey(binned.pc.df, time_bin, cell_id)
  setorder(binned.pc.df, time_bin, cell_id)
  
  return(binned.pc.df)
}
