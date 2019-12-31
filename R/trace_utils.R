library(data.table)

source('utils.R')

zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
}

sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))


read.data.trace = function(ca_img_result_dir, filter_exp_title = NA) {
  print(paste('Processing dir: ', ca_img_result_dir))
  cellmapping_file = file.path(ca_img_result_dir, 'cell_mapping.csv')
  cellmapping.df = read.csv(cellmapping_file)
  traces_file = file.path(ca_img_result_dir, 'traces_and_positions.csv')
  data = fread(traces_file)
  
  if (!is.na(filter_exp_title)) {
    data = data[exp_title == filter_exp_title,]  
  }
  
  data.traces = melt.traces(data)
  data.traces = left_join(data.traces, cellmapping.df, by=c('cell'='cell_no'))
  dir_parts = str_split(ca_img_result_dir, '/')
  animal = dir_parts[[1]][length(dir_parts[[1]]) - 2]
  data.traces$animal = animal
  data.traces = data.table(data.traces)
  data.traces = data.traces[smooth_trans_x >= 0 & smooth_trans_y >= 0, ]
  
  return(data.traces)
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
  trace.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'trace_')]
  events.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'events_')]
  deconv.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'deconvTrace_')]
  melted.df = melt(data, 
       measure = list(trace.measure.vars, events.measure.vars, deconv.measure.vars), 
       value.name = c('trace', 'nevents', 'deconv_trace'))
  
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
  
  data.events = data %>%
    select('date', 'trial_id', 'timestamp', starts_with('events_')) %>%
    gather('cell', 'nevents', starts_with('events_')) %>%
    mutate(cell = str_replace(cell, 'events_[.]?[.]?','')) %>%
    arrange(date, trial_id, cell, timestamp)
  
  data = cbind(data.traces, data.deconv_traces$deconv_trace, data.events$nevents)
  colnames(data)[(ncol(data)-1):ncol(data)] = c('deconv_trace', 'nevents')
  data$cell = as.integer(data$cell)
  
  return(data)
}


# Filter based on velocity calculated for the time window. Filtering is done on epochs
# Slow implementation, use isRunning to faster.
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

