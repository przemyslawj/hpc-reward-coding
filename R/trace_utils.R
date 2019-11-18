library(data.table)
library(foreach)


zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
}

sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))

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