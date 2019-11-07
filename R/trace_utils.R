zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
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