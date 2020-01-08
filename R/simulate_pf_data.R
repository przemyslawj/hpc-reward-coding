source('trace_utils.R')
source('bayesian_decoder.R')


#subset.dirs = caimg_result_dirs[shuffle(caimg_result_dirs)][1:4]
sim.eval.parts = lapply(caimg_result_dirs, function(caimg_result_dir) {
  data.traces = read.data.trace(caimg_result_dir, filter_exp_title = 'trial')
  setorder(data.traces, trial_id, cell_id, timestamp)
  detect.events(data.traces, deconv.threshold=0.1)
  
  running.index = isRunning(data.traces, 2, 3, 500)
  data.traces.run = data.traces[which(running.index), ]
  
  date_str = data.traces$date[1]
  animal_name = data.traces$animal[1]
  timebinned.traces = timebin.traces(data.traces.filtered, timebin.dur.msec=200, xybins=20, trace.var=trace)
  quantile.fractions = c(0.95, 1.0)
  nresponse.bins = length(quantile.fractions)
  binned.traces = bin.responses(timebinned.traces, quantile.fractions)
  binned.traces = binned.traces %>%
    mutate(bin.xy = to_1dim(bin.x, bin.y))
  
  stimulus = binned.traces %>%
    select(time_bin, bin.xy) %>%
    arrange(time_bin) %>%
    distinct()
  
  model.response.matrix = reshape2::acast(binned.traces, cell_id ~ time_bin, value.var='response_bin')
  model.bayes = create.model2(model.response.matrix, stimulus$bin.xy, nresponse.bins, nstim.bins=xybins*xybins)
  
  #create data.frame with simulated response bins based on the model
  simulate.data.frame = function(model.bayes, binned.traces) {
    simulated.binned.traces = data.frame(binned.traces)
    x = runif(nrow(simulated.binned.traces))
    for (i in 1:nrow(simulated.binned.traces)) {
      cell_name = binned.traces$cell_id[i]
      probs = model.bayes$likelihood[binned.traces$bin.xy[i], paste(cell_name),]
      #probs = model.bayes$likelihood[binned.traces$bin.xy[i], 1,]
      response_bins = which(x[i] <= cumsum(probs) )
      simulated.binned.traces$response_bin[i] = response_bins[[1]]
    }
    simulated.binned.traces
  }
  system.time(simulated.binned.traces <- simulate.data.frame(model.bayes, binned.traces))
  
  sim.eval.df = eval.decoder(simulated.binned.traces, nresponse.bins=nresponse.bins, cv=TRUE)
  sim.eval.df
})

sim.decoding.df = do.call('rbind', sim.eval.parts)

summary(sim.decoding.df$error)
