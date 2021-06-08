library(dplyr)
library(data.table)
library(datatrace)
library(permute)

source('distances.R')
source('traces_input.R')
source('locations.R')
source('plotting_params.R')
source('fixed_effects.R')
source('utils.R')

mouse.meta.df = read.mouse.meta(rootdirs)
trials.meta.df = read.trials.meta(rootdirs)
all.locations.df = map_dfr(rootdirs, read_locations)
rewards.df =  all.locations.df %>%
  filter(!is_test) %>%
  left_join(mouse.meta.df, by='animal') %>%
  left_join(trials.meta.df, by=c('date', 'animal', 'location_set'))

nshuffles=30
ncells=30

xybins = 20
nbins = 20
cheeseboard.bin.cm = 120 / xybins

get.max.thr = function(max.fraction=0.25, absolute.min=0.01) {
  function(vals) {  
    max.val = max(c(vals, absolute.min), na.rm = TRUE)
    c( max.val* max.fraction, max.val+1) 
  }
}

get.nvenent.thr = function(min.nevents=0.1) {
  function(vals) { c(min.nevents, max(vals) + 1)}
}

prepare.traces = function(data.traces,
                          binned.var,
                          get.bin.thresholds.fun,
                          filter.running=TRUE,
                          timebin.dur.msec=200) {
  setorder(data.traces, exp_title, trial_id, cell_id, timestamp)
  data.traces = detect.events(data.traces, deconv.threshold=0.1)
  # running speed avg > 4 cm/s in 0.5 s window
  data.traces = add.running.col(data.traces, 3.3, 10)
  data.traces = gauss.smooth.df.var(data.traces, filter.len=20, sigma=4.0)
  data.traces[, `:=` (zscored_deconv_trace = zscore(deconv_trace), 
                      zscored_trace = zscore(trace),
                      zscored_smooth_deconv_trace = zscore(smoothed_deconv_trace)),
              by=.(exp_title, cell_id)]
  
  if (filter.running) {
    data.traces.filtered = data.traces[ x > 0 & y > 0 & is_running, ]
  } else {
    data.traces.filtered = data.traces[x > 0 & y > 0,]
  }
  
  date_str = format(data.traces$date[1])
  animal_name = data.traces$animal[1]
  
  binned.traces = bin.time.space(data.traces.filtered,
                                 nbins.x = xybins,
                                 nbins.y = xybins,
                                 get.bin.thresholds.fun = get.bin.thresholds.fun,
                                 binned.var=binned.var,
                                 timebin.dur.msec=timebin.dur.msec)
  binned.traces = gauss.smooth.df.var(binned.traces, var='nevents', out.var='smoothed_nevents',
                                      filter.len=5, sigma=1.0)
  
  return(binned.traces)
}


discrete.bayes.spatial.decoder = list(
  stim.var=quo(bin.xy),
  value.var=quo(response_bin),
  nstim.bins=xybins^2,
  train.fun = create.discrete.bayes,
  predict.fun = bayesmax,
  error.fun=bin.distance.error(xybins, xybins)
)

mfr.bayes.spatial.decoder = list(
  stim.var=quo(bin.xy),
  value.var=quo(smoothed_deconv_trace),
  nstim.bins=xybins^2,
  train.fun = create.mfr.bayes,
  predict.fun = bayesmax_poisson,
  error.fun=bin.distance.error(xybins, xybins)
)


read.binned.traces = function(caimg_result_dir, 
                              exp_title='trial') {
  data.traces = read.data.trace(caimg_result_dir, filter_exp_title = exp_title)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  prepare.traces(data.traces, 
                 binned.var='zscored_trace',
                 get.bin.thresholds.fun = get.quantiles.fun(c(0.9, 1.0)),
                 timebin.dur.msec = 200)
}
  
eval.decoder.same.day = function(binned.traces,
                                 decoder.config, 
                                 prepare.vars.fun=NULL,
                                 filter.test.traces.fun=NULL,
                                 min.samples=10,
                                 ncells=30, # If ncells == 0, take all cells
                                 nshuffles=20) {
  
  if (nrow(binned.traces) == 0) {
    return(data.frame())
  } 
  
  print('Started evaluating the model')
  if (!is.null(prepare.vars.fun)) {
    binned.traces = prepare.vars.fun(binned.traces)
  }
  cell_ids = binned.traces$cell_id %>% unique
  
  all.shuffles.eval.res = map_dfr(1:nshuffles, function(j) {
    subset.traces = binned.traces
    if (ncells > 0) {
      cell_ids = cell_ids[shuffle(cell_ids)]
      subset.traces = binned.traces[cell_id %in% cell_ids[1:min(ncells, length(cell_ids))],]
    }
    
    eval.df = eval.decoder(subset.traces,
                           nstim.bins=decoder.config$nstim.bins,
                           stim.var=!!decoder.config$stim.var,
                           error.fun=decoder.config$error.fun,
                           train.fun=decoder.config$train.fun,
                           predict.fun=decoder.config$predict.fun,
                           value.var=!!decoder.config$value.var,
                           cv=TRUE,
                           min.samples=min.samples,
                           filter.test.traces.fun=filter.test.traces.fun)
    
    if (nrow(eval.df) > 0) {
      eval.df$animal = binned.traces$animal[1]
      eval.df$date = binned.traces$date[1]
      eval.df$shuffle = j
    }
    eval.df
  })
  
  return(all.shuffles.eval.res)
}


scale.bins2cm = function(decoding.df) {
  decoding.df$error = decoding.df$error * cheeseboard.bin.cm
  if ('random_error' %in% names(decoding.df)) {
    decoding.df$random_error = decoding.df$random_error * cheeseboard.bin.cm
  }
  return(decoding.df)
}

join.meta.dfs = function(df) {
  df = left_join(df, mouse.meta.df, by='animal') %>%
    left_join(trials.meta.df, by=c("animal", "date"))
  df$day_desc = as.factor(df$day_desc)
  return(df)
}

eval.parts = list()
shuffle.configs.df = data.frame(ncells=c(0, 60, 30), nshuffles=c(1, 20, 30))

for (caimg_result_dir in habit_caimg_dirs) {
  binned.traces = read.binned.traces(caimg_result_dir)
  for (config_i in 1:nrow(shuffle.configs.df)) {
    out.name = paste(caimg_result_dir, config_i, sep='_')
    res.df = eval.decoder.same.day(
      binned.traces,
      decoder.config=discrete.bayes.spatial.decoder,
      #decoder.config=mfr.bayes.spatial.decoder,
      min.samples=10,
      ncells=shuffle.configs.df$ncells[config_i],
      nshuffles=shuffle.configs.df$nshuffles[config_i])
    res.df$decoder_ncells = shuffle.configs.df$ncells[config_i]
    eval.parts[[out.name]] = res.df
  }
}

spatial.error.df = do.call('rbind', eval.parts)
spatial.error.df = scale.bins2cm(spatial.error.df) %>% join.meta.dfs()

today_str = format(Sys.Date())
print("Saving env variables")
save.image(file=paste0("data/", today_str, "spatial_decoding.RData"))
# 
# eval.parts = list()
# for (caimg_result_dir in test_caimg_dirs) {
#   eval.parts[[caimg_result_dir]] = eval.decoder.same.day(
#     caimg_result_dir, 
#     decoder.config = discrete.bayes.spatial.decoder,
#     #decoder.config=mfr.bayes.spatial.decoder,
#     min.samples=5,
#     ncells=ncells,
#     nshuffles=nshuffles,
#     exp_title='beforetest')
# }
# 
# beforetest.spatial.error.df = do.call('rbind', eval.parts)
# beforetest.spatial.error.df = beforetest.spatial.error.df %>% 
#   scale.bins2cm %>%
#   join.meta.dfs()
           