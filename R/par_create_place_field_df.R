library(dplyr)
library(plyr)
library(pracma)
library(readr)
library(tidyr)
library(stringr)
library(DT)
library(data.table)
library(tictoc)
library(datatrace)
library(foreach)
library(doSNOW)

setwd("~/mnt_code/cheeseboard_analysis/R")

summarise = dplyr::summarise
summarize = dplyr::summarize

source('locations.R')
source('place_field_utils.R')
source('plotting_params.R')
source('traces_input.R')
source('utils.R')

nbins = 20
timebin.dur.msec = 200
gen_imgs_dir = '/mnt/DATA/Prez/pf_stability/'


add.meta.cols = function(df, animal, date) {
  df$animal = as.factor(rep(animal, nrow(df)))
  df$date = as.factor(rep(date, nrow(df)))
  return(df)
}

traces2pf = function(binned.traces.run) {
  if (nrow(binned.traces.run) == 0) {
    return(list())
  }

  traces2pf.subset = function(subset.binned.traces, subset.name) {
    if (nrow(subset.binned.traces) == 0) {
      return(list())
    }

    subset.result = list()
    animal = binned.traces.run$animal[1]
    date = binned.traces.run$date[1]
    plot.dir.prefix = paste(gen_imgs_dir, animal, format(date), sep='/')

    tic(paste0("spatial info on subset=", subset.name))
    subset.result = calc.spatial.info(subset.binned.traces,
                                      plot.dir=paste(plot.dir.prefix, subset.name, sep='/'),
                                      generate.plots=FALSE,
                                      nshuffles=1000,
									  trace.var='deconv_trace',
                                      timebin.dur.msec=timebin.dur.msec,
                                      nbins=nbins,
                                      min.occupancy.sec=0.7)
    subset.result$df = add.meta.cols(subset.result$df, animal, date)
    toc()
    return(subset.result)
  }

  run.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial'], 'run')

  # Test trials
  #max_test_trial_dur_msec = 240 * 1000
  beforetest.traces = binned.traces.run[exp_title == 'beforetest']# & timestamp <= max_test_trial_dur_msec]
  beforetest.pf = traces2pf.subset(beforetest.traces, 'beforetest')

  aftertest.traces = binned.traces.run[exp_title == 'aftertest']# & timestamp <= max_test_trial_dur_msec]
  aftertest.pf = traces2pf.subset(aftertest.traces, 'aftertest')

  # # Odd vs Even
  # odd.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial' & trial %% 2 == 1,], 'odd')
  # even.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial' & trial %% 2 == 0,], 'even')

  # Early vs late
  half.trial = max(binned.traces.run$trial) / 2
  early.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial' & trial <= half.trial,], 'early')
  late.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial' & trial > half.trial,], 'late')

  return(list(run=run.pf,
              beforetest=beforetest.pf,
              aftertest=aftertest.pf,
              #odd=odd.pf,
              #even=even.pf,
              early=early.pf,
              late=late.pf))
}

#daytraces.pf.list = list()
#for (caimg_result_dir in caimg_result_dirs) {

print(paste('Started processing the traces', Sys.time()))
# balance speed and memory, outfile "" should make the output print display in the console
cl = makeCluster(3, outfile="")
registerDoSNOW(cl)

daytraces.pf.list = foreach(caimg_result_dir=caimg_result_dirs,
                            .packages=c('datatrace', 'dplyr', "data.table", 'tictoc')) %dopar% {
  tic("reading and preprocessing traces")
  data.traces = read.data.trace(caimg_result_dir)
  binned.traces.run = prepare.run.dirtraces(data.traces, nbins)
  toc()

  date_str = format(binned.traces.run$date[1])
  animal_str = binned.traces.run$animal[1]
  day.pf = traces2pf(binned.traces.run)
  list(animal=animal_str, date=date_str, pfval=day.pf)
}

stopCluster(cl)
print(paste('Finished processing the traces', Sys.time()))

run.fields = list()
run.occupancies = list()
early.fields = list()
early.occupancies = list()
late.fields = list()
late.occupancies = list()
beforetest.fields = list()
beforetest.occupancies = list()
aftertest.fields = list()
aftertest.occupancies = list()

for (i in 1:length(daytraces.pf.list)) {
  animal = daytraces.pf.list[[i]]$animal
  date = daytraces.pf.list[[i]]$date

  run.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$run$field
  run.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$run$occupancy

  early.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$early$field
  early.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$early$occupancy

  late.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$late$field
  late.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$late$occupancy

  beforetest.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$beforetest$field
  beforetest.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$beforetest$occupancy

  aftertest.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$aftertest$field
  aftertest.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfval$aftertest$occupancy
}

run.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfval$run$df)
early.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfval$early$df)
late.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfval$late$df)
beforetest.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfval$beforetest$df)
aftertest.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfval$aftertest$df)

print("Saving env variables")
save.image(file="data/deconv_place_field_dfs_smoothed_percentile_95_bin200msec_nbins20_shuffle5sec_occupancy07sec_gaussvar2.RData")

