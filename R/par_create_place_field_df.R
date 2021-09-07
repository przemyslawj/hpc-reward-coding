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
trace.var = 'smoothed_deconv_trace'


add.meta.cols = function(df, animal, date) {
  df$animal = as.factor(rep(animal, nrow(df)))
  df$date = as.factor(rep(date, nrow(df)))
  return(df)
}

traces2pf = function(binned.traces.run) {
  if (nrow(binned.traces.run) == 0) {
    return(list())
  }

  traces2pf.subset = function(subset.binned.traces, subset.name, shuffle.shift.sec = 10) {
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
                                      shuffle.shift.sec = shuffle.shift.sec,
                                      trace.var=trace.var,
                                      timebin.dur.msec=timebin.dur.msec,
                                      nbins=nbins,
                                      min.occupancy.sec=1,
                                      gaussian.var=2)
    subset.result$df = add.meta.cols(subset.result$df, animal, date)
    toc()
    return(subset.result)
  }

  run.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial'], 'run')

  # Test trials
  #max_test_trial_dur_msec = 240 * 1000
  beforetest.traces = binned.traces.run[exp_title == 'beforetest' | exp_title == 'test']# & timestamp <= max_test_trial_dur_msec]
  beforetest.pf = traces2pf.subset(beforetest.traces, 'beforetest', shuffle.shift.sec = 20)

  aftertest.traces = binned.traces.run[exp_title == 'aftertest']# & timestamp <= max_test_trial_dur_msec]
  aftertest.pf = traces2pf.subset(aftertest.traces, 'aftertest', shuffle.shift.sec = 20)

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

print(paste('Started processing the traces', Sys.time()))
# balance speed and memory, outfile "" should make the output print display in the console
cl = makeCluster(2, outfile="")
registerDoSNOW(cl)

daytraces.pf.list = foreach(caimg_result_dir=caimg_result_dirs,
                            .packages=c('datatrace', 'dplyr', "data.table", 'tictoc')) %dopar% {
  tic("reading and preprocessing traces")
  data.traces = read.data.trace(caimg_result_dir)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  binned.traces.run = prepare.run.dirtraces(data.traces, nbins, binned.var=trace.var)
  toc()

  date_str = format(binned.traces.run$date[1])
  animal_str = binned.traces.run$animal[1]
  day.pf = traces2pf(binned.traces.run)
  list(animal=animal_str, date=date_str, pfavl=day.pf)
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

  run.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$run$field
  run.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$run$occupancy

  early.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$early$field
  early.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$early$occupancy

  late.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$late$field
  late.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$late$occupancy

  beforetest.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$beforetest$field
  beforetest.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$beforetest$occupancy

  aftertest.fields[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$aftertest$field
  aftertest.occupancies[[animal]][[date]] = daytraces.pf.list[[i]]$pfavl$aftertest$occupancy
}

run.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfavl$run$df)
early.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfavl$early$df)
late.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfavl$late$df)
beforetest.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfavl$beforetest$df)
aftertest.trials.si = map_dfr(daytraces.pf.list, ~ .x$pfavl$aftertest$df)

print("Saving env variables")
save.image(file="data/2021-01-11-all_pf_smooth_deconv_dfs_percentile_95_bin200msec_nbins20_shuffle20sec_occupancy1sec_gaussvar2.RData")

