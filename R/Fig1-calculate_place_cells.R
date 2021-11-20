library(dplyr)
library(plyr)
library(pracma)
library(readr)
library(tidyr)
library(stringr)
library(DT)
library(data.table)
library(tictoc)
library(pryr) # for memory checks
library(datatrace)


summarise = dplyr::summarise
summarize = dplyr::summarize

source('locations.R')
source('place_field_utils.R')
source('plotting_params.R')
source('traces_input.R')
source('utils.R')

nbins = 20
timebin.dur.msec = 200
min.occupancy.sec = 1.0
trace.var = 'smoothed_deconv_trace'
gen_imgs_dir = file.path(base_dir, .., 'pf_stability')

#all.trials.si = data.frame()
run.trials.si = data.frame()
odd.trials.si = data.frame()
even.trials.si = data.frame()
early.trials.si = data.frame()
late.trials.si = data.frame()
beforetest.trials.si = data.frame()
aftertest.trials.si = data.frame()

#all.fields = list()
#all.occupancies = list()
run.fields = list()
run.occupancies = list()
odd.fields = list()
odd.occupancies = list()
even.fields = list()
even.occupancies = list()
early.fields = list()
early.occupancies = list()
late.fields = list()
late.occupancies = list()
beforetest.fields = list()
beforetest.occupancies = list()
aftertest.fields = list()
aftertest.occupancies = list()

add.meta.cols = function(df, animal, date) {
  df$animal = as.factor(rep(animal, nrow(df)))
  df$date = as.factor(rep(date, nrow(df)))
  return(df)
}


for (caimg_result_dir in caimg_result_dirs) {
  tic("reading and preprocessing traces")
  data.traces = read.data.trace(caimg_result_dir)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  date = data.traces$date[1]
  animal = data.traces$animal[1]

  binned.traces.run = prepare.run.dirtraces(data.traces, nbins)
  toc()

  plot.dir.prefix = file.path(gen_imgs_dir, animal, format(date))

  # tic("spatial info on all trials")
  # all.spatial = calc.spatial.info(binned.traces[exp_title == 'trial'],
  #                                 plot.dir=paste0(plot.dir.prefix, '/all/'),
  #                                 generate.plots=FALSE,
  #                                 nshuffles=1000,
  #                                 timebin.dur.msec=timebin.dur.msec,
  #                                 nbins=nbins,
  #                                 min.occupancy.sec=min.occupancy.sec)
  # all.trials.si = bind_rows(all.trials.si, add.meta.cols(all.spatial$df, animal, date))
  # all.fields[[animal]][[format(date)]] = all.spatial$field
  # all.occupancies[[animal]][[format(date)]] = all.spatial$occupancy
  # toc()

  tic("spatial info on all trials running")
  run.spatial = calc.spatial.info(binned.traces.run[exp_title == 'trial'],
                                  plot.dir=paste0(plot.dir.prefix, '/run/'),
                                  generate.plots=FALSE,
                                  nshuffles=1000,
                                  timebin.dur.msec=timebin.dur.msec,
                                  shuffle.shift.sec = 10,
                                  trace.var=trace.var,
                                  nbins=nbins,
                                  min.occupancy.sec=min.occupancy.sec,
                                  gaussian.var = 2)
  run.trials.si = bind_rows(run.trials.si, add.meta.cols(run.spatial$df, animal, date))
  run.fields[[animal]][[format(date)]] = run.spatial$field
  run.occupancies[[animal]][[format(date)]] = run.spatial$occupancy
  toc()

  # Test trials
  #max_test_trial_dur_msec = 240 * 1000
  beforetest.traces = binned.traces.run[exp_title == 'beforetest'] # & timestamp <= max_test_trial_dur_msec]
  if (nrow(beforetest.traces) > 0) {
    tic("spatial info on beforetest run")
    test.spatial = calc.spatial.info(beforetest.traces,
                                     plot.dir=paste0(plot.dir.prefix, '/beforetest/'),
                                     generate.plots=FALSE,
                                     nshuffles=1000,
                                     timebin.dur.msec=timebin.dur.msec,
                                     shuffle.shift.sec = 20,
                                     trace.var=trace.var,
                                     nbins=nbins,
                                     min.occupancy.sec=min.occupancy.sec,
                                     gaussian.var = 2)
    beforetest.trials.si = bind_rows(beforetest.trials.si, add.meta.cols(test.spatial$df, animal, date))
    beforetest.fields[[animal]][[format(date)]] = test.spatial$field
    beforetest.occupancies[[animal]][[format(date)]] = test.spatial$occupancy
    toc()
  }

  aftertest.traces = binned.traces.run[exp_title == 'aftertest'] # & timestamp <= max_test_trial_dur_msec]
  if (nrow(aftertest.traces) > 0) {
    tic("spatial info on aftertest run")
    test.spatial = calc.spatial.info(aftertest.traces,
                                     plot.dir=paste0(plot.dir.prefix, '/aftertest/'),
                                     generate.plots=FALSE,
                                     nshuffles=1000,
                                     shuffle.shift.sec = 20,
                                     trace.var=trace.var,
                                     timebin.dur.msec=timebin.dur.msec,
                                     nbins=nbins,
                                     min.occupancy.sec=min.occupancy.sec,
                                     gaussian.var = 2)
    aftertest.trials.si = bind_rows(aftertest.trials.si, add.meta.cols(test.spatial$df, animal, date))
    aftertest.fields[[animal]][[format(date)]] = test.spatial$field
    aftertest.occupancies[[animal]][[format(date)]] = test.spatial$occupancy
    toc()
  }

  # # Odd vs Even
  # tic("spatial info on odd trials")
  # odd.spatial = calc.spatial.info(binned.traces.run[exp_title == 'trial' & trial %% 2 == 1,],
  #                                 paste0(plot.dir.prefix, '/odd/'),
  #                                 nshuffles=0,
  #                                 timebin.dur.msec=timebin.dur.msec,
  #                                 nbins=nbins,
  #                                 min.occupancy.sec=min.occupancy.sec,
  #                                 gaussian.var = 2)
  # odd.trials.si = bind_rows(odd.trials.si, add.meta.cols(odd.spatial$df, animal, date))
  # odd.fields[[animal]][[format(date)]] = odd.spatial$field
  # odd.occupancies[[animal]][[format(date)]] = odd.spatial$occupancy
  # toc()
  # 
  # tic("spatial info on even trials")
  # even.spatial = calc.spatial.info(binned.traces.run[exp_title == 'trial' & trial %% 2 == 0,],
  #                                  paste0(plot.dir.prefix, '/even/'),
  #                                  nshuffles=0,
  #                                  timebin.dur.msec=timebin.dur.msec,
  #                                  nbins=nbins,
  #                                  min.occupancy.sec=min.occupancy.sec,
  #                                  gaussian.var = 2)
  # even.trials.si = bind_rows(even.trials.si, add.meta.cols(even.spatial$df, animal, date))
  # even.fields[[animal]][[format(date)]] = even.spatial$field
  # even.occupancies[[animal]][[format(date)]] = even.spatial$occupancy
  # toc()

  # Early vs late
  half.trial = ceiling(max(data.traces$trial) / 2)
  tic("spatial info on early trials")
  early.spatial = calc.spatial.info(binned.traces.run[exp_title == 'trial' & trial <= half.trial,],
                                    paste0(plot.dir.prefix, '/early/'),
                                    nshuffles=0,
                                    trace.var=trace.var,
                                    timebin.dur.msec=timebin.dur.msec,
                                    nbins=nbins,
                                    min.occupancy.sec=min.occupancy.sec,
                                    gaussian.var = 2)
  early.trials.si = bind_rows(early.trials.si, add.meta.cols(early.spatial$df, animal, date))
  early.fields[[animal]][[format(date)]] = early.spatial$field
  early.occupancies[[animal]][[format(date)]] = early.spatial$occupancy
  toc()

  tic("spatial info on late trials")
  late.spatial = calc.spatial.info(binned.traces.run[exp_title == 'trial' & trial > half.trial, ],
                                   paste0(plot.dir.prefix, '/late/'),
                                   nshuffles=0,
                                   timebin.dur.msec=timebin.dur.msec,
                                   trace.var=trace.var,
                                   nbins=nbins,
                                   min.occupancy.sec=min.occupancy.sec,
                                   gaussian.var = 2)
  late.trials.si = bind_rows(late.trials.si, add.meta.cols(late.spatial$df, animal, date))
  late.fields[[animal]][[format(date)]] = late.spatial$field
  late.occupancies[[animal]][[format(date)]] = late.spatial$occupancy
  toc()

  print('Memory used')
  print(mem_used())
}


print("Saving env variables")
save.image(file="data/2021-01-12_pf_smooth_deconv_dfs_percentile_95_shuffle20sec_occupancy1sec.RData")

