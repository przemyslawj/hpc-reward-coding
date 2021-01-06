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

# Plotting
library(ggplot2)
library(cowplot)

summarise = dplyr::summarise
summarize = dplyr::summarize

source('distances.R')
source('fixed_effects.R')
source('locations.R')
source('plotting_params.R')
source('place_field_utils.R')
source('utils.R')
source('traces_input.R')

mouse.meta.df = read.mouse.meta(rootdirs)
trials.meta.df = read.trials.meta(rootdirs)

habit3.days.df = filter(trials.meta.df, day_desc == 'habituation day#3') %>%
  dplyr::select(date, animal) %>%
  dplyr::distinct() %>%
  mutate(caimg_dir = map2(animal, date, ~ find.caimg.dir(caimg_result_dirs, .x, .y)))

beforetest.rewards.df = map_dfr(rootdirs, read_locations) %>%
  filter(is_test, exp_title == 'beforetest') %>%
  left_join(mouse.meta.df, by='animal') %>%
  left_join(dplyr::select(trials.meta.df, -location_set), by=c('date', 'animal'))
beforetest.rewards.df = data.table(beforetest.rewards.df)

perc2dist = 1.2
rew.dist.threshold = 20 / perc2dist

prepare.timebinned.run.traces = function(data.traces, timebin.dur.msec) {
  data.traces = data.traces[exp_title != 'homecage' & exp_title != 'aftertest',]
  
  data.traces = add.running.col(data.traces, 3.3, 10)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  setorder(data.traces, exp_title, trial_id, cell_id, timestamp)
  detect.events(data.traces, deconv.threshold=0.2)
  
  data.traces = gauss.smooth.df.var(data.traces, filter.len=20, sigma=5.0)

  timebinned.run = timebin.traces(data.traces[is_running & x > 0 & y > 0, ],
                                  timebin.dur.msec = timebin.dur.msec)
  return(timebinned.run)
}


stimbin.traces.xy = function(traces, nbins) {
  binned.traces = traces %>% stimbin.traces(x, nbins) %>% stimbin.traces(y, nbins) %>%
    as.data.table()
  binned.traces[, bin.xy := to_1dim(bin.x, bin.y, nbins)]
  binned.traces
}

# Draw timestamps count to match occupancy
# Returns DFs with one row per each timestamp, specyfing its spatial bin and count of timestamps
# that matches occupancy in that bin between the two traces. Returns one DF per the
# input trace.
timestamps.per.bin = function(fst.trace, snd.trace, down_nbins=10) {
  binned.fst.trace = stimbin.traces.xy(fst.trace[cell_id == fst.trace$cell_id[1]], down_nbins)
  fst.mfr = with(binned.fst.trace, create_mfr_model(bin.x, bin.y,
                                                    nbins_x=down_nbins,
                                                    nbins_y=down_nbins,
                                                    trace=trace,
                                                    minOccupancy=1))
  binned.snd.trace = stimbin.traces.xy(snd.trace[cell_id == snd.trace$cell_id[1]], down_nbins)
  snd.mfr = with(binned.snd.trace, create_mfr_model(bin.x, bin.y,
                                                    nbins_x=down_nbins,
                                                    nbins_y=down_nbins,
                                                    trace=trace,
                                                    minOccupancy=1))
  matching.occupancyMap = pmin(fst.mfr$occupancyMap, snd.mfr$occupancyMap)

  fst.timestamps.per.bin = binned.fst.trace[, .(timestamp, bin.xy, bin.x, bin.y)]
  fst.timestamps.per.bin[, down.ntimestamp := map2_int(bin.x, bin.y, ~ as.integer(matching.occupancyMap[.x, .y]))]

  snd.timestamps.per.bin = binned.snd.trace[, .(timestamp, bin.xy, bin.x, bin.y)]
  snd.timestamps.per.bin[, down.ntimestamp := map2_int(bin.x, bin.y, ~ as.integer(matching.occupancyMap[.x, .y]))]

  return(list(df1 = fst.timestamps.per.bin, df2 = snd.timestamps.per.bin, occupancyMap=matching.occupancyMap))
}

# Shuffles the timestamps and returns next batch of downsampled timestamps per each bin
next.down.timestamps = function(shuffle.timestamps) {
  shuffle.timestamps = shuffle.timestamps[sample(nrow(shuffle.timestamps)),]
  shuffle.timestamps[, timestamp_id := rowid(bin.xy)]
  down.timestamps = shuffle.timestamps[timestamp_id <= down.ntimestamp, timestamp]
  return(down.timestamps)
}

nbins = 20
ndownsample_shuffles = 100
min.occupancy.sec = 1.0
timebin.dur.msec = 200
habit3.days.df = as.data.table(habit3.days.df)
down.test.si = data.frame()
down.habit.si = data.frame()

beforetest.fields = list()
beforetest.occupancies = list()
habituation.fields = list()
habituation.occupancies = list()

#test_caimg_dirs = test_caimg_dirs[19]
for (caimg_result_dir in test_caimg_dirs) {
  data.traces = read.data.trace(caimg_result_dir)

  # Beforetest trace
  #max_test_trial_dur_msec = 240 * 1000
  #test.traces = data.traces[exp_title == 'beforetest' & timestamp <= max_test_trial_dur_msec]
  test.traces = data.traces[exp_title == 'beforetest']
  timebinned.test.traces.run = prepare.timebinned.run.traces(test.traces, timebin.dur.msec)

  # Habit trace
  animal_name = data.traces$animal[1]
  date_str = data.traces$date[1]
  habit_caimg_dir = habit3.days.df[animal == animal_name, caimg_dir][[1]]
  habit.data.traces = read.data.trace(habit_caimg_dir)
  timebinned.habit.traces.run = prepare.timebinned.run.traces(habit.data.traces, timebin.dur.msec)

  # Group timestamps by spatial bin
  timestamps.per.bin.res = timestamps.per.bin(timebinned.test.traces.run, timebinned.habit.traces.run)
  shuffle.test.timestamps = timestamps.per.bin.res$df1
  shuffle.habit.timestamps = timestamps.per.bin.res$df2

  get.bin.thresholds.fun = get.quantiles.fun(c(0.95, 1.0))
  binned.habit.traces.run = stimbin.traces.xy(timebinned.habit.traces.run, nbins)
  binned.habit.traces.run = bin.responses(binned.habit.traces.run, get.bin.thresholds.fun, binned.var='smoothed_deconv_trace')
  setkey(binned.habit.traces.run, timestamp, cell_id)

  binned.test.traces.run = stimbin.traces.xy(timebinned.test.traces.run, nbins)
  binned.test.traces.run = bin.responses(binned.test.traces.run, get.bin.thresholds.fun, binned.var='smoothed_deconv_trace')
  setkey(binned.test.traces.run, timestamp, cell_id)

  sample.spatial.info = function(timestamps.df, binned.traces, shuffle_i, nshuffles=200) {
    down.timestamps = next.down.timestamps(timestamps.df)
    down.binned.traces = binned.traces[timestamp %in% down.timestamps,]
    spatial.res = calc.spatial.info(down.binned.traces,
                                    nshuffles=nshuffles,
                                    timebin.dur.msec = timebin.dur.msec,
                                    trace.var = 'smoothed_deconv_trace',
                                    min.occupancy.sec=min.occupancy.sec,
                                    shuffle.shift.sec=10,
                                    nbins=nbins,
                                    gaussian.var=2)
    spatial.res$df$shuffle_i = shuffle_i

    beforetest.rewards.df$current_loc = FALSE
    fields.df = calc.field.peaks.info(binned.traces$date[1],
                                      binned.traces$animal[1],
                                      spatial.res$field,
                                      filter.rews.df(beforetest.rewards.df, 
                                                     timebinned.test.traces.run$date[1], 
                                                     timebinned.test.traces.run$animal[1]),
                                      max.rew.dist.thr = rew.dist.threshold)
    return(left_join(spatial.res$df, fields.df, by=c('cell_id')))
    #return(spatial.res$df)
  }

  habit.trials.si = map_dfr(1:ndownsample_shuffles,
                            ~ sample.spatial.info(shuffle.habit.timestamps, binned.habit.traces.run, .x))
  habit.trials.si$date = binned.habit.traces.run$date[1]
  habit.trials.si$animal = animal_name
  habit.trials.si$date.test = char2date(test.traces$date[1])
  down.habit.si = bind_rows(down.habit.si, habit.trials.si)

  test.trials.si = map_dfr(1:ndownsample_shuffles,
                           ~ sample.spatial.info(shuffle.test.timestamps, binned.test.traces.run, .x))
  # test.trials.si = calc.spatial.info(binned.test.traces.run,
  #                                    nshuffles=1000,
  #                                    timebin.dur.msec=timebin.dur.msec,
  #                                    trace.var = 'smoothed_deconv_trace',
  #                                    min.occupancy.sec=min.occupancy.sec,
  #                                    shuffle.shift.sec=10,
  #                                    nbins=nbins,
  #                                    gaussian.var=2)$df
  #test.trials.si = filter(beforetest.trials.si, animal == animal_name, date == test.traces$date[1])
  test.trials.si$date = char2date(test.traces$date[1])
  test.trials.si$animal = animal_name
  down.test.si = bind_rows(down.test.si, test.trials.si)
}

#down.habit.si$date.test = char2date(down.habit.si$date.test)
print("Saving env variables")
save.image(file="data/2020-12-31-downsampled_smoothed_deconv_bin200msec_nbins20_occupancy1sec_gaussvar2_dist15_shuffle_10s.RData")

