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
gtheme = theme_minimal() + theme_cowplot() +
  theme(text = element_text(size=10), axis.text = element_text(size=8))

summarise = dplyr::summarise
summarize = dplyr::summarize

source('locations.R')
source('utils.R')

nbins = 20
timebin.dur.msec = 200

root_dir07 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-07/'
root_dir08 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-08/'
root_dir01 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2020-01/'
rootdirs = c(root_dir07, root_dir08, root_dir01)
gen_imgs_dir = '/mnt/DATA/Prez/pf_stability/'

caimg_result_dirs = c(
  get.subject.result.dirs(root_dir08, 'E-BL'),
  get.subject.result.dirs(root_dir08, 'E-TR'),
  get.subject.result.dirs(root_dir08, 'F-BL'),
  get.subject.result.dirs(root_dir08, 'F-TL'),
  get.subject.result.dirs(root_dir07, 'B-BL'),
  get.subject.result.dirs(root_dir07, 'C-1R'),
  get.subject.result.dirs(root_dir07, 'D-BR'),
  get.subject.result.dirs(root_dir01, 'G-BR'),
  get.subject.result.dirs(root_dir01, 'K-BR'),
  get.subject.result.dirs(root_dir01, 'L-TL')
)


caimg_result_dirs = Filter(
  function(ca_img_dir) {file.exists(paste(ca_img_dir, 'traces_and_positions.csv', sep='/'))},
  caimg_result_dirs)

calc.spatial.info = function(binned.traces, plot.dir='/tmp/pf_stability/',
                             generate.plots=FALSE, nshuffles=0) {
  binned.traces = data.table(binned.traces)
  cells = unique(binned.traces$cell_id) %>% sort
  pci.df = data.frame()
  fields = list()
  occupancies = list()

  for (cell_name in cells) {
    cell.df = binned.traces[cell_id == cell_name ,]
    pf = cell.spatial.info(cell.df, nbins, nbins, generate.plots, nshuffles, trace.var='trace',
                           bin.hz=1000/timebin.dur.msec,
                           shuffle.shift.sec=5,
                           min.occupancy.sec=1.5,
                           kernel.size = 9,
                           gaussian.var = 2)
    if (length(pf$cell_info) > 0) {
      fields[[format(cell_name)]] = pf$field
      occupancies[[format(cell_name)]] = pf$occupancy
      pci.df = bind_rows(pci.df, pf$cell_info)

      if (!is.na(plot.dir) && generate.plots && !is.na(pf$g)) {
        if (!file.exists(plot.dir)) {
          dir.create(plot.dir, recursive=TRUE)
        }
        ggsave(paste0(plot.dir, 'place_field_cell_', cell_name, '.jpg'), pf$g,
               width=3.5, height=3.0, units='cm', dpi=300)
      }
    }
  }

  return(list(df=pci.df, field=fields, occupancy=occupancies))
}

run.threshold = 2

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
                                      nshuffles=1000)
    subset.result$df = add.meta.cols(subset.result$df, animal, date)
    toc()
    return(subset.result)
  }

  run.pf = traces2pf.subset(binned.traces.run[exp_title == 'trial'], 'run')

  # Test trials
  max_test_trial_dur_msec = 240 * 1000
  beforetest.traces = binned.traces.run[exp_title == 'beforetest' & timestamp <= max_test_trial_dur_msec]
  beforetest.pf = traces2pf.subset(beforetest.traces, 'beforetest')

  aftertest.traces = binned.traces.run[exp_title == 'aftertest' & timestamp <= max_test_trial_dur_msec]
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

read.traces.for.dir = function(caimg_result_dir) {
  data.traces = read.data.trace(caimg_result_dir)
  data.traces = data.traces[exp_title != 'homecage',]
  detect.events(data.traces, deconv.threshold=0.2)

  setorder(data.traces, trial_id, cell_id, timestamp)
  running.index = isRunning(data.traces, 2, 4, 500)

  response.bin.quantiles = c(0.95, 1.0)
  binned.traces.run = bin.time.space(data.traces[which(running.index) & x >= 0 & y >= 0, ],
                                     nbins.x = nbins,
                                     nbins.y = nbins,
                                     get.bin.thresholds.fun = get.quantiles.fun(response.bin.quantiles),
                                     binned.var='trace',
                                     timebin.dur.msec=timebin.dur.msec)
  return(binned.traces.run)
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
  binned.traces.run = read.traces.for.dir(caimg_result_dir)
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
save.image(file="data/par_place_field_dfs_smoothed_percentile_95_shuffle5sec_occupancy1_5sec.RData")

