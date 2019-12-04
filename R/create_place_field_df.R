library(dplyr)
library(plyr)
library(pracma)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(DT)
library(data.table)
library(tictoc)
library(pryr) # for memory checks

# Plotting
library(ggplot2)
library(cowplot)
gtheme = theme_minimal() + theme_cowplot() +
  theme(text = element_text(size=18), axis.text = element_text(size=15))

summarise = dplyr::summarise
summarize = dplyr::summarize

source('locations.R')
source('place_field.R')
source('trace_utils.R')
source('utils.R')

root_dir = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-08/'
gen_imgs_dir = '/mnt/DATA/Prez/pf_stability/'

caimg_result_dirs = c(
  get.subject.result.dirs(root_dir, 'E-BL'),
  get.subject.result.dirs(root_dir, 'E-TR'),
  get.subject.result.dirs(root_dir, 'F-BL'),
  get.subject.result.dirs(root_dir, 'F-TL')
)


caimg_result_dirs = Filter(
  function(ca_img_dir) {file.exists(paste(ca_img_dir, 'traces_and_positions.csv', sep='/'))},
  caimg_result_dirs)
#caimg_result_dirs = c('/mnt/DATA/Prez/cheeseboard-down/down_2/2019-08/habituation/2019-08-29/caiman/F-BL/filtered/')


calc.spatial.info = function(data.traces, plot.dir='/tmp/pf_stability/',
                             generate.plots=FALSE, nshuffles=0) {
  frame.hz = 20
  data.traces = data.table(data.traces)
  cells = unique(data.traces$cell_id) %>% sort
  pci.df = data.frame()
  fields = list()
  occupancies = list()

  for (cell_name in cells) {
    cell.df = data.traces[cell_id == cell_name & smooth_trans_x >= 0 & smooth_trans_y >= 0,]
    pf = cell.spatial.info(cell.df, generate.plots, nshuffles, trace.col='trace')

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

  return(list(df=pci.df, field=fields, occupancy=occupancies))
}

run.threshold = 2

all.trials.si = data.frame()
run.trials.si = data.frame()
odd.trials.si = data.frame()
even.trials.si = data.frame()
early.trials.si = data.frame()
late.trials.si = data.frame()
test.trials.si = data.frame()

all.fields = list()
all.occupancies = list()
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
test.fields = list()
test.occupancies = list()

add.meta.cols = function(df, animal, date) {
  df$animal = as.factor(rep(animal, nrow(df)))
  df$date = as.factor(rep(date, nrow(df)))
  return(df)
}


for (ca_img_result_dir in caimg_result_dirs) {
  print(paste('Processing dir: ', ca_img_result_dir))
  cellmapping_file = file.path(ca_img_result_dir, 'cell_mapping.csv')
  cellmapping.df = read.csv(cellmapping_file)
  traces_file = file.path(ca_img_result_dir, 'traces_and_positions.csv')
  data = fread(traces_file)
  data = data[exp_title != 'homecage',]

  data.traces = melt.traces(data)
  data.traces = left_join(data.traces, cellmapping.df, by=c('cell'='cell_no'))
  dir_parts = str_split(ca_img_result_dir, '/')
  animal = dir_parts[[1]][length(dir_parts[[1]]) - 2]
  data.traces$animal = rep(animal, nrow(data.traces))
  date = data.traces$date[1]
  data.traces = data.table(data.traces)
  data.traces.run = data.traces[velocity >= run.threshold,]

  print('Analysing spatial information')
  plot.dir.prefix = paste(gen_imgs_dir, animal, format(date), sep='/')

  #tic("spatial info on all trials")
  #all.spatial = calc.spatial.info(data.traces[exp_title == 'trial'],
  #                                plot.dir=paste0(plot.dir.prefix, '/all/'),
  #                                generate.plots=TRUE,
  #                                nshuffles=100)
  #all.trials.si = bind_rows(all.trials.si, add.meta.cols(all.spatial$df, animal, date))
  #all.fields[[animal]][[format(date)]] = all.spatial$field
  #all.occupancies[[animal]][[format(date)]] = all.spatial$occupancy
  #toc()

  tic("spatial info on all trials running")
  run.spatial = calc.spatial.info(data.traces.run[exp_title == 'trial'],
                                  plot.dir=paste0(plot.dir.prefix, '/run/'),
                                  generate.plots=TRUE,
                                  nshuffles=100)
  run.trials.si = bind_rows(run.trials.si, add.meta.cols(run.spatial$df, animal, date))
  run.fields[[animal]][[format(date)]] = run.spatial$field
  run.occupancies[[animal]][[format(date)]] = run.spatial$occupancy
  toc()

  test.traces = data.traces.run[exp_title == 'test']
  if (nrow(test.traces) > 0) {
    tic("spatial info on test run")
    test.spatial = calc.spatial.info(test.traces,
                                     plot.dir=paste0(plot.dir.prefix, '/testrun/'),
                                     generate.plots=TRUE,
                                     nshuffles=100)
    test.trials.si = bind_rows(test.trials.si, add.meta.cols(test.spatial$df, animal, date))
    test.fields[[animal]][[format(date)]] = test.spatial$field
    test.occupancies[[animal]][[format(date)]] = test.spatial$occupancy
    toc()
  }

  # Odd vs Even
  tic("spatial info on odd trials")
  odd.spatial = calc.spatial.info(data.traces.run[exp_title == 'trial' & trial %% 2 == 1,],
                                  paste0(plot.dir.prefix, '/odd/'),
                                  nshuffles=0)
  odd.trials.si = bind_rows(odd.trials.si, add.meta.cols(odd.spatial$df, animal, date))
  odd.fields[[animal]][[format(date)]] = odd.spatial$field
  odd.occupancies[[animal]][[format(date)]] = odd.spatial$occupancy
  toc()

  tic("spatial info on even trials")
  even.spatial = calc.spatial.info(data.traces.run[exp_title == 'trial' & trial %% 2 == 0,],
                                   paste0(plot.dir.prefix, '/even/'),
                                   nshuffles=0)
  even.trials.si = bind_rows(even.trials.si, add.meta.cols(even.spatial$df, animal, date))
  even.fields[[animal]][[format(date)]] = even.spatial$field
  even.occupancies[[animal]][[format(date)]] = even.spatial$occupancy
  toc()

  # Early vs late
  half.trial = ceiling(max(data.traces$trial) / 2)
  tic("spatial info on early trials")
  early.spatial = calc.spatial.info(data.traces.run[exp_title == 'trial' & trial <= half.trial,],
                                    paste0(plot.dir.prefix, '/early/'),
                                    nshuffles=0)
  early.trials.si = bind_rows(early.trials.si, add.meta.cols(early.spatial$df, animal, date))
  early.fields[[animal]][[format(date)]] = early.spatial$field
  early.occupancies[[animal]][[format(date)]] = early.spatial$occupancy
  toc()

  tic("spatial info on late trials")
  late.spatial = calc.spatial.info(data.traces.run[exp_title == 'trial' & trial > half.trial, ],
                                   paste0(plot.dir.prefix, '/late/'),
                                   nshuffles=0)
  late.trials.si = bind_rows(late.trials.si, add.meta.cols(late.spatial$df, animal, date))
  late.fields[[animal]][[format(date)]] = late.spatial$field
  late.occupancies[[animal]][[format(date)]] = late.spatial$occupancy
  toc()

  print('Memory used')
  print(mem_used())
}


print("Saving env variables")
save.image(file="place_field_dfs_short.RData")

