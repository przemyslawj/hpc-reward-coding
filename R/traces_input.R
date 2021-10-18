library(data.table)
library(dplyr)
library(stringr)

source('locations.R')
source('utils.R')

root_dir07 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-07/'
root_dir08 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-08/'
root_dir01 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2020-01/'
root_dir10 = '/mnt/DATA/Prez/cheeseboard-down/down_2/2020-10/'
rootdirs = c(root_dir07, root_dir08, root_dir01, root_dir10)
#rootdirs = c(root_dir10)
tracking_rootdirs = stringr::str_replace(rootdirs, '-down/down_2', '')
#tracking_rootdirs = rootdirs

animal_names = c('O-TR', 'R-TR', 'M-BR', 'N-BR', 'N-BL',
                 'E-BL', 'E-TR', 'F-BL', 'F-TL',
                 'G-BR', 'K-BR', 'B-BL', 'D-BR')

caimg_result_dirs = c(
  get.subject.result.dirs(root_dir10, 'O-TR'),
  get.subject.result.dirs(root_dir10, 'R-TR'),
  get.subject.result.dirs(root_dir10, 'M-BR'),
  get.subject.result.dirs(root_dir10, 'N-BR'),
  get.subject.result.dirs(root_dir10, 'N-BL'),
  get.subject.result.dirs(root_dir08, 'E-BL'),
  get.subject.result.dirs(root_dir08, 'E-TR'),
  get.subject.result.dirs(root_dir08, 'F-BL'),
  get.subject.result.dirs(root_dir08, 'F-TL'),
  get.subject.result.dirs(root_dir01, 'G-BR'),
  get.subject.result.dirs(root_dir01, 'K-BR'),
  #get.subject.result.dirs(root_dir01, 'L-TL'), # misplaced lens - targeting ca3
  get.subject.result.dirs(root_dir07, 'B-BL'),
  #get.subject.result.dirs(root_dir07, 'C-1R'),
  get.subject.result.dirs(root_dir07, 'D-BR')
)

caimg_result_dirs = Filter(
  function(ca_img_dir) {file.exists(paste(ca_img_dir, 'traces_and_positions.csv', sep='/'))},
  caimg_result_dirs)

locations.df = map_dfr(rootdirs, read_locations) 
test.days.df = locations.df %>%
  filter(exp_title == 'beforetest') %>%
  filter(!animal %in% c('A-BL', 'C-1R', 'L-TL', 'P-BR')) %>%
  dplyr::select(animal, date) %>%
  dplyr::distinct()

test_caimg_dirs = map2(test.days.df$animal, test.days.df$date, ~ find.caimg.dir(caimg_result_dirs, .x, .y))
test_caimg_dirs = Filter(function(x) !is.na(x), test_caimg_dirs) 
habit_caimg_dirs = Filter(function(x) str_detect(x, 'habit'), caimg_result_dirs)
learning_caimg_dirs = Filter(function(x) str_detect(x, 'learning'), caimg_result_dirs)

prepare.run.dirtraces = function(data.traces, 
                                 nbins, 
                                 binned.var='trace',
                                 event.deconv.threshold=0.2) {
  data.traces = data.traces[exp_title != 'homecage',]
  detect.events(data.traces, deconv.threshold=event.deconv.threshold)

  # running speed avg > 4 cm/s in 0.5 s window
  data.traces = add.running.col(data.traces, 3.3, 10)
  data.traces = gauss.smooth.df.var(data.traces, filter.len=10, sigma=1.5)

  response.bin.quantiles = c(0.95, 1.0)
  binned.traces.run = bin.time.space(data.traces[is_running & x > 0 & y > 0, ],
                                     nbins.x = nbins,
                                     nbins.y = nbins,
                                     get.bin.thresholds.fun = get.quantiles.fun(response.bin.quantiles),
                                     binned.var=binned.var,
                                     timebin.dur.msec=timebin.dur.msec)
  return(binned.traces.run)
}

