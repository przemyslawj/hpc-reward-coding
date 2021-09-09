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
source('decoder_utils.R')

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
           