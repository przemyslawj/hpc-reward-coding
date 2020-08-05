library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(rlang)
library(data.table)

# Plotting
library(ggplot2)
library(cowplot)
source('plotting_params.R')

source('locations.R')
source('utils.R')
source('traces_input.R')
source('place_field_utils.R')

summarise = dplyr::summarise
summarize = dplyr::summarize

perc2dist = 1.2
goal.cell.max.dist = 20 / perc2dist

mouse.meta.df = read.mouse.meta(rootdirs)
trials.meta.df = read.trials.meta(rootdirs)

snd.moved.location.df = all.locations.df %>%
  filter(!is_test) %>%
  add_prev_locations(prev.loc.set.diff=1) %>% 
  filter(!current_loc, prev_loc_set==2) %>%
  ungroup() %>%
  dplyr::select(-date) %>%
  dplyr::distinct()

snd.moved.location.id = snd.moved.location.df %>%
  dplyr::mutate(loc_id = paste0(animal, '_', location_ordinal)) %>%
  dplyr::pull(loc_id)

thd.moved.location.df = all.locations.df %>%
  filter(!is_test, location_set == 2) %>% 
  dplyr::select(-date) %>%
  dplyr::distinct() %>% 
  dplyr::mutate(loc_id = paste0(animal, '_', location_ordinal),
                current_loc=FALSE) %>%
  filter(!(loc_id %in% snd.moved.location.id)) 

fourth.moved.location.df = all.locations.df %>%
  filter(!is_test, location_set == 3) %>% 
  dplyr::select(-date) %>%
  dplyr::distinct() %>% 
  dplyr::mutate(loc_id = paste0(animal, '_', location_ordinal),
                current_loc=FALSE) %>%
  filter(!(loc_id %in% thd.moved.location.df)) 
# DFs need to have columns:
# - implant
# - animal
# - day_desc
# - cell_id
# - signif.si
# - mindist.var

join.chars = function(animal, cell_id) {
  stringr::str_c(animal, cell_id, sep='_')
}

trial.field.dist2current.rew = trial.peaks.at.current.rew %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))
trial.field.dist2fst.moved = trial.peaks.at.fst.moved
trial.field.dist2snd.moved = create.peaks.df.at.locs(run.fields, snd.moved.location.df, filter(trial.days.df, exp=='habituation'))
trial.field.dist2thd.moved = create.peaks.df.at.locs(run.fields, thd.moved.location.df, filter(trial.days.df, exp=='habituation'))
trial.field.dist2moved.loc = bind_rows(
    dplyr::mutate(trial.field.dist2fst.moved, moved_loc = 'fst'),
    dplyr::mutate(trial.field.dist2snd.moved, moved_loc = 'snd'),
    dplyr::mutate(trial.field.dist2thd.moved, moved_loc = 'thd'),
  ) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))

beforetest.field.dist2current.rew = beforetest.peaks.at.current.rew %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))
beforetest.field.dist2fst.moved = beforetest.peaks.fst.moved
beforetest.field.dist2snd.moved = create.peaks.df.at.locs(beforetest.fields, snd.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2thd.moved = create.peaks.df.at.locs(beforetest.fields, thd.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2fourth.moved = create.peaks.df.at.locs(beforetest.fields, fourth.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2moved.loc = bind_rows(
    mutate(beforetest.field.dist2fst.moved, moved_loc = 'fst'),
    mutate(beforetest.field.dist2snd.moved, moved_loc = 'snd'),
    mutate(beforetest.field.dist2thd.moved, moved_loc = 'thd'),
    mutate(beforetest.field.dist2thd.moved, moved_loc = 'fourth')
  ) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))

mindist.var = quo(min.rew.dist)

# mindist.var = quo(min.rew.dist)
# trial.field.dist2current.rew = filter(run.trials.pc2rew.dist, location_set.activity == location_set.reward) %>%
#   dplyr::rename(date=date.activity, day_desc=day_desc.activity)
# trial.field.dist2fst.moved = merged.run.trials.pc %>%
#   dplyr::group_by(date, day_desc, animal) %>%
#   dplyr::mutate(
#     min.rew.dist=calc.min.rew.dist(filter.rews.df(fst.moved.location.df, date[1], animal[1]),
#                                    field.max.x, field.max.y)$rew.dist) %>%
#   ungroup()
# trial.field.dist2snd.moved = merged.run.trials.pc %>%
#   dplyr::group_by(date, day_desc, animal) %>%
#   dplyr::mutate(
#     min.rew.dist=calc.min.rew.dist(filter.rews.df(snd.moved.location.df, date[1], animal[1]),
#                                    field.max.x, field.max.y)$rew.dist) %>%
#   ungroup()
# 
# beforetest.field.dist2current.rew = pc.test.df %>% ungroup()
# beforetest.field.dist2fst.moved = pc.test.df %>%
#   dplyr::group_by(date, day_desc, animal) %>%
#   dplyr::mutate(
#     min.rew.dist=calc.min.rew.dist(filter.rews.df(fst.moved.location.df, date[1], animal[1]),
#                                    field.max.x, field.max.y)$rew.dist) %>%
#   ungroup()
# 
# beforetest.field.dist2snd.moved = pc.test.df %>%
#   dplyr::group_by(date, day_desc, animal) %>%
#   dplyr::mutate(
#     min.rew.dist=calc.min.rew.dist(filter.rews.df(snd.moved.location.df, date[1], animal[1]),
#                                    field.max.x, field.max.y)$rew.dist) %>%
#   ungroup()
# 
beforetest.field.dist2moved.loc = dplyr::mutate(beforetest.field.dist2moved.loc,
                                                day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))
beforetest.field.dist2current.rew = dplyr::mutate(beforetest.field.dist2current.rew,
                                                  day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))


# Extracts first letters and numbers, e.g. learning1 day#1 -> l1d1
get.shortname.from.day.desc = function(day_desc) {
  if (length(day_desc) > 1) {
    base::warning('Function not vectorized')
  }
  desc.words = stringr::str_split(day_desc, ' ')[[1]]
  fst.letters = stringr::str_sub(desc.words, 1, 1)
  fst.numbers = stringr::str_match(desc.words, '\\d')[,1]
  fst.numbers[is.na(fst.numbers)] = ''
  paste0(stringr::str_c(fst.letters, fst.numbers), collapse='')
}

get.stage.from.day.desc = function(day_desc) {                                                                                                                                   
  if (stringr::str_starts(day_desc, 'habituation')) {                              
    return('H1')                                                                   
  }                                                                                
  if (stringr::str_starts(day_desc, 'learning1 day#6 test')) {                     
    return('L1_test')                                                              
  }                                                                                
  if (stringr::str_starts(day_desc, 'learning1')) {                                
    return('L1')                                                                   
  }                                                                                
  if (stringr::str_starts(day_desc, 'learning2')) {                                
    return('L2')                                                                   
  }                                                                                
  if (stringr::str_starts(day_desc, 'learning3')) {                                
    return('L3')                                                                   
  }                                                                                
}   

get.activity.pattern.desc = function(clust_char) {
  if (stringr::str_detect(clust_char, '^\\d\\d222$')) {
    return('1-place_cell')
  }
  if (stringr::str_detect(clust_char, '^\\d\\d2\\d\\d$')) {
    return('2-reward_cell')
  }
  return('other')
}

plot.cluster.assignments = function(clust.df, 
                                    clust.distance=50, 
                                    pattern.desc.fun=get.activity.pattern.desc, 
                                    ...) {
  xvars = enquos(...)
  xvar_names = lapply(xvars, quo_name)
  
  max.ncells = unique(clust.df$cell_id) %>% max
  clust.df$animal = as.factor(clust.df$animal)
  clust.df = dplyr::mutate(clust.df,
                           cell_id=(as.numeric(animal) - 1) * (max.ncells + 1) + cell_id)
  
  my.clust.df = dplyr::select(clust.df, cell_id, exp, my.clust) %>%
    spread(exp, my.clust)
  if (! (all(levels(clust.df$exp) %in% colnames(my.clust.df)))) {
    return(ggplot())
  }

  ncells = unique(clust.df$cell_id) %>% length
  my.clust.df = dplyr::arrange(my.clust.df, !!!xvars) 
  my.clust.df$ordinal = 1:nrow(my.clust.df)
  
  assignment.groups = my.clust.df %>% 
    dplyr::group_by(!!!xvars) %>%
    dplyr::summarise(min.ordinal = min(ordinal),
                     max.ordinal = max(ordinal),
                     pct = round((max.ordinal - min.ordinal + 1)/ ncells * 100)) %>%
    dplyr::mutate(clust_char = paste0(!!!xvars)) %>%
    dplyr::mutate(group_desc = map_chr(clust_char, pattern.desc.fun))
  
  assignment.groups$group_desc = as.factor(assignment.groups$group_desc)
  assignment.groups$group_ordinal = 1:nrow(assignment.groups)
  
  assignment.points = gather(assignment.groups, key='exp', value='clust', 
                             unlist(xvar_names))
  
  assignment.points = assignment.points %>%
    arrange(clust, group_ordinal) %>%
    ddply(.(exp), mutate,
                cum.pct=cumsum(pct))
  
  assignment.points$group_ordinal = as.factor(assignment.points$group_ordinal)
  assignment.points$exp = factor(assignment.points$exp, levels=unlist(xvar_names))
  assignments.polypoints = assignment.points %>%
    mutate(cum.pct.start = cum.pct - pct) %>%
    gather(cum.pct, cum.pct.start, key='src', value='y') %>%
    arrange(clust_char, exp, src)
  
  assignments.polypoints$exp = as.numeric(assignments.polypoints$exp)
  assignments.polypoints = bind_rows(
    dplyr::mutate(assignments.polypoints, exp=exp-0.05),
    dplyr::mutate(assignments.polypoints, exp=exp+0.05)
  )
  ordered.polypoints = bind_rows(
    filter(assignments.polypoints, src=='cum.pct.start') %>% arrange(exp),
    filter(assignments.polypoints, src=='cum.pct') %>% arrange(desc(exp)))
  
  yintercept = max((filter(ordered.polypoints, clust==1))$y) + clust.distance/2
  ordered.polypoints %>%
    ggplot() +
    geom_hline(aes(yintercept = yintercept), linetype='dashed') +
    geom_polygon(aes(x=exp, 
                     y=ifelse(clust==1, y, y + clust.distance), 
                     #y=cum.pct,
                     fill=group_desc,
                     group=group_ordinal), alpha=0.7) +
    scale_fill_grey(start=0.0, end=0.8) +
    scale_x_continuous(breaks=seq_along(xvar_names), labels = xvar_names) +
    theme(axis.line.x = element_blank(),
          line = element_line(size = px2pt * 0.5),
          legend.position = 'top') +
    geom_text(x=1.5, y=10, label=paste0('ncells=',ncells)) +
    geom_text(x=2, y=-5, label='inactive') +
    geom_text(x=2, y=206, label='active') +
    xlab('') + ylab('Cells (%)') 
  
}

filter.cell.present.ntimes = function(df, ntimes, cast.var, value.var) {
  cast.var.arg = enquo(cast.var)
  value.var.arg = enquo(value.var)
  wide.df = df %>%
    dplyr::rename(cast.var = !!cast.var.arg) %>% 
    reshape2::dcast(animal + cell_id ~ cast.var, 
                    value.var=quo_name(value.var.arg)) 
  wide.df$times.present = apply(wide.df[,3:ncol(wide.df)], 1, 
                                FUN=function(vals) {sum(!is.na(vals))})
  wide.df %>%
    dplyr::filter(times.present >= ntimes) %>% 
    dplyr::select(-times.present) %>%
    reshape2::melt(id.vars=c('animal', 'cell_id'), 
                   variable.name=quo_name(cast.var.arg), 
                   value.name=quo_name(value.var.arg)) %>%
    dplyr::mutate(cell.present = !is.na(!!value.var.arg))
}

#########################################################################################
# Cluster cells by activity at translocated reward, median by the learning stage
#########################################################################################

# Group days by stages and show active change schematic on groups
min.nstages = 4

stage.peaks.df = trial.field.dist2fst.moved %>%
  dplyr::mutate(exp=map_chr(day_desc, get.stage.from.day.desc)) %>%
  group_by(animal, cell_id, exp) %>%
  dplyr::summarise(med.rew.dist=median(!!mindist.var, na.rm=TRUE),
                   ndays.stage=n()) %>%
  dplyr::mutate(active = as.integer(med.rew.dist <= goal.cell.max.dist))
stage.peaks.filtered = filter.cell.present.ntimes(stage.peaks.df, min.nstages, exp, active)

stage.peaks.filtered$my.clust = as.integer(stage.peaks.filtered$active) + 1
stage.peaks.filtered = left_join(stage.peaks.filtered, mouse.meta.df)

get.stage.activity.pattern.desc = function(clust_char) {
  if (stringr::str_detect(clust_char, '^\\d222$')) {
    return('1-place_cell')
  }
  if (stringr::str_detect(clust_char, '^\\d2\\d\\d$')) {
    return('2-reward_memory_cell')
  }
  return('other')
}

plot.cluster.assignments(filter(stage.peaks.filtered, 
                                #cell_id %in% stable.animal.cells$cell_id,
                                #animal=='G-BR'),
                                implant=='vCA1'),
                         clust.distance = 100,
                         pattern.desc.fun = get.stage.activity.pattern.desc,
                         H1,
                         L1,
                         L2,
                         L3
) +
  scale_x_continuous(labels=c('habituation','learning 1', 'learning 2', 'learning 3')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title='Reward activated place cells at first translocated location\n(stage median)')

 ########################################################################################################
# For the cells active at the translocated reward, were they active at reward before or after?
# View of how the peaks in the place field of reward-active cell followed the reward.

get.animal_cells = function(dist.df, cell.predicate) {
  dist.df %>%
    filter(!!cell.predicate) %>%
    filter(signif.si) %>%
    dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
    dplyr::pull(animal_cell)
}

get.rew.active.animal_cells = function(dist.df, dist.thr = goal.cell.max.dist) {
  get.animal_cells(dist.df, expr(!!mindist.var <= goal.cell.max.dist))
}

#############################################################################
# Distance to current reward for place cells active at fst translocated reward
#############################################################################
cells.transloc.rew.active.l15 = get.rew.active.animal_cells(
  filter(trial.field.dist2fst.moved, day_desc == 'learning1 day#5'))

selected.cells.peaks2current.rew.dist = bind_rows(
  dplyr::filter(trial.field.dist2current.rew, day_desc %in% c('learning2 day#2', 'learning3 day#2')),
  dplyr::filter(trial.field.dist2moved.loc, moved_loc == 'fst', day_desc %in% c('learning1 day#5', 'habituation day#3'))
) %>% 
  dplyr::mutate(active = !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l15)

filtered.cells.peaks2current.rew.dist = filter.cell.present.ntimes(selected.cells.peaks2current.rew.dist, 2, day_desc, !!mindist.var)
filtered.cells.peaks2current.rew.dist$compared_location = 'current reward'

filtered.cells.peaks2fst.rew.dist = 
  dplyr::filter(trial.field.dist2moved.loc, moved_loc == 'fst',
                day_desc %in% c('habituation day#3', 'learning1 day#5',
                                 'learning2 day#2', 'learning3 day#2')) %>%
  dplyr::mutate(active = !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l15) %>%
  filter.cell.present.ntimes(2, day_desc, !!mindist.var)
filtered.cells.peaks2fst.rew.dist$compared_location = 'first reward location'
  
left_join(bind_rows(filtered.cells.peaks2current.rew.dist, filtered.cells.peaks2fst.rew.dist), mouse.meta.df) %>%
  dplyr::mutate(unique_cell_id=join.chars(animal, cell_id)) %>%
  ggplot(aes(x=day_desc, y=!!mindist.var * perc2dist)) +
  geom_path(mapping=aes(group=unique_cell_id), alpha=0.2) +
  geom_violin(alpha=0.6) +
  stat_summary(fun.y=median, geom="point", size=2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(implant ~ compared_location) +
  geom_hline(yintercept = goal.cell.max.dist * perc2dist, linetype = 'dashed') +
  xlab('') + ylab('Distance to location (cm)') +
  scale_x_discrete(labels=c('habituation','learning1 day5','learning 2', 'learning 3'))

#############################################################################################
# Distance to current reward for place cells active at fst translocated reward - test probes
#############################################################################################
cells.transloc.rew.active.l2d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, day_desc == 'learning2 day#1 test', moved_loc == 'fst'))

# Baseline cells: inactive at the current rewards and not present at the future reward
cells.transloc.rew.inactive.l2d1t = get.animal_cells(
  filter(beforetest.field.dist2moved.loc, 
         day_desc == 'learning2 day#1 test',
         moved_loc %in% c('fst', 'snd', 'thd')) %>%
    dplyr::group_by(animal, cell_id) %>%
    dplyr::summarise(!!mindist.var := min(!!mindist.var), signif.si=any(signif.si)),
  expr(!!mindist.var > goal.cell.max.dist))

cells.transloc.rew.active.l3d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, day_desc == 'learning3 day#1 test', moved_loc == 'snd'))
cells.transloc.rew.inactive.l3d1t = get.animal_cells(
  filter(beforetest.field.dist2moved.loc, 
         day_desc == 'learning3 day#1 test',
         moved_loc %in% c('snd', 'thd', 'fourth')) %>%
    dplyr::group_by(animal, cell_id) %>%
    dplyr::summarise(!!mindist.var := min(!!mindist.var), signif.si=any(signif.si)),
  expr(!!mindist.var > goal.cell.max.dist))

# Distance to current reward for place cells active at translocated reward
selected.cells.peaks2current.rew.dist = bind_rows(
  dplyr::filter(trial.field.dist2moved.loc, day_desc == 'habituation day#3', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2moved.loc, day_desc == 'learning2 day#1 test', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2current.rew, day_desc %in% c('learning3 day#1 test', 'learning3 day#3 test'))
) %>% 
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l2d1t)

filtered.cells.peaks2current.rew.dist = filter.cell.present.ntimes(
  selected.cells.peaks2current.rew.dist, 2, day_desc, !!mindist.var) %>%
  dplyr::mutate(active.at.rew = !!mindist.var <= goal.cell.max.dist)
filtered.cells.peaks2current.rew.dist$compared_location = 'current reward'

filtered.cells.peaks2fst.rew.dist = bind_rows(
  dplyr::filter(trial.field.dist2moved.loc, day_desc == 'habituation day#3', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2moved.loc, 
                day_desc %in% c('learning2 day#1 test', 'learning3 day#1 test', 'learning3 day#3 test'),
                moved_loc == 'fst')) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l2d1t) %>%
  filter.cell.present.ntimes(2, day_desc, !!mindist.var) %>%
  dplyr::mutate(active.at.loc =  !!mindist.var <= goal.cell.max.dist)
filtered.cells.peaks2fst.rew.dist$compared_location = 'first reward location'

pc.beforetest.change.df = dplyr::mutate(pc.beforetest.change.df, 
                                        animal_cell = join.chars(animal, cell_id))
remapped.cells.l3d1 = filter(pc.beforetest.change.df,
                             day_desc.fst=='learning2 day#1',
                             day_desc.snd=='learning3 day#1',
                             n > 1,
                             field.cor < 0.3) %>%
  dplyr::pull(animal_cell)
remapped.cells.l3d3 = filter(pc.beforetest.change.df,
                             day_desc.fst=='learning3 day#1',
                             day_desc.snd=='learning3 day#3',
                              n > 1,
                              field.cor < 0.3) %>%
  dplyr::pull(animal_cell)
remapped.cells.fromL1toL3 = filter(pc.beforetest.change.df,
                             day_desc.fst=='learning2 day#1',
                             day_desc.snd=='learning3 day#3',
                             n > 1,
                             field.cor < 0.3) %>%
  dplyr::pull(animal_cell)
# remapped.cells.l3d1 = filter(filtered.cells.peaks2fst.rew.dist, !active.at.loc, day_desc == 'learning3 day#1 test') %>%
#   dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
#   dplyr::pull(animal_cell)
# 
# # This is wrong
# remapped.cells.l3d3 = filter(filtered.cells.peaks2fst.rew.dist, !active.at.loc, day_desc == 'learning3 day#3 test') %>%
#   dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
#   dplyr::pull(animal_cell)

left_join(bind_rows(filtered.cells.peaks2current.rew.dist, filtered.cells.peaks2fst.rew.dist), mouse.meta.df) %>%
  dplyr::mutate(unique_cell_id=join.chars(animal, cell_id)) %>%
  ggplot(aes(x=day_desc, y=!!mindist.var * perc2dist)) +
  geom_path(mapping=aes(group=unique_cell_id), alpha=0.5) +
  #geom_violin(alpha=0.6) +
  #stat_summary(fun.y=median, geom="point", size=2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(implant ~ compared_location) +
  geom_hline(yintercept = goal.cell.max.dist * perc2dist, linetype = 'dashed') +
  xlab('') + ylab('Distance to location (cm)') +
  scale_x_discrete(labels=c('habituation', 'learning1 probe', 'learning2 probe', 'learning3 probe'))
ggsave('/tmp/translocated_distance.svg', units = 'cm', width=8, height=8)

# DF with distance to the current reward, filtered to cells that remapped and 
# were reward-active during the previous beforetest
dist2current.rew.active = bind_rows(
  filter(beforetest.field.dist2current.rew, 
         day_desc == 'learning3 day#1 test',
         animal_cell %in% cells.transloc.rew.active.l2d1t) %>%
      dplyr::mutate(remapped = animal_cell %in% remapped.cells.l3d1),
  filter(beforetest.field.dist2current.rew, 
         day_desc == 'learning3 day#3 test',
         animal_cell %in% cells.transloc.rew.active.l3d1t) %>%
      dplyr::mutate(remapped = animal_cell %in% remapped.cells.l3d3))

# Summary DF with distance to the current reward filtered to the cells that remapped
# and were reward-active during the previous beforetest 
animal.learning2test.rew.following = dist2current.rew.active %>%
  filter(remapped) %>%
  dplyr::mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>% 
  group_by(implant, animal, day_desc) %>% 
  dplyr::summarise(at.rew = sum(is.at.rew), 
                   at.rew.pct = mean(is.at.rew) * 100, 
                   quantile30.dist = unname(quantile(!!mindist.var, 0.3)[1]),
                   med.dist = median(!!mindist.var), 
                   n=n(), .groups='drop') %>%
  arrange(animal, day_desc)

# Summary DF for distance to current reward, filtered to cells that remapped 
# and were reward-active two beforetests ago
animal.learning2test.rew.following2 = bind_rows(
  filter(beforetest.field.dist2current.rew,
         day_desc == 'learning3 day#3 test',
         animal_cell %in% get.rew.active.animal_cells(
            filter(beforetest.field.dist2moved.loc, 
                   day_desc == 'learning2 day#1 test', 
                   (moved_loc == 'fst' | moved_loc == 'snd'))),
         animal_cell %in% remapped.cells.fromL1toL3)
) %>%
  dplyr::mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>%
  group_by(implant, animal, day_desc) %>%
  dplyr::summarise(at.rew = sum(is.at.rew), 
                   at.rew.pct = mean(is.at.rew) * 100, 
                   med.dist = median(!!mindist.var), 
                   quantile30.dist = unname(quantile(!!mindist.var, 0.3)[1]),
                   n=n(), .groups='drop') %>%
  arrange(animal, day_desc)

# DF with distance to the current reward, filtered to cells that
# were reward-inactive during the previous beforetest
dist2current.rew.inactive = bind_rows(
    filter(beforetest.field.dist2current.rew,
           day_desc == 'learning3 day#1 test',
           animal_cell %in% cells.transloc.rew.inactive.l2d1t) %>%
        dplyr::mutate(remapped=animal_cell %in% remapped.cells.l3d1),
    filter(beforetest.field.dist2current.rew,
           day_desc == 'learning3 day#3 test',
           animal_cell %in% cells.transloc.rew.inactive.l3d1t,
           animal_cell %in% remapped.cells.l3d3) %>%
        dplyr::mutate(remapped=animal_cell %in% remapped.cells.l3d3))

dist2.current.rew.df = bind_rows(
  dplyr::mutate(dist2current.rew.active, group='active'),
  dplyr::mutate(dist2current.rew.inactive, group='inactive'))

# Cells pooled: the rew.active vs inactive have the same distance to the reward
wilcox.test(min.rew.dist ~ group,
            dist2.current.rew.df,
            subset = implant=='dCA1' & remapped,
            paired=FALSE,
            alternative='less')

rew.dist.model = lmerTest::lmer(log2(2 + min.rew.dist) ~ group + (1 + day_desc | animal),
                                dist2.current.rew.df,
                                REML=TRUE,
                                subset=implant=='dCA1' & remapped)
plot.model.diagnostics(rew.dist.model,
                       subset(dist2.current.rew.df, implant=='dCA1' & remapped)$animal,
                       subset(dist2.current.rew.df, implant=='dCA1' & remapped)$group)
summary(rew.dist.model)

dist2.current.rew.df %>%
  filter(remapped) %>%
  ggplot(aes(x=min.rew.dist * perc2dist, group=group)) +
  #geom_histogram(aes(fill=group, y=..count..), color='grey80') +
  geom_density(aes(linetype=group)) +
  facet_grid(. ~ implant) +
  geom_vline(xintercept = goal.cell.max.dist, linetype='dotted') +
  coord_flip() +
  xlab('Distance to closer of the current\nreward locations (cm)') + ylab('Density') + labs(linetype='') +
  gtheme
ggsave('/home/prez/tmp/cheeseboard/reward_active_dist_to_current.svg', width=9, height=4.7, units='cm')

beforetest.animal.peaks.at.rew.pct = dist2current.rew.inactive %>%
  filter(remapped) %>%
  dplyr::mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::group_by(animal, implant, day_desc) %>%
  dplyr::summarise(at.rew = sum(is.at.rew), 
                   at.rew.pct = mean(is.at.rew) * 100, 
                   med.dist = median(!!mindist.var), 
                   quantile30.dist = unname(quantile(!!mindist.var, 0.3)[1]),
                   n=n(), .groups='drop') %>%
  arrange(animal, day_desc)

# Reward-active and -inactive place cells equally stable
pc.beforetest.change.df = as.data.table(pc.beforetest.change.df)
pc.beforetest.change.df$rew.active.fst = FALSE
setDT(pc.beforetest.change.df)[day_desc.fst=='learning2 day#1' & day_desc.snd=='learning3 day#1',
                               rew.active.fst := (animal_cell %in% cells.transloc.rew.active.l2d1t)]
setDT(pc.beforetest.change.df)[day_desc.fst=='learning3 day#1' & day_desc.snd=='learning3 day#3',
                               rew.active.fst := (animal_cell %in% cells.transloc.rew.active.l3d1t)]
#wilcox.test(field.cor ~ rew.active.fst,
#            data=pc.beforetest.change.df,
#            subset=implant=='dCA1' & n > 1 & !is.na(rew.active.fst) & signif.si.fst,
#            alternative='greater')
field.stability.model = lmerTest::lmer(field.cor ~  rew.active.fst + (1 + day_desc.fst | animal), 
                                       data=pc.beforetest.change.df, 
                                       subset=implant=='dCA1' & n > 1 & !is.na(rew.active.fst) & signif.si.fst)
summary(field.stability.model)

pc.beforetest.change.df %>%
  filter(n>1, !is.na(field.cor), signif.si.fst) %>%
  ggplot() +
  geom_histogram(aes(x=field.cor, 
                     y=stat(ncount))) +
                     #y=stat(count / sum(count)))) +
  facet_grid(implant ~ rew.active.fst) +
  ggtitle('Field stability not different for reward-active and -inactive place cells') +
  gtheme + xlab('Place field correlation') + ylab('Cell count')

cells.comparison.mat = bind_rows(
  mutate(animal.learning2test.rew.following, cell_group='reward_cells'), 
  mutate(beforetest.animal.peaks.at.rew.pct, cell_group='other_cells')) %>%
  mutate(animal_beforetest_day = paste(animal, day_desc, sep='_') %>% str_replace(' test', '')) %>%
  # only keep trials when learnt
  #filter(animal_beforetest_day %in% kept.testtrials$animal_beforetest_day) %>%
  tidyr::pivot_wider(id_cols=c(implant, animal, day_desc), 
                     names_from=cell_group, 
                     values_from=c(at.rew, at.rew.pct, med.dist, quantile30.dist, n),
                     names_sep='.') %>%
  filter(!is.na(n.reward_cells), !is.na(n.other_cells))

cells.comparison.mat2 = bind_rows(
  mutate(animal.learning2test.rew.following2, cell_group='reward_cells'), 
  mutate(beforetest.animal.peaks.at.rew.pct, cell_group='other_cells')) %>%
  mutate(animal_beforetest_day = paste(animal, day_desc, sep='_') %>% str_replace(' test', '')) %>% 
  # only keep trials when learnt
  #filter(animal_beforetest_day %in% kept.testtrials$animal_beforetest_day) %>%
  tidyr::pivot_wider(id_cols=c(implant, animal, day_desc), 
                     names_from=cell_group, 
                     values_from=c(at.rew, at.rew.pct, med.dist, quantile30.dist, n),
                     names_sep='.') 

print('Estimate for % of the cells following the reward')
cells.comparison.mat.summary = group_by(cells.comparison.mat, implant) %>%
  dplyr::summarise(n.rew_cells = sum(n.reward_cells, na.rm=TRUE),
                   n.rew_cells_followed = sum(at.rew.reward_cells, na.rm=TRUE),
                   n.rew_cells_notfollowed = n.rew_cells - n.rew_cells_followed,
                   n.other_cells = sum(n.other_cells, na.rm=TRUE),
                   n.other_cells_followed = sum(at.rew.other_cells, na.rm=TRUE),
                   n.other_cells_notfollowed = n.other_cells - n.other_cells_followed,
                   rew_cells_followed.pct = n.rew_cells_followed / n.rew_cells,
                   other_cells_followed.pct= n.other_cells_followed / n.other_cells)
cells.comparison.mat.summary

dca1.following.mat = matrix(unlist(cells.comparison.mat.summary[1, c(3,4,6,7)]), ncol=2)
dimnames(dca1.following.mat) = list(c('at rew 2', 'not at rew 2'), c('at rew 1', 'not at rew 1'))
fisher.test(dca1.following.mat, alternative='greater')

print('T test on reward cells percent - higher number of cells moves to reward than in general population')
cells.following.pct = cells.comparison.mat %>%
  group_by(implant, animal) %>%
  dplyr::summarise(n.rew_cells = sum(n.reward_cells, na.rm=TRUE),
                   n.rew_cells_followed = sum(at.rew.reward_cells, na.rm=TRUE),
                   n.other_cells = sum(n.other_cells, na.rm=TRUE),
                   n.other_cells_followed = sum(at.rew.other_cells, na.rm=TRUE),
                   at.rew.pct.reward_cells = n.rew_cells_followed / n.rew_cells,
                   at.rew.pct.other_cells = n.other_cells_followed / n.other_cells) %>% 
  filter(n.rew_cells > 0) %>% # calculate percentages only when pop large enough
  dplyr::select(implant, animal, at.rew.pct.reward_cells, at.rew.pct.other_cells) %>%
  tidyr::pivot_longer(-c(implant, animal), 
                      names_to='cell_group', 
                      names_transform = list(cell_group= ~ str_remove(.x, 'at.rew.pct.')),
                      values_to='at.rew.pct')

t.test(at.rew.pct ~ cell_group,
       data=cells.following.pct,
       subset=implant=='dCA1',
       paired=TRUE,
       alternative = "less")

# Consider changing to mixed-effects
t.test(subset(cells.comparison.mat, implant=='dCA1')$med.dist.reward_cells, 
       subset(cells.comparison.mat, implant=='dCA1')$med.dist.other_cells,
       paired=TRUE,
       alternative = "less")

# dCA1: BF ~= 1/3 - insufficient/modest evidence, strong evidence when BF<1/10
# vCA1: BF ~= 1/5 - modest evidence
library(BayesFactor)
ttestBF(subset(cells.comparison.mat, implant=='vCA1')$med.dist.reward_cells, 
        subset(cells.comparison.mat, implant=='vCA1')$med.dist.other_cells,
        paired=TRUE,
        rscale=sqrt(2)/2,
        nullInterval=c(-Inf, 0))

# Test on quantiles
t.test(subset(cells.comparison.mat, implant=='dCA1')$quantile30.dist.reward_cells, 
       subset(cells.comparison.mat, implant=='dCA1')$quantile30.dist.other_cells,
       paired=TRUE,
       alternative = "less")

print('The count of reward following cells not different in dCA1 than in the vCA1')
learning2test.rew.following.summary = animal.learning2test.rew.following %>%
  dplyr::ungroup() %>%
  group_by(implant) %>%
  dplyr::summarise(implant.active=sum(at.rew), n=sum(n), implant.notactive=n-implant.active, mean(at.rew.pct)) %>% 
  arrange(implant)
reward.cells.comparison = 
  matrix(c(learning2test.rew.following.summary$implant.active[1],
           learning2test.rew.following.summary$implant.active[2],
           learning2test.rew.following.summary$implant.notactive[1],
           learning2test.rew.following.summary$implant.notactive[2]),
  nrow = 2,
  dimnames = list(c("dCA1", "vCA1"), c("reward", "non-reward") ))
fisher.test(reward.cells.comparison)

################
# Cluster plot
#################

get.reward.pattern.desc = function(clust_char) {
  if (stringr::str_detect(clust_char, '^\\d21$')) {
    return('2-place_cell')
  }
  if (stringr::str_detect(clust_char, '^\\d22$')) {
    return('1-reward_cell')
  }
  return('other')
}

selected.cells.peaks2current.rew.dist.active = filtered.cells.peaks2current.rew.dist %>%
  dplyr::filter(day_desc != 'learning3 day#3 test') %>% 
  dplyr::mutate(active=!!mindist.var <= goal.cell.max.dist,
                exp=map_chr(day_desc, get.shortname.from.day.desc)) %>%
  filter.cell.present.ntimes(3, exp, active)
selected.cells.peaks2current.rew.dist.active$my.clust = as.integer(selected.cells.peaks2current.rew.dist.active$active) + 1
selected.cells.peaks2current.rew.dist.active = left_join(selected.cells.peaks2current.rew.dist.active, mouse.meta.df, by=c('animal'))
plot.cluster.assignments(filter(selected.cells.peaks2current.rew.dist.active, 
                                #animal=='F-TL'),
                                implant=='dCA1'),
                         clust.distance = 100,
                         pattern.desc.fun = get.reward.pattern.desc,
                         hd3, 
                         l2d1t, 
                         l3d1t
                       ) +
  labs(title='Changes in translocated active cells location')

#TODOs regenerate the peaks and present the figures..
# name the groups: group_desc
# place_cell, reward cell, reward silenced cell, reward_move activated
