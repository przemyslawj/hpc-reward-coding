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

summarise = dplyr::summarise
summarize = dplyr::summarize

mouse.meta.df = read.mouse.meta(rootdirs)
trials.meta.df = read.trials.meta(rootdirs)


snd.moved.location.df = all.locations.df %>%
  filter(!is_test) %>%
  add_prev_locations(prev.loc.set.diff=1) %>% 
  filter(!current_loc, prev_loc_set==2) %>%
  ungroup() %>%
  dplyr::select(-date) %>%
  dplyr::distinct()

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
trial.field.dist2moved.loc = bind_rows(
    dplyr::mutate(trial.field.dist2fst.moved, moved_loc = 'fst'),
    dplyr::mutate(trial.field.dist2snd.moved, moved_loc = 'snd')
  ) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))

beforetest.field.dist2current.rew = beforetest.peaks.at.rew %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id))
beforetest.field.dist2fst.moved = beforetest.peaks.fst.moved
beforetest.field.dist2snd.moved = create.peaks.df.at.locs(beforetest.fields, snd.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2moved.loc = bind_rows(
    mutate(beforetest.field.dist2fst.moved, moved_loc = 'fst'),
    mutate(beforetest.field.dist2snd.moved, moved_loc = 'snd')
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

perc2dist = 1.2
goal.cell.max.dist = 20 / perc2dist

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
                           cell_id=(as.numeric(animal) - 1) * max.ncells + cell_id)
  
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

get.rew.active.animal_cells = function(dist.df, dist.thr = goal.cell.max.dist) {
  dist.df %>%
    filter(!!mindist.var <= dist.thr) %>%
    filter(signif.si) %>%
    dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
    dplyr::pull(animal_cell)
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
cells.transloc.rew.active.l3d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, day_desc == 'learning3 day#1 test', moved_loc == 'snd'))

# Distance to current reward for place cells active at translocated reward
selected.cells.peaks2current.rew.dist = bind_rows(
  dplyr::filter(trial.field.dist2moved.loc, day_desc == 'habituation day#3', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2moved.loc, day_desc == 'learning2 day#1 test', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2current.rew, day_desc %in% c('learning3 day#1 test', 'learning3 day#3 test'))
) %>% 
  dplyr::mutate(active = !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l2d1t)

filtered.cells.peaks2current.rew.dist = filter.cell.present.ntimes(selected.cells.peaks2current.rew.dist, 2, day_desc, !!mindist.var)
filtered.cells.peaks2current.rew.dist$compared_location = 'current reward'

filtered.cells.peaks2fst.rew.dist = bind_rows(
  dplyr::filter(trial.field.dist2moved.loc, day_desc == 'habituation day#3', moved_loc == 'fst'),
  dplyr::filter(beforetest.field.dist2moved.loc, 
                day_desc %in% c('learning2 day#1 test', 'learning3 day#1 test', 'learning3 day#3 test'),
                moved_loc == 'fst')) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active =  !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active.l2d1t) %>%
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
  scale_x_discrete(labels=c('habituation', 'learning1 probe', 'learning2 probe', 'learning3 probe'))

animal.learning2test.rew.following = bind_rows(
    filter(beforetest.field.dist2current.rew, 
           day_desc == 'learning3 day#1 test',
           animal_cell %in% cells.transloc.rew.active.l2d1t),
    filter(beforetest.field.dist2current.rew, 
           day_desc == 'learning3 day#3 test',
           animal_cell %in% cells.transloc.rew.active.l3d1t)
  ) %>%
  mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>% 
  group_by(animal, day_desc, cell_id) %>%
  dplyr::summarise(npresent = sum(!is.na(!!mindist.var)), 
                   nactive = sum(is.at.rew)) %>% 
  filter(npresent == 1) %>% # include only cells found on both days
  left_join(mouse.meta.df, by='animal') %>% 
  group_by(implant, animal, day_desc) %>% 
  dplyr::summarise(at.rew = sum(nactive), at.rew.pct = mean(nactive), n=n()) %>%
  arrange(animal, day_desc)
# TODO include only trials when learnt

beforetest.animal.peaks.at.rew.pct = beforetest.field.dist2current.rew %>%
  filter(signif.si) %>%
  filter(animal != 'A-BL', animal != 'L-TL') %>%
  filter(day_desc %in% c('learning3 day#1 test', 'learning3 day#3 test')) %>%
  dplyr::mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>%
  dplyr::group_by(animal, implant, day_desc) %>%
  dplyr::summarise(at.rew = sum(is.at.rew), at.rew.pct = mean(is.at.rew), n=n()) %>%
  arrange(animal, day_desc)

print('T test on reward cells percent - higher number of cells moves to reward than in general population')
cells.following.pct = bind_rows(
  mutate(animal.learning2test.rew.following, cell_group='reward_cells'), 
  mutate(beforetest.animal.peaks.at.rew.pct, cell_group='all_PCs')) %>%
  reshape2::dcast(implant + animal + day_desc ~ cell_group, value.var='at.rew.pct') %>% 
  filter(!is.na(reward_cells)) %>%
  reshape2::melt(id.vars=c('implant', 'animal', 'day_desc'), value.name='at.rew.pct', variable.name='cell_group') 
  #replace_na(list('at.rew.pct'=0))
wilcox.test(at.rew.pct ~ cell_group,
            data=cells.following.pct,
            subset=implant=='dCA1' & animal != 'C-1R',
            paired=TRUE)

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
