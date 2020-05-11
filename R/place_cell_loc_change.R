library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)

# Plotting
library(ggplot2)
library(cowplot)
source('plotting_params.R')

summarise = dplyr::summarise
summarize = dplyr::summarize


px2pt = 1/(ggplot2::.pt*72.27/96)


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

beforetest.peaks.fst.moved = dplyr::mutate(beforetest.peaks.fst.moved, 
                                           day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))
beforetest.peaks.at.rew = dplyr::mutate(beforetest.peaks.at.rew,
                                        day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))

selected.peaks.df = bind_rows(
  #filter(trial.peaks.at.fst.moved, day_desc=='habituation day#1'),
  filter(trial.peaks.at.fst.moved, day_desc=='habituation day#3'),
  filter(trial.peaks.at.fst.moved, day_desc=='learning1 day#1'),
  filter(trial.peaks.at.fst.moved, day_desc=='learning1 day#5'),
  #filter(beforetest.peaks.fst.moved, day_desc=='learning2 day#1'),
  filter(trial.peaks.at.fst.moved, day_desc=='learning2 day#1'),
  #filter(trial.peaks.at.fst.moved, day_desc=='learning2 day#2'),
  filter(trial.peaks.at.fst.moved, day_desc=='learning3 day#2')
) %>%
  dplyr::mutate(active=as.integer(rew.peaks.count > 0),
                exp=map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::select(animal, day_desc, cell_id, exp, active, peak2rew.mindist)

ndays = selected.peaks.df$exp %>% unique %>% length

peakatrew.clusters = filter.cell.present.ntimes(selected.peaks.df, ndays, exp, active)
peakatrew.clusters$my.clust = as.integer(peakatrew.clusters$active) + 1

peakatrew.clusters = left_join(peakatrew.clusters, mouse.meta.df, by=c('animal'))
plot.cluster.assignments(filter(peakatrew.clusters, 
                                animal=='F-TL'),
                                #implant=='vCA1'),
                         clust.distance = 100,
                         pattern.desc.fun = get.activity.pattern.desc,
                         #hd1, 
                         hd3, 
                         l1d1,
                         l1d5, 
                         #l2d1t, 
                         l2d1, 
                         #l2d2
                         l3d2) +
  labs(title='Reward activated place cells at first translocated location')


#######################################################################
# Plot change of the distance to first translocated reward zone
#######################################################################

peakatrew.mindist = filter.cell.present.ntimes(trial.peaks.at.fst.moved, 6, day_desc, peak2rew.mindist)
peakatrew.mindist = left_join(peakatrew.mindist, mouse.meta.df)
peakatrew.mindist$animal = as.factor(peakatrew.mindist$animal)
peakatrew.mindist %>%
  #filter(implant=='dCA1') %>%
  filter(animal=='D-BR') %>%
  dplyr::mutate(unique_cell_id=(as.numeric(animal) - 1)*1000 + cell_id) %>%
  ggplot(aes(x=day_desc, y=peak2rew.mindist, group=unique_cell_id)) +
  geom_point(alpha=0.1) +
  #geom_smooth(formula=y ~ poly(x, 2), method='lm', se=FALSE, alpha=0.2) +
  geom_path(alpha=0.2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('') + ylab('Zone distance (%)')


#########################################################################################
# Show change in distance to the first translocated reward location, filtered for stable cells
#########################################################################################

day.ordinal.rew.changed = 9
trial.peaks.at.fst.moved$day_desc = as.factor(trial.peaks.at.fst.moved$day_desc)
trial.peaks.at.fst.moved = trial.peaks.at.fst.moved %>%
  dplyr::mutate(day_ordinal = as.numeric(day_desc),
                day_ordinal = day_ordinal + ifelse(day_ordinal >= day.ordinal.rew.changed, 1, 0),
                next_day_ordinal = day_ordinal + 1)

l2d1.beforetest = filter(beforetest.peaks.fst.moved, day_desc=='learning2 day#1') %>%
  dplyr::mutate(day_ordinal=day.ordinal.rew.changed,
                day_desc = 'learning1 day#6 test',
                next_day_ordinal = day_ordinal + 1) %>%
  left_join(mouse.meta.df) %>%
  left_join(dplyr::select(pc.test.df, cell_id, animal, date, signif.si))
joined.test.trial.peaks = bind_rows(trial.peaks.at.fst.moved, l2d1.beforetest)

joined.nextday.peaks =  joined.test.trial.peaks %>%
  left_join(joined.test.trial.peaks, 
            by=c('implant'='implant', 'animal'='animal', 'cell_id'='cell_id', 'next_day_ordinal'='day_ordinal'),
            suffix=c('.fst', '.snd')) %>%
  filter(!is.na(peak2rew.mindist.snd))

# How many times does a cell cross from close to far
cell.activity.changes = joined.nextday.peaks %>%
  dplyr::mutate(changed.rew.activity=rew.peaks.count.fst != rew.peaks.count.snd ) %>%
  dplyr::group_by(implant, animal, cell_id) %>%
  dplyr::summarise(nday.pairs=n(),
                   ndays=unique(c(day_desc.fst, day_desc.snd)) %>% length,
                   nchanged.rew.activity=sum(changed.rew.activity),
                   nactive=unique(c(ifelse(rew.peaks.count.fst > 0, day_desc.fst, '-1'),
                                    ifelse(rew.peaks.count.snd > 0, day_desc.snd, '-1'))) %>% 
                            setdiff('-1') %>% length
                   ) 

stable.animal.cells = cell.activity.changes %>%
  filter(animal=='K-BR') %>%
  filter(nchanged.rew.activity / nday.pairs < 0.3)
  #filter(nday.pairs >= 6) %>%
  #filter(nchanged.rew.activity <= 2) 
  
perc2dist = 1.8

joined.nextday.peaks %>%
  #filter(animal=='K-BR') %>%
  filter(implant=='vCA1') %>%
  dplyr::mutate(stable_cell=(cell_id %in% stable.animal.cells$cell_id)) %>%
  filter(signif.si.fst) %>%
  filter(day_desc.snd != 'learning3 day#3') %>%
  ggplot(aes(x=day_desc.fst, y=peak2rew.mindist.fst * perc2dist)) +
  geom_segment(aes(xend=day_desc.snd, yend=peak2rew.mindist.snd * perc2dist),#, color=stable_cell),
               alpha=0.2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('') + ylab('Distance (cm)')

# TODO: show if the cells further away are in proximity of the other reward?


#########################################################################################
# Plot change of cell activity at translocated reward, median by the learning stage
#########################################################################################

# Group days by stages and show active change schematic on groups
rew.dist.threshold = 20
min.nstages = 4

stage.peaks.df = joined.test.trial.peaks %>%
  dplyr::mutate(exp=map_chr(day_desc, get.stage.from.day.desc)) %>%
  group_by(animal, cell_id, exp) %>%
  dplyr::summarise(med.rew.dist=median(peak2rew.mindist, na.rm=TRUE),
                   ndays.stage=n()) %>%
  dplyr::mutate(active = as.integer(med.rew.dist <= rew.dist.threshold))
stage.peaks.filtered = filter.cell.present.ntimes(stage.peaks.df, min.nstages, exp, active)

stage.peaks.filtered$my.clust = as.integer(stage.peaks.filtered$active) + 1
stage.peaks.filtered = left_join(stage.peaks.filtered, mouse.meta.df)

get.stage.activity.pattern.desc = function(clust_char) {
  if (stringr::str_detect(clust_char, '^\\d2222$')) {
    return('1-place_cell')
  }
  if (stringr::str_detect(clust_char, '^\\d22\\d\\d$')) {
    return('2-reward_memory_cell')
  }
  return('other')
}

plot.cluster.assignments(filter(stage.peaks.filtered, 
                                #cell_id %in% stable.animal.cells$cell_id,
                                animal=='K-BR'),
                                #implant=='vCA1'),
                         clust.distance = 100,
                         pattern.desc.fun = get.stage.activity.pattern.desc,
                         H1,
                         L1,
                         #L1_test,
                         L2,
                         L3
) +
  scale_x_continuous(labels=c('habituation','learning 1', 'learning 2', 'learning 3')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title='Reward activated place cells at first translocated location\n(stage median)')

########################################################################################################
# For the cells active at the translocated reward, were they active at reward before or after?
# View of how the peaks in the place field of reward-active cell followed the reward.

join.chars = function(animal, cell_id) {
  stringr::str_c(animal, cell_id, sep='_')
}

cells.transloc.rew.active = trial.peaks.at.fst.moved %>%
  filter(day_desc == 'learning1 day#5', rew.peaks.count > 0, signif.si) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  dplyr::pull(animal_cell)

#############################################################################
# Distance to current reward for place  (cells active at translocated reward)
#############################################################################
filtered.cells.peaks.at.moved.rew = bind_rows(
  dplyr::filter(trial.peaks.at.current.rew, day_desc %in% c('learning2 day#2', 'learning3 day#2')),
  dplyr::filter(trial.peaks.at.fst.moved, day_desc %in% c('learning1 day#5', 'habituation day#3')),
  dplyr::filter(beforetest.peaks.fst.moved, day_desc %in% c('learning2 day#1 test')),
  dplyr::filter(beforetest.peaks.at.rew, day_desc %in% c('learning3 day#1 test', 'learning3 day#3 test'))
) %>% 
  dplyr::mutate(exp = map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active = rew.peaks.count > 0) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active)

peakatrew.mindist = filter.cell.present.ntimes(filtered.cells.peaks.at.moved.rew, 2, exp, peak2rew.mindist)
peakatrew.mindist$compared_location = 'current reward'

peakat.prev.rew.mindist = bind_rows(
  dplyr::filter(trial.peaks.at.fst.moved, day_desc %in% c('habituation day#3',
                                                          'learning1 day#5',
                                                          'learning2 day#2',
                                                          'learning3 day#2'
                                                          )),
  dplyr::filter(beforetest.peaks.fst.moved, day_desc %in% c('learning2 day#1 test',
                                                            'learning3 day#1 test',
                                                            'learning3 day#3 test'))) %>%
  dplyr::mutate(exp = map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active = rew.peaks.count > 0) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active) %>%
  filter.cell.present.ntimes(2, exp, peak2rew.mindist)
peakat.prev.rew.mindist$compared_location = 'translocated reward'
  
left_join(bind_rows(peakatrew.mindist, peakat.prev.rew.mindist), mouse.meta.df) %>%
  #filter(animal=='D-BR') %>%
  dplyr::mutate(unique_cell_id=join.chars(animal, cell_id)) %>%
  ggplot(aes(x=exp, y=peak2rew.mindist * perc2dist)) +
  geom_path(mapping=aes(group=unique_cell_id), alpha=0.2) +
  geom_violin(alpha=0.6) +
  stat_summary(fun.y=median, geom="point", size=2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(implant ~ compared_location) +
  #facet_wrap(implant + animal ~ compared_location, ncol=4) +
  geom_hline(yintercept = rew.dist.threshold * perc2dist, linetype = 'dashed') +
  xlab('') + ylab('Distance to reward (cm)') +
  scale_x_discrete(labels=c('habituation','learning1 day5','learning1 probe', 'learning 2', 'learning2 probe', 'learning 3', 'learning3 probe'))

###
# Only probe tests
#####
cells.transloc.rew.active = beforetest.peaks.fst.moved %>%
  filter(day_desc == 'learning2 day#1 test', rew.peaks.count > 0, signif.si) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  dplyr::pull(animal_cell)

# Distance to current reward for place  (cells active at translocated reward)
probe.peaks.at.moved.rew = bind_rows(
  dplyr::filter(beforetest.peaks.fst.moved, day_desc %in% c('learning2 day#1 test')),
  dplyr::filter(beforetest.peaks.at.rew, day_desc %in% c('learning3 day#1 test', 'learning3 day#3 test'))
)  %>% 
  dplyr::mutate(exp = map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active = rew.peaks.count > 0) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active)

filtered.cells.peaks.at.moved.rew = bind_rows(
  dplyr::filter(trial.peaks.at.fst.moved, day_desc %in% c(#'learning1 day#5', 
                                                          'habituation day#3')),
  probe.peaks.at.moved.rew
) %>% 
  dplyr::mutate(exp = map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active = rew.peaks.count > 0) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active)

peakatrew.mindist = filter.cell.present.ntimes(filtered.cells.peaks.at.moved.rew, 2, exp, peak2rew.mindist)
peakatrew.mindist$compared_location = 'current reward'

peakat.prev.rew.mindist = bind_rows(
  dplyr::filter(trial.peaks.at.fst.moved, day_desc %in% c(#'learning1 day#5',
                                                          'habituation day#3')),
  dplyr::filter(beforetest.peaks.fst.moved, day_desc %in% c('learning2 day#1 test',
                                                            'learning3 day#1 test',
                                                            'learning3 day#3 test'))) %>%
  dplyr::mutate(exp = map_chr(day_desc, get.shortname.from.day.desc)) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id),
                active = rew.peaks.count > 0) %>%
  dplyr::filter(animal_cell %in% cells.transloc.rew.active) %>%
  filter.cell.present.ntimes(2, exp, peak2rew.mindist)
peakat.prev.rew.mindist$compared_location = 'translocated reward'

left_join(bind_rows(peakatrew.mindist, peakat.prev.rew.mindist), mouse.meta.df) %>%
  #filter(animal=='D-BR') %>%
  dplyr::mutate(unique_cell_id=join.chars(animal, cell_id)) %>%
  ggplot(aes(x=exp, y=peak2rew.mindist * perc2dist)) +
  geom_path(mapping=aes(group=unique_cell_id), alpha=0.2) +
  geom_violin(alpha=0.6) +
  stat_summary(fun.y=median, geom="point", size=2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #facet_wrap(implant + animal ~ compared_location, ncol=4) +
  facet_grid(implant ~ compared_location) +
  geom_hline(yintercept = rew.dist.threshold * perc2dist, linetype = 'dashed') +
  xlab('') + ylab('Distance to reward (cm)') +
  scale_x_discrete(labels=c('habituation', 'learning1 probe', 'learning2 probe', 'learning3 probe'))




################
# Cluster plot
#################

get.reward.pattern.desc = function(clust_char) {
  if (stringr::str_detect(clust_char, '^\\d211$')) {
    return('2-place_cell')
  }
  if (stringr::str_detect(clust_char, '^\\d222$')) {
    return('1-reward_cell')
  }
  return('other')
}

filtered.cells.peaks.at.moved.rew.active = filter.cell.present.ntimes(filtered.cells.peaks.at.moved.rew, 4, exp, active)
filtered.cells.peaks.at.moved.rew.active$my.clust = as.integer(filtered.cells.peaks.at.moved.rew.active$active) + 1
filtered.cells.peaks.at.moved.rew.active = left_join(filtered.cells.peaks.at.moved.rew.active, mouse.meta.df, by=c('animal'))
plot.cluster.assignments(filter(filtered.cells.peaks.at.moved.rew.active, 
                                #animal=='F-TL'),
                                implant=='dCA1'),
                         clust.distance = 100,
                         pattern.desc.fun = get.reward.pattern.desc,
                         #hd1, 
                         hd3, 
                         #l1d1,
                         #l1d5, 
                         l2d1t, 
                         #l2d1,
                         l3d1t, 
                         #l2d2
                         l3d3t) +
  labs(title='Changes in translocated active cells location')

# TODO: 
# - plot with minimum distance to reward
# - plot using the translocated location before and current locations after

#TODOs regenerate the peaks and present the figures..
# name the groups: group_desc
# place_cell, reward cell, reward silenced cell, reward_move activated
# quantify 
# Should use only cells with signif spatial information in the quantification?
