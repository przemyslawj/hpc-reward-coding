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
source('fixed_effects.R')

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
  filter(!(loc_id %in% thd.moved.location.df$loc_id)) 
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

add.beforetest.ordinal = function(df, day.var=day_desc, ordinal.var=loc_set_ordinal) {
  day.var = enquo(day.var)
  ordinal.var = enquo(ordinal.var)
  group_vars = quos(animal, !!day.var)
  loc.ordinal.info = df %>%
    dplyr::select(!!!group_vars) %>%
    dplyr::distinct() %>%
    dplyr::arrange(!!!group_vars) %>%
    dplyr::group_by(animal) %>%
    dplyr::mutate(!!ordinal.var:=dplyr::row_number())
  left_join(df, loc.ordinal.info, by=c('animal', quo_name(day.var)))
}

trial.field.dist2current.rew = trial.peaks.at.current.rew %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  add.beforetest.ordinal()

trial.field.dist2fst.moved = trial.peaks.at.fst.moved
trial.field.dist2snd.moved = create.peaks.df.at.locs(run.fields, snd.moved.location.df, filter(trial.days.df, exp=='habituation'))
trial.field.dist2thd.moved = create.peaks.df.at.locs(run.fields, thd.moved.location.df, filter(trial.days.df, exp=='habituation'))
trial.field.dist2moved.loc = bind_rows(
    dplyr::mutate(trial.field.dist2fst.moved, moved_loc = 'fst'),
    dplyr::mutate(trial.field.dist2snd.moved, moved_loc = 'snd'),
    dplyr::mutate(trial.field.dist2thd.moved, moved_loc = 'thd'),
  ) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  add.beforetest.ordinal()

beforetest.field.dist2current.rew = beforetest.peaks.at.current.rew %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  add.beforetest.ordinal()

beforetest.field.dist2fst.moved = beforetest.peaks.fst.moved
beforetest.field.dist2snd.moved = create.peaks.df.at.locs(beforetest.fields, snd.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2thd.moved = create.peaks.df.at.locs(beforetest.fields, thd.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2fourth.moved = create.peaks.df.at.locs(beforetest.fields, fourth.moved.location.df, trial.days.df, cell.db = beforetest.cell.db)
beforetest.field.dist2moved.loc = bind_rows(
    mutate(beforetest.field.dist2fst.moved, moved_loc = 'fst'),
    mutate(beforetest.field.dist2snd.moved, moved_loc = 'snd'),
    mutate(beforetest.field.dist2thd.moved, moved_loc = 'thd'),
    mutate(beforetest.field.dist2fourth.moved, moved_loc = 'fourth')
  ) %>%
  dplyr::mutate(animal_cell = join.chars(animal, cell_id)) %>%
  add.beforetest.ordinal()

mindist.var = quo(min.rew.dist)

beforetest.field.dist2moved.loc = dplyr::mutate(beforetest.field.dist2moved.loc,
                                                day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))
beforetest.field.dist2current.rew = dplyr::mutate(beforetest.field.dist2current.rew,
                                                  day_desc = ifelse(endsWith(day_desc, 'test'), day_desc, paste(day_desc, 'test')))




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

#############################################################################################
# Distance to current reward for place cells active at fst translocated reward - test probes
#############################################################################################
cells.transloc.rew.active.l2d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, loc_set_ordinal==1, moved_loc=='fst'))

# Baseline cells: inactive at the current rewards and not present at the subsequent reward locations
cells.transloc.rew.inactive.l2d1t = get.animal_cells(
  filter(beforetest.field.dist2moved.loc, 
         loc_set_ordinal==1,
         moved_loc %in% c('fst', 'snd')) %>%
    dplyr::group_by(animal, cell_id) %>%
    dplyr::summarise(!!mindist.var := min(!!mindist.var), signif.si=any(signif.si)),
  expr(!!mindist.var > goal.cell.max.dist))

cells.transloc.rew.active.l3d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, loc_set_ordinal==2, moved_loc=='snd'))
cells.transloc.rew.inactive.l3d1t = get.animal_cells(
  filter(beforetest.field.dist2moved.loc, 
         loc_set_ordinal==2,
         moved_loc %in% c('snd', 'thd')) %>%
    dplyr::group_by(animal, cell_id) %>%
    dplyr::summarise(!!mindist.var := min(!!mindist.var), signif.si=any(signif.si)),
  expr(!!mindist.var > goal.cell.max.dist))

cells.transloc.rew.active.l4d1t = get.rew.active.animal_cells(
  filter(beforetest.field.dist2moved.loc, loc_set_ordinal==3, moved_loc=='thd'))
cells.transloc.rew.inactive.l4d1t = get.animal_cells(
  filter(beforetest.field.dist2moved.loc, 
         loc_set_ordinal==3,
         moved_loc %in% c('thd', 'fourth')) %>%
    dplyr::group_by(animal, cell_id) %>%
    dplyr::summarise(!!mindist.var := min(!!mindist.var), signif.si=any(signif.si)),
  expr(!!mindist.var > goal.cell.max.dist))


pc.beforetest.change.df = dplyr::mutate(pc.beforetest.change.df, 
                                        animal_cell = join.chars(animal, cell_id))
remap.corr.thr = 0.3
remapped.cells.1to2 = filter(pc.beforetest.change.df,
                             loc_set_ordinal.fst==1,
                             loc_set_ordinal.snd==2,
                             n > 1,
                             field.cor < remap.corr.thr) %>%
  dplyr::pull(animal_cell)
remapped.cells.2to3 = filter(pc.beforetest.change.df,
                             loc_set_ordinal.fst==2,
                             loc_set_ordinal.snd==3,
                             n > 1,
                             field.cor < remap.corr.thr) %>%
  dplyr::pull(animal_cell)
remapped.cells.3to4 = filter(pc.beforetest.change.df,
                             loc_set_ordinal.fst==3,
                             loc_set_ordinal.snd==4,
                             n > 1,
                             field.cor < remap.corr.thr) %>%
  dplyr::pull(animal_cell)
remapped.cells.fromL1toL3 = filter(pc.beforetest.change.df,
                             loc_set_ordinal.fst==1,
                             loc_set_ordinal.snd==3,
                             n > 1,
                             field.cor < remap.corr.thr) %>%
  dplyr::pull(animal_cell)


# DF with distance to the current reward, filtered to cells that 
# were reward-active during the previous beforetest
dist2current.rew.active = bind_rows(
  filter(beforetest.field.dist2current.rew, 
         loc_set_ordinal==2,
         animal_cell %in% cells.transloc.rew.active.l2d1t) %>%
      dplyr::mutate(remapped = animal_cell %in% remapped.cells.1to2),
  filter(beforetest.field.dist2current.rew, 
         loc_set_ordinal==3,
         animal_cell %in% cells.transloc.rew.active.l3d1t) %>%
      dplyr::mutate(remapped = animal_cell %in% remapped.cells.2to3),
  filter(beforetest.field.dist2current.rew, 
         loc_set_ordinal==4,
         animal_cell %in% cells.transloc.rew.active.l4d1t) %>%
    dplyr::mutate(remapped = animal_cell %in% remapped.cells.3to4))

dist2prev.rew.active = bind_rows(
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='fst',
         loc_set_ordinal==2,
         animal_cell %in% cells.transloc.rew.active.l2d1t),
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='snd',
         loc_set_ordinal==3,
         animal_cell %in% cells.transloc.rew.active.l3d1t),
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='thd',
         loc_set_ordinal==4,
         animal_cell %in% cells.transloc.rew.active.l4d1t) 
)

dist2prev.rew.inactive = bind_rows(
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='fst',
         loc_set_ordinal==2,
         animal_cell %in% cells.transloc.rew.inactive.l2d1t),
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='snd',
         loc_set_ordinal==3,
         animal_cell %in% cells.transloc.rew.inactive.l3d1t),
  filter(beforetest.field.dist2moved.loc, 
         moved_loc=='thd',
         loc_set_ordinal==4,
         animal_cell %in% cells.transloc.rew.inactive.l4d1t) 
)


dist2current.rew.active.trials = bind_rows(
  filter(trial.field.dist2current.rew, 
         day_desc=='learning2 day#1',
         animal_cell %in% cells.transloc.rew.active.l2d1t),
  filter(trial.field.dist2current.rew, 
         day_desc=='learning3 day#1',
         animal_cell %in% cells.transloc.rew.active.l3d1t),
  filter(trial.field.dist2current.rew, 
         day_desc=='learning4 day#1',
         animal_cell %in% cells.transloc.rew.active.l4d1t))


# Summary DF with distance to the current reward filtered to the cells that remapped
# and were reward-active during the previous beforetest 
animal.learning2test.rew.following = dist2current.rew.active %>%
  #filter(remapped) %>%
  dplyr::mutate(is.at.rew = !!mindist.var <= goal.cell.max.dist) %>% 
  group_by(implant, animal, day_desc) %>% 
  dplyr::summarise(at.rew = sum(is.at.rew), 
                   at.rew.pct = mean(is.at.rew) * 100, 
                   quantile30.dist = unname(quantile(!!mindist.var, 0.3)[1]),
                   med.dist = median(!!mindist.var), 
                   mean.dist = mean(!!mindist.var),
                   n=n(), .groups='drop') %>%
  arrange(animal, day_desc)

# DF with distance to the current reward, filtered to cells that
# # were reward-inactive during the previous beforetest
dist2current.rew.inactive = bind_rows(
  filter(beforetest.field.dist2current.rew,
         loc_set_ordinal==2,
         animal_cell %in% cells.transloc.rew.inactive.l2d1t) %>%
    dplyr::mutate(remapped=animal_cell %in% remapped.cells.1to2),
  filter(beforetest.field.dist2current.rew,
         loc_set_ordinal==3,
         animal_cell %in% cells.transloc.rew.inactive.l3d1t) %>%
    dplyr::mutate(remapped=animal_cell %in% remapped.cells.2to3),
  filter(beforetest.field.dist2current.rew,
         loc_set_ordinal==4,
         animal_cell %in% cells.transloc.rew.inactive.l4d1t) %>%
    dplyr::mutate(remapped=animal_cell %in% remapped.cells.3to4))

dist2current.rew.inactive.trials = bind_rows(
  filter(trial.field.dist2current.rew,
         day_desc == 'learning2 day#1',
         animal_cell %in% cells.transloc.rew.inactive.l2d1t),
  filter(trial.field.dist2current.rew,
         day_desc == 'learning3 day#1',
         animal_cell %in% cells.transloc.rew.inactive.l3d1t),
  filter(trial.field.dist2current.rew,
         day_desc == 'learning4 day#1',
         animal_cell %in% cells.transloc.rew.inactive.l4d1t))

dist2current.rew.df = bind_rows(
  dplyr::mutate(dist2current.rew.active, group='active'),
  dplyr::mutate(dist2current.rew.inactive, group='inactive'))

dist2current.rew.trials.df = bind_rows(
  dplyr::mutate(dist2current.rew.active.trials, group='active'),
  dplyr::mutate(dist2current.rew.inactive.trials, group='inactive'))

trialday.closer2rew = filter(trial.field.dist2current.rew, day_ordinal==2 | day_ordinal == 1) %>%
  filter(npeaks > 0) %>%
  group_by(implant, animal, cell_id, location_set) %>%
  dplyr::summarise(min.rew.dist.trial = min(min.rew.dist),
                   signif.si = signif.si[which.min(min.rew.dist)])

test.trial.activity.df = left_join(
  dist2current.rew.active,
  trialday.closer2rew
    %>% mutate(compared_loc_ordinal=location_set),
  suffix=c('.test', '.trial'),
  by=c('implant'='implant', 'animal'='animal', 'cell_id'='cell_id', 
          'location_set'='compared_loc_ordinal'))  %>%
  rename(min.rew.dist.test = min.rew.dist) %>%
  filter(!is.na(min.rew.dist.trial))

test.trial.activity.df %>%
  group_by(implant, min.rew.dist.test <= goal.cell.max.dist) %>%
  dplyr::summarise(mean.dist = mean(min.rew.dist.trial * perc2dist), 
                   rew.cells = sum(min.rew.dist.trial <= goal.cell.max.dist) / n(),
                   med.dist = median(min.rew.dist.trial * perc2dist), n())


eval.val.current.rew = function(dist2current.rew.df, 
                                implant_loc = 'vCA1', 
                                val.expr=log(min.rew.dist * perc2dist),
                                yval.inverse=exp) {
  val.expr = enexpr(val.expr)
  dist2current.rew.df = mutate(dist2current.rew.df, val=!!val.expr)
  dist2current.rew.df$implant = as.factor(dist2current.rew.df$implant)
  dist2current.rew.df$group = as.factor(dist2current.rew.df$group)
  rew.dist.model = lmerTest::lmer(val ~ group + (1 | animal),
                                  dist2current.rew.df,
                                  REML=TRUE,
                                  subset=implant==implant_loc)
  g = plot.model.diagnostics(rew.dist.model,
                             subset(dist2current.rew.df, implant==implant_loc)$animal,
                             subset(dist2current.rew.df, implant==implant_loc)$group)
  print(g)
  print(summary(rew.dist.model))
  print(anova(rew.dist.model, refit=FALSE, ddf='Satterthwaite'))
  
  models = create.bayes.lm.pair(filter(dist2current.rew.df, implant==implant_loc),
                                val ~ 1 + group + animal,
                                val ~ 1 + animal,
                                whichRandom = 'animal',
                                iterations = 100000)
  print(models$full / models$null)
  calc.pair.95CI(models$full, show.percent.change = FALSE,
                 ytransform = yval.inverse,
                 pair.vars = c('group-active', 'group-inactive'))
}
eval.val.current.rew(dist2current.rew.df, 'dCA1')
eval.val.current.rew(dist2current.rew.df, 'vCA1')

# Field size not different in vCA1 between active and not-active group
eval.val.current.rew(dist2current.rew.df, 'vCA1', field.size.50 / nbins / nbins * 100, identity)

# How many cells remapped in each group
dist2current.rew.df %>%
  group_by(implant, group) %>%
  dplyr::summarise(persistent_field=1-mean(remapped),
                   sum(min.rew.dist <= goal.cell.max.dist),
                   ncells=n())
# 
# dist2current.rew.df %>%
#   group_by(implant, group) %>%
#   ggplot() +
#   geom_density(aes(x=min.rew.dist * perc2dist, group=group, color=group)) +
#   facet_wrap(. ~ implant) +
#   scale_color_manual(values=c(side.two.coulours[1], '#999999')) +
#   xlab('Distance to reward (cm)') + 
#   ylab('Cell density') + 
#   labs(linetype='', title='Field distance to reward after translocation') +
#   xlim(c(0,100))+
#   gtheme

dist2current.rew.df %>%
  group_by(implant, group) %>%
  ggplot(aes(x=min.rew.dist * perc2dist, group=group, color=group)) +
  stat_ecdf(geom='step') +
  facet_grid(implant ~ .) +
  scale_color_manual(values=c(side.two.coulours[1], '#999999')) +
  xlab('Distance to reward (cm)') + 
  ylab('Cell CDF') + 
  xlim(c(0,100)) +
  scale_y_continuous(breaks=c(0, 0.5, 1.0)) + 
  gtheme + theme(legend.position = 'none')
figure.ggsave('Figs4B-reward_active_dist_to_current_reward.pdf', 
              height=5.4, width=4.1)

dist2current.rew.df %>%
  group_by(implant, group, animal) %>%
  mutate(line_id=paste(group,animal)) %>%
  ggplot(aes(x=min.rew.dist * perc2dist, group=line_id, color=group)) +
  stat_ecdf(geom='step') +
  facet_wrap(implant ~ animal) +
  scale_color_manual(values=c(side.two.coulours[1], '#999999')) +
  xlab('Distance to reward (cm)') + 
  ylab('Cell CDF') + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0)) + 
  scale_x_continuous(breaks=c(0, 50, 100), limits = c(0,100)) + 
  gtheme + theme(legend.position = 'none')

figure.ggsave('Figs4B-reward_active_dist_to_current_reward-per-animal.pdf', 
       height=9.8, width=7.2)
