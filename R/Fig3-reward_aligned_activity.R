source('reward_approach_activity_utils.R')

learning_result_dirs = Filter(function(x) {str_detect(x, 'learning')}, caimg_result_dirs)

mouse.meta.df = read.mouse.meta(rootdirs)
trials.meta.df = read.trials.meta(rootdirs)


reward.aligned.pop.traces.list = list()
reward.aligned.pop.traces.by.pc.list = list()
beforetest.aligned.pop.traces.list = list()
beforetest.aligned.pop.traces.by.pc.list = list()
nonreward.aligned.pop.traces.list = list()
immobility.aligned.pop.traces.list = list()
all.rew.arriving.bouts.list = list()
all.cells.fr.outside.seq.list = list()
all.cors.df.list = list()
place.cell.db.min = as.data.table(place.cell.db)[,.(animal, date, cell_id, is.pc, spatial.information)]
beforetest.cell.db.min = as.data.table(beforetest.cell.db)[, .(animal, date, cell_id, is.pc, spatial.information)]
  
i = 0
for (caimg_result_dir in learning_result_dirs) {
  i = i + 1
  binned.traces = prepare.binned.traces(caimg_result_dir, timebin.dur.msec=timebin.dur.msec)
  behav.traces = get.behav.traces(binned.traces)

  # Running bouts
  mobility.bouts = get.mobility.bouts(behav.traces, timebin.dur.msec=timebin.dur.msec,
                                      mobility.expr=is_running)
  reward.arriving.bouts = mobility.bouts[rew_at_end > 0,]
  all.rew.arriving.bouts.list[[i]] = reward.arriving.bouts

  day.reward.aligned.pop.traces = get.event.aligned.pop.traces(binned.traces, reward.arriving.bouts)
  trial.traces.pctype = binned.traces[exp_title == 'trial'][place.cell.db.min, on=.(animal, date, cell_id), nomatch=0]
  day.reward.aligned.pop.traces.pc = get.event.aligned.pop.traces(trial.traces.pctype[(is.pc)], reward.arriving.bouts)[, is.pc:=TRUE]
  day.reward.aligned.pop.traces.nonpc = get.event.aligned.pop.traces(trial.traces.pctype[(!is.pc)], reward.arriving.bouts)[, is.pc:=FALSE]
  reward.aligned.pop.traces.by.pc.list[[i]] = rbindlist(list(day.reward.aligned.pop.traces.pc,
                                                             day.reward.aligned.pop.traces.nonpc))

  rew.aligned.traces = binned.traces[exp_title=='trial',
                                     center.timestamps.around.events(.SD, reward.arriving.bouts, exp_title[1],
                                                                     padding.after.event=2000/timebin.dur.msec),
                                     by=.(exp_title, trial_id, cell_id)][aligned_event_id >= 0,]

  # Find max activity timestamp per cell
  rew.aligned.traces = rew.aligned.traces[timestamp_from_end <= 0 & timestamp_from_start >= 0,
    `:=` (timestamp_to_max = timestamp_from_start - timestamp_from_start[which.max(smoothed_deconv_trace)],
          zscored_smooth_deconv_trace.peak = max(zscored_smooth_deconv_trace, na.rm=TRUE)),
    by=.(exp_title, trial_id, cell_id, aligned_event_id)]
  rew.aligned.traces[timestamp_to_max >= -6000 & timestamp_to_max <= 6000,
                     .(deconv_trace.mean = mean(smoothed_deconv_trace),
                       deconv_trace.sem = sem(smoothed_deconv_trace),
                       timestamp_from_end = timestamp_from_end[1],
                       zscored_smooth_deconv_trace.peak = mean(zscored_smooth_deconv_trace.peak, na.rm=TRUE)),
                     by=.(date, animal, exp_title, cell_id, timestamp_to_max)] -> day.seq.centered.activity
  day.seq.centered.activity[, total.activity := sum(deconv_trace.mean),
                            by=.(date, animal, exp_title, cell_id)]
  day.seq.centered.activity[, .(timestamp.var = sum((timestamp_to_max/1000)^2 * deconv_trace.mean / total.activity),
                                zscored_smooth_deconv_trace.peak = mean(zscored_smooth_deconv_trace.peak),
                                peak.timestamp_from_end = timestamp_from_end[which(timestamp_to_max == 0)]),
                            by=.(date, animal, exp_title, cell_id)] -> day.seq.cell.timestmap.var
  #seq.cell.timestamp.var.list[[i]] = day.seq.cell.timestmap.var

  if (plot.sequences) {
    gs = rew.aligned.traces[, plot.approaches(.SD, trial_id), by=.(exp_title, trial_id)]
  }

  reward.aligned.pop.traces.list[[i]] = day.reward.aligned.pop.traces[aligned_event_id >= 0,]

  nonreward.arriving.bouts = mobility.bouts[rew_at_end == 0 & rew_dist_at_end >= 20,]
  day.nonreward.aligned.pop.traces = get.event.aligned.pop.traces(binned.traces, nonreward.arriving.bouts)
  nonreward.aligned.pop.traces.list[[i]] = day.nonreward.aligned.pop.traces[aligned_event_id >= 0,]

  # Immobility periods at reward
  immobility.bouts = get.mobility.bouts(behav.traces, timebin.dur.msec=timebin.dur.msec,
                                        mobility.expr=!is_running, min.bout.len.msec=8000)
  reward.immobility.bouts = immobility.bouts[rew_at_end > 0,]
  day.immobility.aligned.pop.traces = get.event.aligned.pop.traces(binned.traces, reward.immobility.bouts)
  immobility.aligned.pop.traces.list[[i]] = day.immobility.aligned.pop.traces[aligned_event_id >= 0,]

  # Activity on approaches to reward
  beforetest.behav.traces = behav.traces[exp_title == 'beforetest']
  if (nrow(beforetest.behav.traces) > 0) {
    beforetest.reward.approaches = get.reward.approaches(beforetest.behav.traces,
                                                         timebin.dur.msec=timebin.dur.msec,
                                                         min.bout.len.msec=4000)[rew_dist_at_end < 15]
    beforetest.traces = binned.traces[exp_title == 'beforetest']
    day.beforetest.aligned.pop.traces = get.event.aligned.pop.traces(beforetest.traces,
                                                                     beforetest.reward.approaches)
    #day.beforetest.aligned.pop.traces[, dist_approached_rew := ifelse(rew_at_end==1, dist_reward0, dist_reward1)]
    beforetest.aligned.pop.traces.list[[i]] = day.beforetest.aligned.pop.traces

    beforetest.traces.pctype = beforetest.traces[beforetest.cell.db.min, on=.(animal, date, cell_id), nomatch=0]
    day.beforetest.aligned.pop.traces.pc = get.event.aligned.pop.traces(beforetest.traces.pctype[(is.pc)],
                                                                        beforetest.reward.approaches)[, is.pc:=TRUE]
    day.beforetest.aligned.pop.traces.nonpc = get.event.aligned.pop.traces(beforetest.traces.pctype[(!is.pc)],
                                                                           beforetest.reward.approaches)[, is.pc:=FALSE]
    beforetest.aligned.pop.traces.by.pc.list[[i]] = rbindlist(list(day.beforetest.aligned.pop.traces.pc,
                                                               day.beforetest.aligned.pop.traces.nonpc))

    if (plot.sequences) {
      beforetest.aligned.traces = binned.traces[exp_title=='beforetest',
                                                center.timestamps.around.events(.SD, beforetest.reward.approaches, exp_title[1],
                                                                                padding.after.event=2000/timebin.dur.msec),
                                                by=.(exp_title, trial_id, cell_id)][aligned_event_id >= 0,]
      gs = plot.approaches(beforetest.aligned.traces, beforetest.aligned.traces$trial_id[1])

      #beforetest.aligned.traces[, dist_approached_rew := ifelse(rew_at_end==1, dist_reward0, dist_reward1)]
      dist.bin.width = 10
      dist.traces = stimbin.traces(beforetest.aligned.traces,
                                   dist_approached_rew, 100 / dist.bin.width, max.width = 100)
      binned.dist.traces = dist.traces[timestamp_from_start >= 0 & timestamp_from_end <= 3000,
                                       .(zscored_smooth_deconv_trace = mean(zscored_smooth_deconv_trace),
                                         zscored_deconv_trace = mean(zscored_deconv_trace),
                                         zscored_trace = mean(zscored_trace),
                                         dist_approached_rew = dist.bin.width/2 + (bin.dist_approached_rew - 1) * dist.bin.width),
                                       by=.(animal, date, trial, aligned_event_id, cell_id, bin.dist_approached_rew)]
      gs = plot.approaches.dist(binned.dist.traces, beforetest.aligned.traces$trial_id[1], dist.bin.width=dist.bin.width)
    }
  }

}


join.meta.df = function(df) {
  df %>%
    left_join(mouse.meta.df, by='animal') %>%
    left_join(trials.meta.df, by=c('date', 'animal')) %>%
    as.data.table()
}

reward.aligned.pop.traces = join.meta.df(rbindlist(reward.aligned.pop.traces.list))
reward.aligned.pop.traces.by.pc = join.meta.df(rbindlist(reward.aligned.pop.traces.by.pc.list))
all.rew.arriving.bouts = join.meta.df(rbindlist(all.rew.arriving.bouts.list))
nonreward.aligned.pop.traces = join.meta.df(rbindlist(nonreward.aligned.pop.traces.list))
immobility.aligned.pop.traces = join.meta.df(rbindlist(immobility.aligned.pop.traces.list))
beforetest.aligned.pop.traces = join.meta.df(rbindlist(beforetest.aligned.pop.traces.list))
beforetest.aligned.pop.traces.by.pc = join.meta.df(rbindlist(beforetest.aligned.pop.traces.by.pc.list))
#   seq.cell.timestamp.var = join.meta.df(rbindlist(seq.cell.timestamp.var.list))

save.image(file="data/rew_aligned_pop_activity12-200ms.RData")

calc.acceleration = function(velocity, dt.ms=timebin.dur.msec) {
  dt.bins = ceiling(timebin.dur.msec / dt.ms)
  # Skip the first reading of velocity which often is wrongly high
  velocity = velocity[2:length(velocity)]
  velocity.indexes1 = seq(length(velocity) - dt.bins)
  velocity.indexes2 = velocity.indexes1 + dt.bins
  acc.vec = (velocity[velocity.indexes2] - velocity[velocity.indexes1]) * (1000/dt.ms)
  window.length = max(3, dt.bins * 2- 1)
  smooth.acc.vec = acc.vec
  # Prepend the first acceleration value
  c(0, rep(smooth.acc.vec[1], length(velocity) - length(smooth.acc.vec)), smooth.acc.vec)
}

reward.aligned.pop.traces = reward.aligned.pop.traces %>% 
  arrange(exp_title, trial_id, animal, date, trial, implant, event_seq_no, day_ordinal, location_set, timestamp) %>%
  group_by(exp_title, trial_id, animal, date, trial, implant, event_seq_no, day_ordinal, location_set) %>%
  dplyr::mutate(acceleration = calc.acceleration(velocity, timebin.dur.msec)) %>%
  ungroup()

reward.aligned.pop.traces.by.pc = reward.aligned.pop.traces.by.pc %>% 
  arrange(exp_title, trial_id, animal, date, trial, implant, event_seq_no, day_ordinal, location_set, is.pc, timestamp) %>%
  group_by(exp_title, trial_id, animal, date, trial, implant, event_seq_no, day_ordinal, location_set, is.pc) %>%
  dplyr::mutate(acceleration = calc.acceleration(velocity, timebin.dur.msec)) %>%
  ungroup()

early.learning.days = c('learning1 day#1')
late.learning.days = c('learning1 day#5', 'learning2 day#2', 'learning3 day#2', 'learning4 day#2')

add.early.late.col = function(df, filter.early.and.late=FALSE) {
  df = df %>%
    dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
                  is.late.learning = day_desc %in% late.learning.days) %>%
    dplyr::mutate(bout=ifelse(is.early.learning, 'early', ifelse(is.late.learning, 'late', NA)))
  if (filter.early.and.late) {
    df = filter(df, (is.early.learning | is.late.learning))
  }
  df
}

add.proximal.distal.col = function(df, filter.proximal.distal=TRUE) {
  df = df %>%
    dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
    dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width))
  if (filter.proximal.distal) {
    df = filter(df, distal.timestamp | proximal.timestamp)
  }
  df
}

reward.aligned.pop.traces.earlylate = reward.aligned.pop.traces %>%
  filter(exp_title == 'trial', aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  add.early.late.col(filter.early.and.late = TRUE)

reward.aligned.pop.traces.earlylate.by.pc = reward.aligned.pop.traces.by.pc %>%
  filter(exp_title == 'trial', aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  add.early.late.col(filter.early.and.late = TRUE)

nonreward.aligned.pop.traces.filtered = filter(nonreward.aligned.pop.traces, exp_title == 'trial',
       aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  dplyr::mutate(bout='non-reward')

# Examples of apopulation activity at approach

#example.animal = 'D-BR'
#example.animal = 'O-TR'
example.animal = 'F-BL'

example.reward.aligned.pop.traces = reward.aligned.pop.traces %>%
  filter(exp_title == 'trial', aligned_event_id >= 0) %>%
  filter(timestamp_from_end >= -5000 & timestamp_from_end <= 2000) %>%
  filter(animal == example.animal) %>%
  filter(day_desc %in% c('learning1 day#1', 'learning1 day#3', 'learning1 day#5')) %>%
  filter(location_set == 1)


g1 = example.reward.aligned.pop.traces %>%
  group_by(exp, exp_title, implant, day_desc, day_ordinal, location_set, timestamp_from_end) %>%
  dplyr::summarise(mean.val = mean(pmin(30, velocity)),
                   sem.val = sem(pmin(30, velocity)),
                   .groups = 'drop') %>%
  ggplot(aes(x=timestamp_from_end/1000)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_ribbon(aes(ymin=mean.val-sem.val, ymax=mean.val+sem.val), color='#555555', fill='#555555') +
  facet_wrap(. ~ day_desc, ncol=5) +
  xlab('Timestamp (s)') + ylab('Velocity') +
  gtheme


g2 = example.reward.aligned.pop.traces %>%
  ggplot(aes(x=timestamp_from_end / 1000, y=as.integer(event_seq_no))) +
  geom_raster(aes(fill=pmin(1, zscored_smooth_deconv_trace.mean))) +
  scale_fill_gradient(low='white', high='#333333') +
  geom_vline(xintercept = 0, linetype='dashed') +
  facet_wrap(. ~ day_desc, ncol=5) +
  xlab('Timestamp (s)') + ylab('Approach') + labs(fill='Population z-score') +
  scale_y_reverse() +
  gtheme +
  theme(legend.position = 'top')

g3 = example.reward.aligned.pop.traces %>%
  group_by(exp, exp_title, implant, day_desc, day_ordinal, location_set, timestamp_from_end) %>%
  dplyr::summarise(mean.val = mean(zscored_smooth_deconv_trace.mean),
                   sem.val = sem(zscored_smooth_deconv_trace.mean),
                   .groups = 'drop') %>%
  ggplot(aes(x=timestamp_from_end/1000)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_ribbon(aes(ymin=mean.val-sem.val, ymax=mean.val+sem.val), color='#555555', fill='#555555') +
  facet_wrap(. ~ day_desc, ncol=5) +
  xlab('Timestamp (s)') + ylab('Population z-score') +
  gtheme

plot_grid(g2, g3, g1, ncol=1, rel_heights = c(1.8, 1, 1))
figure.ggsave('FigS3AB-examaple_rew_aligned_activity.pdf',
              width=9, height=8.5)

#####################################################
## Beforetest trials
#####################################################

traces.ordinal.by.max.time = function(traces_v, timestamps_v, traces.grouping_v) {
  df = data.frame(traces=traces_v,
                  timestamps=timestamps_v,
                  traces.grouping = traces.grouping_v) 
  max.group.timestamp = df %>%
    dplyr::arrange(traces.grouping, timestamps) %>%
    dplyr::group_by(traces.grouping) %>%
    dplyr::mutate(smooth_traces = pracma::movavg(traces, 3)) %>%
    dplyr::summarise(maxval.timestamp = timestamps_v[which.max(traces)]) %>%
    arrange(desc(maxval.timestamp)) %>%
    dplyr::mutate(ordinal = dplyr::row_number())
  df %>% left_join(max.group.timestamp, by = 'traces.grouping') %>%
    pull(ordinal)
}

beforetest.aligned.pop.traces.by.pc %>%
  filter(implant== 'dCA1') %>%
  filter(timestamp_from_end >= -5000, timestamp_from_end <= 5000, timestamp_from_start >= 0) %>%
  dplyr::mutate(unique_event_id = paste(trial_id, event_seq_no)) %>%
  dplyr::group_by(implant, is.pc, animal) %>%
  dplyr::mutate(ordinal = traces.ordinal.by.max.time(cells.active.pct, 
                                                     timestamp_from_end, unique_event_id),
                max.val = max(cells.active.pct),
                norm.active.pct = cells.active.pct / max.val) %>%
  dplyr::group_by(implant, is.pc) %>%
  dplyr::mutate(animal_ordinal = traces.ordinal.by.group(animal)) %>%
  ggplot(aes(x=timestamp_from_end/1000, 
             y=animal_ordinal * 60 + ordinal,
             #y=ordinal, 
             fill = cells.active.pct)) +
  geom_tile(width=timebin.dur.msec*2/1000) +
  facet_wrap(implant ~ is.pc) +
  gtheme +
  geom_vline(xintercept = 0, linetype='dashed', color='#eeeeee') +
  scale_fill_gradientn(colours = jet.colours(7)) +
  #scale_fill_viridis_c(option = 'B') +
  gtheme

beforetest.aligned.pop.summary.by.pc = beforetest.aligned.pop.traces.by.pc %>%
  filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000, timestamp_from_end <= 2000) %>%
  group_by(implant, is.pc, timestamp_from_end) %>%
  dplyr::summarise(
    zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
    cells.active.pct.sem=sem(cells.active.pct),
    velocity.mean = mean(pmin(40, velocity)),
    velocity.sem = sem(pmin(40, velocity)),
    .groups='drop') %>%
  as.data.table()


beforetest.aligned.pop.summary.by.pc %>%
  ggplot(aes(x=timestamp_from_end / 1000)) +
  geom_ribbon(aes(ymin=cells.active.pct.mean-cells.active.pct.sem,
                  ymax=cells.active.pct.mean+cells.active.pct.sem,
                  group=is.pc, fill=is.pc), alpha=0.75) +
  geom_vline(xintercept=0, linetype='dashed') +
  facet_wrap(. ~ implant, scales='free') +
  scale_fill_manual(values=rev(side.two.coulours)) +
  labs(title='Approach to reward location (beforetest)') +
  ylab('Cells active (%)') +
  xlab('Time (s)') +
  ylim(c(5, 21)) +
  gtheme +
  theme(legend.position='top')
figure.ggsave('Fig3B-right-beforetest_reward_approach.pdf', height=5.1, width=6.4)

# Plot per animal
beforetest.aligned.pop.summary.by.pc.animal = beforetest.aligned.pop.traces.by.pc %>%
  filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000, timestamp_from_end <= 2000) %>%
  group_by(implant, animal, is.pc, timestamp_from_end) %>%
  dplyr::summarise(
    zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
    cells.active.pct.sem=sem(cells.active.pct),
    velocity.mean = mean(pmin(40, velocity)),
    velocity.sem = sem(pmin(40, velocity)),
    nsessions=n(),
    .groups='drop') %>%
  as.data.table()


beforetest.aligned.pop.summary.by.pc.animal %>%
  # dplyr::arrange(implant, is.pc, animal, timestamp_from_end) %>%
  # dplyr::group_by(implant, is.pc, animal) %>%
  # dplyr::mutate(cells.active.pct.norm = (cells.active.pct.mean / cells.active.pct.mean[1])) %>%
  mutate(is.pc_animal = paste(animal, is.pc)) %>%
  ggplot(aes(x=timestamp_from_end / 1000)) +
  geom_line(aes(y=cells.active.pct.mean,
                group=is.pc_animal, colour=is.pc), alpha=0.5) +
  #geom_ribbon(aes(ymin=cells.active.pct.mean-cells.active.pct.sem,
  #                ymax=cells.active.pct.mean+cells.active.pct.sem,
  #                group=is.pc_animal, fill=is.pc), alpha=0.2) +
  geom_vline(xintercept=0, linetype='dashed') +
  facet_grid(implant ~ is.pc, scales='free') +
  scale_colour_manual(values=rev(side.two.coulours)) +
  scale_fill_manual(values=rev(side.two.coulours)) +
  labs(title='Approach to reward location (beforetest)') +
  ylab('Cells active (%)') +
  xlab('Time (s)') +
  #ylim(c(5, 21)) +
  gtheme +
  theme(legend.position='top')
figure.ggsave('Fig3B-right-beforetest_reward_approach_per_animal.pdf', height=7.4, width=6.5)

beforetest.aligned.pop.summary.by.pc.animal %>%
  #filter(implant == 'dCA1') %>%
  dplyr::arrange(implant, is.pc, animal, timestamp_from_end) %>%
  dplyr::group_by(implant, is.pc, animal) %>%
  dplyr::mutate(cells.active.pct.norm = (cells.active.pct.mean / cells.active.pct.mean[1])) %>%
  #dplyr::group_by(implant, is.pc, animal) %>%
  #dplyr::mutate(cells.active.pct.norm = (cells.active.pct.mean) / max(cells.active.pct.mean)) %>%
  dplyr::group_by(implant, is.pc) %>%
  dplyr::mutate(animal_ordinal = traces.ordinal.by.group(animal)) %>%
  ggplot(aes(x=timestamp_from_end/1000, 
             y=animal_ordinal,
             fill = cells.active.pct.norm)) +
  geom_tile(width=timebin.dur.msec*1.1/1000) +
  facet_wrap(implant ~ is.pc) +
  gtheme +
  #geom_vline(xintercept = 0, linetype='dashed', color='#eeeeee') +
  geom_vline(xintercept = 0, linetype='dashed', color='#333333') +
  #scale_fill_gradientn(colours = jet.colours(7)) +
  scale_fill_gradient2(midpoint=1.0, high = "red", low = "blue") +
  #scale_fill_gradient(low='black', high='white') +
  #scale_fill_viridis_c(option = 'C') +
  gtheme
figure.ggsave('Fig3B-right-beforetest_reward_approach_per_animal.pdf', height=7.4, width=8.9)


# Stats on bouts arriving at reward in beforetest
timebin.width = 1000
beforetest.aligned.pop.traces.by.pc.binned = beforetest.aligned.pop.traces.by.pc %>%
  add.proximal.distal.col() %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial_id, day_desc, is.pc, timestamp_from_end_bin, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
beforetest.aligned.pop.traces.by.pc.binned$proximal.timestamp = as.factor(beforetest.aligned.pop.traces.by.pc.binned$proximal.timestamp)
beforetest.aligned.pop.traces.by.pc.binned$is.pc = as.factor(beforetest.aligned.pop.traces.by.pc.binned$is.pc)

beforetest.aligned.pop.traces.by.pc.binned.per.animal = beforetest.aligned.pop.traces.by.pc %>%
  add.proximal.distal.col() %>%
  #dplyr::group_by(implant, date, animal, exp, exp_title, trial_id, day_desc, is.pc, timestamp_from_end_bin) %>%
  dplyr::group_by(implant, animal, exp, exp_title, is.pc, timestamp_from_end_bin) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
beforetest.aligned.pop.traces.by.pc.binned.per.animal$proximal.timestamp =
  as.factor(beforetest.aligned.pop.traces.by.pc.binned.per.animal$proximal.timestamp)
beforetest.aligned.pop.traces.by.pc.binned.per.animal$is.pc =
  as.factor(beforetest.aligned.pop.traces.by.pc.binned.per.animal$is.pc)

implant.loc = 'dCA1'
beforetest.aligned.pop.traces.by.pc.binned %>%
  filter(implant == implant.loc) %>%
  ggplot(aes(x=proximal.timestamp, y=cells.active.pct, group=animal)) +
  geom_jitter(size=0.05, shape=1, alpha=0.35, width=0.2, height=0, color='#444444') +
  geom_line(data=subset(beforetest.aligned.pop.traces.by.pc.binned.per.animal, implant==implant.loc),
            size=0.5, color='#333333') +
  facet_wrap(implant ~ is.pc) +
  ylim(c(0, 35)) +
  gtheme +
  xlab('Time to reward location') + scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  ylab('Active cells %')
figure.ggsave(paste0('Fig3F-approach-stats-beforetest-',implant.loc, '.pdf'),
              width=4.5, height=4.4)

m.beforetest = lmer.test.print(subset(beforetest.aligned.pop.traces.by.pc.binned, implant == implant.loc),
                               cells.active.pct,
                               fixed.effects = is.pc * proximal.timestamp,
                               randef.str = '(1 | animal)',
                               diagnostics.groupvar=is.pc)
pairwise.post.hoc(m.beforetest, factor.interaction = c('is.pc:proximal.timestamp'))

printf('N approaches = ')
beforetest.aligned.pop.traces.by.pc.binned %>%
  ungroup() %>%
  dplyr::select(implant, trial_id, aligned_event_id) %>%
  dplyr::distinct() %>%
  dplyr::count(implant)

models = create.bayes.lm.pair(subset(beforetest.aligned.pop.traces.by.pc.binned, implant == implant.loc & is.pc=='TRUE'),
                              formula.full = cells.active.pct ~ 1 + proximal.timestamp + animal,
                              formula.null = cells.active.pct ~ 1 + animal,
                              whichRandom = 'animal')
models$full / models$null
calc.pair.95CI(models$full, pair.vars = c('proximal.timestamp-0', 'proximal.timestamp-1'))
group_by(beforetest.aligned.pop.traces.by.pc.binned, implant, is.pc, proximal.timestamp) %>%
  dplyr::summarise(mean(cells.active.pct), sem(cells.active.pct), mean(zscored_smooth_deconv_trace.mean), sem(zscored_smooth_deconv_trace.mean), n())


#####################################################
# Arriving at reward during learning
#####################################################
reward.pop.traces.aligned.comparison = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  group_by(implant, timestamp_from_end, bout) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_trace.mean),
                   zscored_trace.mean=mean(zscored_trace.mean),
                   velocity.mean=mean(velocity, na.rm=TRUE),
                   velocity.sem=sem(velocity),
                   cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop')

reward.pop.traces.aligned.comparison %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=bout)) +
  geom_ribbon(aes(ymin=zscored_deconv_trace.mean-zscored_deconv_trace.sem,
                 ymax=zscored_deconv_trace.mean+zscored_deconv_trace.sem,
                 fill=bout), alpha=0.6) +
  # geom_ribbon(aes(ymin=velocity.mean-velocity.sem,
  #                 ymax=velocity.mean+velocity.sem,
  #                 fill=bout), alpha=0.6) +
  facet_wrap(. ~ implant, scales='free') +
  scale_fill_manual(values=three.colours) +
  scale_color_manual(values=three.colours) +
  scale_alpha_manual(values=c(0.6, 0.6, 0.6)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  labs(title='Learning-changed activity reward approach') +
  xlab('Time (s)') +
  xlim(c(-5, 2)) +
  ylab('deconvolved dF/F (z-score)') +
  gtheme +
  theme(legend.position = 'top')
figure.ggsave('Fig3BE_left-reward_approach.pdf',
              height=6.5, width=8)

# Per animal
reward.pop.traces.aligned.comparison.per.animal = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  group_by(implant, animal, timestamp_from_end, bout) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_trace.mean),
                   zscored_trace.mean=mean(zscored_trace.mean),
                   velocity.mean=mean(velocity, na.rm=TRUE),
                   cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop')

three.colours.activity = c('#e8ad91', '#d75b43', '#999999')
reward.pop.traces.aligned.comparison.per.animal %>%
  mutate(bout_animal = paste(bout, animal)) %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=bout_animal)) +
  geom_line(aes(y=zscored_deconv_trace.mean,
                color=bout), alpha=0.6) +
  geom_ribbon(aes(ymin=zscored_deconv_trace.mean-zscored_deconv_trace.sem,
                  ymax=zscored_deconv_trace.mean+zscored_deconv_trace.sem,
                  fill=bout), alpha=0.3) +
  facet_grid(implant ~ bout, scales='free') +
  scale_fill_manual(values=three.colours.activity) +
  scale_color_manual(values=three.colours.activity) +
  scale_alpha_manual(values=c(0.6, 0.6, 0.6)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  labs(title='Learning-changed activity reward approach') +
  xlab('Time (s)') +
  xlim(c(-5, 2)) +
  ylab('deconvolved dF/F (z-score)') +
  gtheme +
  theme(legend.position = 'top')
figure.ggsave('Fig3BE_left-reward_approach_per_animal.pdf',
              height=7.3, width=9)


# Stats on activity during bouts at reward
timebin.width = 1000
reward.aligned.pop.traces.earlylate.binned = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  add.proximal.distal.col(filter.proximal.distal = TRUE) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial, trial_id, day_desc,
                  proximal.timestamp, bout, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
reward.aligned.pop.traces.earlylate.binned$bout = as.factor(reward.aligned.pop.traces.earlylate.binned$bout)
reward.aligned.pop.traces.earlylate.binned$proximal.timestamp = as.factor(reward.aligned.pop.traces.earlylate.binned$proximal.timestamp)

reward.aligned.pop.traces.earlylate.binned.per.day = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  add.proximal.distal.col(filter.proximal.distal = TRUE) %>%
  dplyr::group_by(implant, animal, exp, exp_title, proximal.timestamp, bout) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2,
                aday=paste(animal, day_desc, sep='_'))
reward.aligned.pop.traces.earlylate.binned.per.day$bout =
  as.factor(reward.aligned.pop.traces.earlylate.binned.per.day$bout)
reward.aligned.pop.traces.earlylate.binned.per.day$proximal.timestamp =
  as.factor(reward.aligned.pop.traces.earlylate.binned.per.day$proximal.timestamp)

print('N approaches = ')
reward.aligned.pop.traces.earlylate.binned %>%
  ungroup() %>%
  dplyr::select(implant, bout, trial_id, aligned_event_id) %>%
  dplyr::distinct() %>%
  dplyr::count(implant, bout)


implant.loc = 'dCA1'
reward.aligned.pop.traces.earlylate.binned %>%
  filter(implant == implant.loc) %>%
  ggplot(aes(x=proximal.timestamp, y=zscored_deconv_trace.mean)) +
  geom_jitter(size=0.05, shape=1, alpha=0.3, width=0.3, height=0, color='#444444') +
  geom_line(data=subset(reward.aligned.pop.traces.earlylate.binned.per.day, implant==implant.loc),
            mapping=aes(group=aday),
            size=0.5, color='#333333') +
  facet_wrap(implant ~ bout) +
  ylim(c(-0.1, 0.85)) +
  gtheme +
  xlab('Time to reward location') +
  scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  ylab('Population z-score')
figure.ggsave(paste0('FigS3C-approach_stats_activity_',implant.loc, '.pdf'),
              width=5.6, height=4.4)

m.earlylate = lmer.test.print(
  subset(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc),
  log(1.0 + zscored_deconv_trace.mean),
  fixed.effects = bout * proximal.timestamp,
  diagnostics.groupvar=bout)
group_by(subset(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc), bout, proximal.timestamp) %>%
  dplyr::summarise(mean(zscored_deconv_trace.mean), sem(zscored_deconv_trace.mean), n())
pairwise.post.hoc(m.earlylate, factor.interaction = c('bout:proximal.timestamp'))

# Significantly lower velocity at reward
m.velocity = lmer.test.print(
  subset(reward.aligned.pop.traces.earlylate.binned, bout=='late'),
  log(1.0 + velocity),
  fixed.effects = proximal.timestamp,
  randef.str = '(1 | animal)',
  diagnostics.groupvar=animal)
group_by(subset(reward.aligned.pop.traces.earlylate.binned, bout=='late'), proximal.timestamp) %>%
  dplyr::summarise(mean(velocity), sem(velocity), n())

models = create.bayes.lm.pair(
  filter(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc, bout == 'non-reward') %>%
    mutate(val = log(1.0 + zscored_deconv_trace.mean)),
  formula.full = zscored_deconv_trace.mean ~ 1 + proximal.timestamp + animal,
  formula.null = zscored_deconv_trace.mean ~ 1 + animal,
  whichRandom = 'animal', iterations = 50000)
models$full / models$null
calc.pair.95CI(models$full, pair.vars = c('proximal.timestamp-FALSE', 'proximal.timestamp-TRUE'),
               ytransform = function(x) {exp(x) - 1 }, show.percent.change = FALSE)

cross.corr.vals = function(trace, velocity){
  lag.cor = ccf(trace, velocity, lag.max=1000/timebin.dur.msec * 1 * 0.6)
  max.lag.index = which.max(lag.cor$acf)
  min.lag.index = which.min(lag.cor$acf)
  dplyr::tibble(max.cor.lag.ms=lag.cor$lag[max.lag.index] * timebin.dur.msec,
                min.cor.lag.ms=lag.cor$lag[min.lag.index] * timebin.dur.msec,
                max.cor=lag.cor$acf[max.lag.index],
                min.cor=lag.cor$acf[min.lag.index])
}

# Find minimum or maximum correlation between velocity and population activity
reward.aligned.pop.traces.late.vel.cor = reward.aligned.pop.traces.earlylate %>%
  filter(bout=='late') %>%
  filter(timestamp_from_end <= 0, timestamp_from_end >= -5000) %>%
  group_by(implant, animal, day_desc, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   velocity.mean=mean(velocity, na.rm=TRUE),
                   acceleration.mean=mean(acceleration, na.rm=TRUE),
                   .groups='drop') %>%
  arrange(implant, animal, day_desc, timestamp_from_end) %>%
  group_by(implant, animal, day_desc) %>%
  #arrange(implant, animal, exp, exp_title, date, trial, trial_id, event_seq_no, timestamp_from_end) %>%
  #group_by(implant, animal, exp, exp_title, date, trial, trial_id, event_seq_no) %>%
  dplyr::summarise(cross.corr.vals(zscored_deconv_trace.mean, velocity.mean)) %>%
  mutate(cor.lag.ms = ifelse(implant=='dCA1', min.cor.lag.ms, max.cor.lag.ms),
         cor.val=ifelse(implant=='dCA1', min.cor, max.cor))

cross.corr.vals.line = function(trace, velocity){
  lag.cor = ccf(trace, velocity, 
                lag.max=1000/timebin.dur.msec * 1 * 0.6,
                plot=FALSE)
  max.lag.index = which.max(lag.cor$acf)
  min.lag.index = which.min(lag.cor$acf)
  data.frame(lag.ms=lag.cor$lag * timebin.dur.msec,
             acf=lag.cor$acf)
}

reward.aligned.pop.traces.late.vel.cor.line = reward.aligned.pop.traces.earlylate %>%
  filter(bout=='late') %>%
  filter(timestamp_from_end <= 0, timestamp_from_end >= -5000) %>%
  group_by(implant, animal, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   velocity.mean=mean(acceleration, na.rm=TRUE),
                   acceleration.mean=mean(acceleration, na.rm=TRUE),
                   .groups='drop') %>%
  #arrange(implant, animal, day_desc, is.pc, timestamp_from_end) %>%
  #group_by(implant, animal, day_desc, is.pc) %>%
  arrange(implant, animal, timestamp_from_end) %>%
  group_by(implant, animal) %>%
  dplyr::summarise(cross.corr.vals.line(zscored_deconv_trace.mean, acceleration.mean))

reward.aligned.pop.traces.late.vel.cor.mean.line = reward.aligned.pop.traces.late.vel.cor.line %>%
  group_by(implant, lag.ms) %>%
  dplyr::summarise(acf.mean=mean(acf), .groups='drop')

reward.aligned.pop.traces.late.vel.cor.mean.line %>%
  group_by(implant) %>%
  dplyr::summarise(max.lag.ms = lag.ms[which.max(acf.mean)],
                   max.acf = max(acf.mean),
                   min.lag.ms = lag.ms[which.min(acf.mean)],
                   min.acf = min(acf.mean))

reward.aligned.pop.traces.late.vel.cor.line %>%
  mutate(animal_line=paste(implant, animal)) %>%
  ggplot(aes(x=-lag.ms/1000, y=acf)) +
  geom_line(aes(group=animal_line), color='#777777', size=0.3) +
  geom_line(data=reward.aligned.pop.traces.late.vel.cor.mean.line,
            mapping=aes(y=acf.mean),
            color='#333333', size=0.8) +
  geom_vline(xintercept = 0, linetype='dashed', color='#333333') +
  geom_hline(yintercept = 0, linetype='dashed', color='#333333') +
  scale_x_continuous(breaks=seq(-0.4,0.4,0.4)) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.5)) +
  facet_grid(implant ~ .) +
  xlab('Acceleration lag (s)') + ylab('Cross-correlation') +
  gtheme
figure.ggsave('FigS3-acceleration-cor.pdf', width=4.0, height=4.4)

reward.aligned.pop.traces.late.vel.cor.mean = reward.aligned.pop.traces.late.vel.cor %>%
  group_by(implant) %>%
  dplyr::summarise(mean.lag.s = mean(cor.lag.ms / 1000),
                   sem.lag.s = sem(cor.lag.ms / 1000),
                   mean.cor.val = mean(cor.val))

reward.aligned.pop.traces.late.vel.cor %>%
  ggplot(aes(x=cor.lag.ms/1000, y=cor.val)) +
  geom_vline(xintercept = 0, linetype='dashed', color='#333333') +
  geom_hline(yintercept = 0, linetype='dashed', color='#333333') +
  geom_point(shape=1, color='#777777') +
  geom_point(data=reward.aligned.pop.traces.late.vel.cor.mean,
             mapping=aes(x=mean.lag.s, y=mean.cor.val),
             shape=3, color='#333333', size=2, stroke=1.5) +
  facet_grid(. ~ implant) +
  xlab('Lag (s)') + ylab('Correlation') +
  scale_x_continuous(breaks=c(-0.5,0,0.5)) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.5)) +
  gtheme

##############################################################
# Active cells as mice arrived at reward during late learning
##############################################################
# Activity at reward bouts in place vs non-place cells, trials after learnt
reward.aligned.pop.traces.earlylate.by.pc %>%
  group_by(implant, is.pc, is.early.learning, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   cells.active.pct.mean=mean(cells.active.pct),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop') %>%
  filter(!is.early.learning) %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=is.pc)) +
  geom_ribbon(aes(
    ymin=cells.active.pct.mean-cells.active.pct.sem,
    ymax=cells.active.pct.mean+cells.active.pct.sem,
    fill=is.pc), alpha=0.75) +
  facet_wrap(. ~ implant, scales='free') +
  scale_fill_manual(values=side.two.coulours) +
  xlim(c(-5, 2)) +
  ylim(c(5, 25)) +
  labs(title='Approach to rewarded location learnt') +
  ylab('Cells active (%)') +
  xlab('Time (s)') +
  gtheme +
  theme(legend.position='top')
figure.ggsave('Fig3BE_middle-learnt_reward_approach.svg',
              height=6.5, width=8)

# Per animal
reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(ncells>1) %>%
  group_by(implant, animal, is.pc, is.early.learning, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   cells.active.pct.mean=mean(cells.active.pct),
                   cells.active.pct.sem=sem(cells.active.pct), 
                   nsessions=n(),
                   .groups='drop') %>%
  filter(!is.early.learning) %>%
  mutate(animal_pc=paste(animal, is.pc)) %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=animal_pc)) +
  geom_line(aes(y=cells.active.pct.mean,
                colour=is.pc), 
            alpha=0.5) +
  geom_ribbon(aes(
    ymin=pmax(0, cells.active.pct.mean-cells.active.pct.sem),
    ymax=cells.active.pct.mean+cells.active.pct.sem,
    fill=is.pc), alpha=0.2) +
  geom_vline(xintercept = 0, linetype='dashed') +
  facet_grid(implant ~ is.pc, scales='free') +
  scale_colour_manual(values=rev(side.two.coulours)) +
  scale_fill_manual(values=rev(side.two.coulours)) +
  xlim(c(-5, 2)) +
  #ylim(c(0, 35)) +
  labs(title='Approach to rewarded location learnt') +
  ylab('Cells active (%)') +
  xlab('Time (s)') +
  gtheme +
  theme(legend.position='top')
figure.ggsave('Fig3BE_middle-learnt_reward_approach_peranimal.pdf', height=7.4, width=6.5)

# Stats on bouts arriving at reward in late learning
timebin.width = 1000
reward.aligned.pop.traces.by.pc.binned.per.day = reward.aligned.pop.traces.by.pc %>%
  filter(exp_title == 'trial', aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  add.early.late.col(filter.early.and.late = FALSE) %>%
  add.proximal.distal.col(filter.proximal.distal = TRUE) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  #dplyr::group_by(implant, date, animal, exp, exp_title, day_desc, is.late.learning, is.pc, proximal.timestamp, timestamp_from_end_bin) %>%
  dplyr::group_by(implant, animal, exp, exp_title, is.late.learning, is.pc, proximal.timestamp, timestamp_from_end_bin) %>%
  dplyr::summarise_all(mean, na.rm=TRUE) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2,
                aday=paste(animal, day_desc, sep='_'))
reward.aligned.pop.traces.by.pc.binned.per.day$is.pc =
  as.factor(reward.aligned.pop.traces.by.pc.binned.per.day$is.pc)
reward.aligned.pop.traces.by.pc.binned.per.day$proximal.timestamp =
  as.factor(reward.aligned.pop.traces.by.pc.binned.per.day$proximal.timestamp)

reward.aligned.pop.traces.late.by.pc.binned.per.day = filter(reward.aligned.pop.traces.by.pc.binned.per.day, is.late.learning)

reward.aligned.pop.traces.late.by.pc.binned = filter(reward.aligned.pop.traces.earlylate.by.pc, is.late.learning) %>%
  add.proximal.distal.col(filter.proximal.distal = TRUE) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial, trial_id, day_desc, is.pc,
                  proximal.timestamp, timestamp_from_end_bin, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
reward.aligned.pop.traces.late.by.pc.binned$is.pc = as.factor(reward.aligned.pop.traces.late.by.pc.binned$is.pc)
reward.aligned.pop.traces.late.by.pc.binned$proximal.timestamp = as.factor(reward.aligned.pop.traces.late.by.pc.binned$proximal.timestamp)

implant.loc = 'dCA1'
reward.aligned.pop.traces.late.by.pc.binned %>%
  filter(implant == implant.loc) %>%
  ggplot(aes(x=proximal.timestamp, y=cells.active.pct)) +
  geom_jitter(size=0.05, shape=1, alpha=0.35, width=0.25, height=0, color='#444444') +
  geom_line(data=subset(reward.aligned.pop.traces.late.by.pc.binned.per.day, implant==implant.loc),
            mapping=aes(group=aday),
            size=0.4, color='#333333') +
  facet_wrap(implant ~ is.pc) +
  ylim(c(0, 42)) +
  gtheme +
  xlab('Time to reward location') +
  scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  ylab('Active cells %')
figure.ggsave(paste0('FigS3D-approach-stats-latelearning-', implant.loc, '.pdf'),
              width=4.0, height=4.4)

m.pc.proximal = lmer.test.print(
  subset(reward.aligned.pop.traces.late.by.pc.binned, implant==implant.loc),
  cells.active.pct,
  diagnostics.groupvar = is.pc,
  fixed.effects = is.pc * proximal.timestamp)
group_by(subset(reward.aligned.pop.traces.late.by.pc.binned, implant==implant.loc), is.pc, proximal.timestamp) %>%
  dplyr::summarise(mean(cells.active.pct), sem(cells.active.pct), n())
pairwise.post.hoc(m.pc.proximal, factor.interaction = c('is.pc:proximal.timestamp'))

models = create.bayes.lm.pair(subset(reward.aligned.pop.traces.late.by.pc.binned, implant==implant.loc & is.pc=='TRUE'),
                              formula.full = cells.active.pct ~ 1 + proximal.timestamp + animal,
                              formula.null = cells.active.pct ~ 1 + animal,
                              whichRandom = 'animal')
models$full / models$null
calc.pair.95CI(models$full,
               pair.vars = c("proximal.timestamp-FALSE", "proximal.timestamp-TRUE"),
               show.percent.change = TRUE)


# Activity at reward bouts in place vs non-place cells, trials after learnt, as function of distance
dist.bin.width = 2
ndist.bins = 20
reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(implant == 'vCA1') %>%
  filter(timestamp_from_end >= -5000, timestamp_from_start >= 0, timestamp_from_end <= 0) %>%
  mutate(dist_approached_rew_bin = round(dist_approached_rew / dist.bin.width)) %>%
  group_by(implant, is.pc, is.early.learning, dist_approached_rew_bin, animal) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   cells.active.pct.mean=mean(cells.active.pct),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop') %>%
  filter(!is.early.learning) %>%
  dplyr::mutate(line_id=paste(implant, is.pc, animal)) %>%
  filter(dist_approached_rew_bin <= ndist.bins) %>%
  arrange(implant, is.pc, animal, desc(dist_approached_rew_bin)) %>%
  group_by(implant, is.pc, animal, line_id) %>%
  dplyr::mutate(cells.active.pct.norm = cells.active.pct.mean / cells.active.pct.mean[1] * 100) %>%
  ggplot(aes(x=dist_approached_rew_bin * dist.bin.width, group=line_id)) +
  # geom_ribbon(aes(
  #   ymin=cells.active.pct.mean-cells.active.pct.sem,
  #   ymax=cells.active.pct.mean+cells.active.pct.sem,
  #   fill=is.pc), alpha=0.3) +
  geom_line(aes(y=cells.active.pct.norm, color=is.pc)) +
  facet_wrap(implant ~ is.pc, scales='free') +
  scale_fill_manual(values=rev(side.two.coulours)) +
  scale_color_manual(values=rev(side.two.coulours)) +
  scale_x_reverse(limit=c(ndist.bins * dist.bin.width,0)) +
  #xlim(c(0, 40)) +
  #ylim(c(5, 33)) +
  #ylim(c(5, 28)) +
  labs(title='Approach to rewarded location learnt') +
  ylab('Cells active (%)') +
  xlab('Distance (cm)') +
  gtheme +
  theme(legend.position='top')
figure.ggsave('FigS3E-learnt_reward_approach_over_distance-vca1b.pdf', height=5.5, width=6.6)
#figure.ggsave('FigS3E-learnt_reward_approach_over_distance.pdf', height=5.2, width=6.6)

traces.ordinal.by.group = function(grouping.factor) {
  as.factor(grouping.factor) %>% as.numeric()
}

reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(implant== 'vCA1') %>%
  filter(timestamp_from_end >= -10000, timestamp_from_start >= 0, timestamp_from_end < 0) %>%
  mutate(dist_approached_rew_bin = round(dist_approached_rew / dist.bin.width),
         unique_event_id = paste(trial_id, event_seq_no)) %>%
  dplyr::group_by(implant, is.pc, animal) %>%
  dplyr::mutate(ordinal = traces.ordinal.by.max.time(cells.active.pct, dist_approached_rew_bin, unique_event_id),
                max.val = max(cells.active.pct),
                norm.active.pct = cells.active.pct / max.val) %>%
  dplyr::group_by(implant, is.pc) %>%
  dplyr::mutate(animal_ordinal = traces.ordinal.by.group(animal)) %>%
  ggplot(aes(x=dist_approached_rew_bin * dist.bin.width, y=animal_ordinal * 80 + ordinal, fill = norm.active.pct)) +
  geom_tile(width=dist.bin.width * 2) +
  facet_wrap(implant ~ is.pc, scales = 'free_y') +
  gtheme +
  scale_fill_gradientn(colours = jet.colours(7)) +
  scale_x_reverse(limit=c(40,0)) +
  gtheme

reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(timestamp_from_end >= -5000, timestamp_from_start >= 0, timestamp_from_end < 0) %>%
  mutate(dist_approached_rew_bin = round(dist_approached_rew / dist.bin.width)) %>%
  group_by(implant, is.pc, is.early.learning, dist_approached_rew_bin, animal) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   cells.active.pct.mean=mean(cells.active.pct),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop') %>%
  filter(!is.early.learning) %>%
  #dplyr::mutate(line_id=paste(implant, is.pc, animal)) %>%
  dplyr::group_by(implant, is.pc) %>%
  #dplyr::mutate(ordinal = dplyr::cur_group_id()) %>%
  #dplyr::mutate(ordinal = traces.ordinal.by.group(animal)) %>%
  dplyr::mutate(ordinal = traces.ordinal.by.max.time(cells.active.pct.mean, dist_approached_rew_bin, animal),
                norm.active.pct = cells.active.pct.mean / max(cells.active.pct.mean)) %>%
  ggplot(aes(x=dist_approached_rew_bin * dist.bin.width, y=ordinal, fill = norm.active.pct)) +
  geom_tile(width=dist.bin.width) +
  facet_wrap(implant ~ is.pc, scales = 'free_y') +
  gtheme +
  scale_fill_gradientn(colours = jet.colours(7)) +
  scale_x_reverse(limit=c(40,0)) +
  gtheme

# Correlation with acceleration
reward.aligned.pop.traces.late.vel.cor = reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(is.late.learning) %>%
  filter(timestamp_from_end <= 0, timestamp_from_end >= -5000) %>%
  group_by(implant, animal, day_desc, is.pc, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   velocity.mean=mean(velocity, na.rm=TRUE),
                   acceleration.mean=mean(acceleration, na.rm=TRUE),
                   .groups='drop') %>%
  arrange(implant, animal, day_desc, is.pc, timestamp_from_end) %>%
  group_by(implant, animal, day_desc, is.pc) %>%
  #dplyr::summarise(cross.corr.vals.line(zscored_deconv_trace.mean, acceleration.mean))
  dplyr::summarise(cross.corr.vals(zscored_deconv_trace.mean, acceleration.mean)) %>%
  mutate(cor.lag.ms = ifelse(implant=='dCA1', min.cor.lag.ms, max.cor.lag.ms),
         cor.val=ifelse(implant=='dCA1', min.cor, max.cor))

# Mean line of cross-correlation 
reward.aligned.pop.traces.late.vel.cor.line = reward.aligned.pop.traces.earlylate.by.pc %>%
  filter(is.late.learning) %>%
  filter(timestamp_from_end <= 0, timestamp_from_end >= -5000) %>%
  group_by(implant, animal, is.pc, timestamp_from_end) %>%
  dplyr::summarise(zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   velocity.mean=mean(velocity, na.rm=TRUE),
                   acceleration.mean=mean(acceleration, na.rm=TRUE),
                   .groups='drop') %>%
  #arrange(implant, animal, day_desc, is.pc, timestamp_from_end) %>%
  #group_by(implant, animal, day_desc, is.pc) %>%
  arrange(implant, animal, is.pc, timestamp_from_end) %>%
  group_by(implant, animal, is.pc) %>%
  dplyr::summarise(cross.corr.vals.line(zscored_deconv_trace.mean, acceleration.mean))

reward.aligned.pop.traces.late.vel.cor.mean.line = reward.aligned.pop.traces.late.vel.cor.line %>%
  group_by(implant, is.pc, lag.ms) %>%
  dplyr::summarise(acf.mean=mean(acf), .groups='drop')

reward.aligned.pop.traces.late.vel.cor.mean.line %>%
  group_by(implant, is.pc) %>%
  dplyr::summarise(max.lag.ms = lag.ms[which.max(acf.mean)],
                   max.acf = max(acf.mean),
                   min.lag.ms = lag.ms[which.min(acf.mean)],
                   min.acf = min(acf.mean))

reward.aligned.pop.traces.late.vel.cor.line %>%
  mutate(animal_line=paste(implant, animal)) %>%
  ggplot(aes(x=lag.ms/1000, y=acf, color=is.pc)) +
  geom_line(aes(group=animal_line), #color='#333333', 
            alpha=0.5,
            size=0.3) +
  geom_line(data=reward.aligned.pop.traces.late.vel.cor.mean.line,
            mapping=aes(y=acf.mean),
            #color='#333333', 
            size=0.7) +
  geom_vline(xintercept = 0, linetype='dashed', color='#333333') +
  geom_hline(yintercept = 0, linetype='dashed', color='#333333') +
  scale_x_continuous(breaks=seq(-0.4,0.4,0.4)) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.5)) +
  scale_color_manual(values=rev(side.two.coulours)) +
  facet_grid(implant ~ is.pc) +
  xlab('Lag (s)') + ylab('Cross-correlation') +
  gtheme + theme(legend.position = 'none')
figure.ggsave('FigS3-acceleration-cor.pdf', width=6.0, height=4.9)

reward.aligned.pop.traces.late.vel.cor.mean = reward.aligned.pop.traces.late.vel.cor %>%
  group_by(implant, is.pc) %>%
  dplyr::summarise(mean.lag.s = mean(cor.lag.ms / 1000),
                   sem.lag.s = sem(cor.lag.ms / 1000),
                   mean.cor.val = mean(cor.val))

reward.aligned.pop.traces.late.vel.cor %>%
  ggplot(aes(x=cor.lag.ms/1000, y=cor.val)) +
  geom_vline(xintercept = 0, linetype='dashed', color='#333333') +
  geom_hline(yintercept = 0, linetype='dashed', color='#333333') +
  geom_point(shape=1, color='#777777') +
  geom_point(data=reward.aligned.pop.traces.late.vel.cor.mean,
             mapping=aes(x=mean.lag.s, y=mean.cor.val),
             shape=3, color='#333333', size=2, stroke=1.5) +
  facet_grid(implant ~ is.pc) +
  xlab('Lag (s)') + ylab('Correlation') +
  scale_x_continuous(breaks=c(-0.5,0,0.5)) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.5)) +
  gtheme

#figure.ggsave('S3-acceleration-cor.pdf', width=6.0, height=4.4)


##############################################
### Correlation of ramping with performance
##############################################
# Day mean difference between proximal and timestamp timestamps at reward
reward.proximity.diff.by.pc.binned.per.day = reward.aligned.pop.traces.by.pc %>%
  add.proximal.distal.col(filter.proximal.distal = TRUE) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial, trial_id, event_seq_no, day_desc,
                  is.pc, day_ordinal, location_set, exp_day_ordinal, proximal.timestamp) %>%
  dplyr::summarise_all(mean, .groups='drop') %>%
  arrange(implant, date, animal, exp, exp_title, trial, trial_id, aligned_event_id, day_desc, 
             is.pc, day_ordinal, location_set, exp_day_ordinal, timestamp_from_end) %>%
  pivot_wider(id_cols=c(implant, date, animal, exp, exp_title, trial, trial_id, event_seq_no, day_desc, 
                        is.pc, day_ordinal, location_set, exp_day_ordinal),
              values_from=c(zscored_smooth_deconv_trace.mean, cells.active.pct),
              values_fn=list) %>%
  filter(map_lgl(zscored_smooth_deconv_trace.mean_, ~ length(.x) > 1)) %>%
  mutate(zscored_smooth_deconv_trace.diff = map_dbl(zscored_smooth_deconv_trace.mean_, diff),
         cells.active.pct.diff = map_dbl(cells.active.pct_, diff)) %>%
  dplyr::select(-c(cells.active.pct_, zscored_smooth_deconv_trace.mean_)) %>%
  group_by(implant, date, animal, exp, exp_title, day_desc, is.pc) %>%
  dplyr::summarise(across(.fns = mean), .groups='drop')

# Learning performance data
df = read.csv(file.path(base_dir, 'trial_stats.csv'))
df$date = as.Date(df$date)

d = filter(df, exp_title == 'trial')
d.day = d %>%
  group_by(animal, date, exp, exp_title, day_ordinal, day_desc) %>%
  dplyr::summarise(mean.logdist = log(total_dist*perc2dist/100) %>% mean)

test.trial.df = filter(df, exp_title == 'beforetest') %>%
  mutate(dist_0_120s = dist_0_30s + dist_30_60s + dist_60_90s + dist_90_120s,
         crossings_0_120s = crossings_0_30s + crossings_30_60s + crossings_60_90s + crossings_90_120s,
         norm_crossings_120s = crossings_0_120s / dist_0_120s * perc2dist * 100)

# Correlation with day-performance
## Activity diff corr with performance
library(effects)
lmer.test.activity.diff = function(df) {
  m = lmerTest::lmer(cells.active.pct.diff ~ mean.logdist + (1 | animal),
                     data = df,
                     REML = TRUE)
  print(summary(m))
  print(anova(m, refit=FALSE, ddf='Satterthwaite'))
  
  models = create.bayes.lm.pair(df,
                                formula.full = cells.active.pct.diff ~ 1 + mean.logdist + animal,
                                formula.null = cells.active.pct.diff ~ 1 + animal,
                                whichRandom = 'animal')
  print(models$full / models$null)
  samples = BayesFactor::posterior(models$full, iterations = 10000, columnFilter="^animal$")
  samples[, 'mean.logdist'] %>% quantile(c(0.5, 0.025, 0.975)) %>% print()
  list(model=m,
       effect=force(Effect(c('mean.logdist'), m)))
}
# Workaround for issue with effects:Effect in the function
reward.proximity.diff.perf.cor.by.pc.day = 
  left_join(reward.proximity.diff.by.pc.binned.per.day, d.day) %>%
  # Remove days with incomplete day-trials
  filter(!day_desc %in% c('learning3 day#3', 'learning4 day#3') )

df = reward.proximity.diff.perf.cor.by.pc.day
m.rew.prox.day.perf.cor.vca1.nonpc = lmer.test.activity.diff(
  filter(reward.proximity.diff.perf.cor.by.pc.day, implant=='vCA1', is.pc == FALSE))
m.rew.prox.day.perf.cor.vca1.pc = lmer.test.activity.diff(
  filter(reward.proximity.diff.perf.cor.by.pc.day, implant=='vCA1', is.pc == TRUE))
m.rew.prox.day.perf.cor.dca1.nonpc = lmer.test.activity.diff(
  filter(reward.proximity.diff.perf.cor.by.pc.day, implant=='dCA1', is.pc == FALSE))
m.rew.prox.day.perf.cor.dca1.pc = lmer.test.activity.diff(
  filter(reward.proximity.diff.perf.cor.by.pc.day, implant=='dCA1', is.pc == TRUE))

m.rew.prox.perf.cor.effects = bind_rows(
  as.data.frame(m.rew.prox.day.perf.cor.vca1.nonpc$effect) %>% mutate(implant='vCA1', is.pc=FALSE),
  as.data.frame(m.rew.prox.day.perf.cor.vca1.pc$effect) %>% mutate(implant='vCA1', is.pc=TRUE),
  as.data.frame(m.rew.prox.day.perf.cor.dca1.nonpc$effect) %>% mutate(implant='dCA1', is.pc=FALSE),
  as.data.frame(m.rew.prox.day.perf.cor.dca1.pc$effect) %>% mutate(implant='dCA1', is.pc=TRUE),
)


reward.proximity.diff.perf.cor.by.pc.day %>%
  ggplot(aes(x=-mean.logdist)) +
  geom_ribbon(data=m.rew.prox.perf.cor.effects,
              mapping=aes(ymin=lower, ymax=upper),
              alpha=0.2, fill='#333333') +
  geom_line(data=m.rew.prox.perf.cor.effects,
            mapping=aes(y=fit, color=is.pc)) +
  geom_point(aes(y=cells.active.pct.diff, 
                 color=is.pc), shape=1, size=1) +
  geom_hline(yintercept = 0, linetype='dashed') +
  xlab('Trial performance (-log mean run distance (m))') +
  ylab('Active cells difference (%)') +
  scale_color_manual(values=rev(side.two.coulours)) +
  facet_grid(implant ~  is.pc, scales='free') +
  gtheme

figure.ggsave('Fig3C-active_cells_diff_cor_with_performance.pdf', height=6.8, width=8.5)


#####################################
# Comparison of runnning behaviour
#####################################
# Plot running velocity during approach during learning and during beforetest
# dCA1 vs vCA1
reward.aligned.pop.traces.earlylate %>%
  filter(bout == 'late') %>%
  filter(timestamp_from_end <= 2000, timestamp_from_end >= -5000) %>%
  filter(velocity <= 50) %>% # Skip incorrectly high velocity
  group_by(exp_title, implant, animal, date, day_desc, trial_id, timestamp, event_seq_no, timestamp_from_end) %>%
  dplyr::slice(1) %>% 
  group_by(exp_title,  implant, animal, timestamp_from_end) %>%
  dplyr::summarise(velocity.mean = mean(velocity),
                   velocity.sem = sem(velocity)) %>%
  ggplot(aes(x=timestamp_from_end/1000, group=animal)) +
  geom_ribbon(aes(ymin=velocity.mean-velocity.sem,
                  ymax=velocity.mean+velocity.sem),
              alpha=0.3) +
  geom_line(aes(y=velocity.mean), color='#666666', size=0.2) +
  geom_vline(xintercept=0, linetype='dashed', color='#333333') +
  geom_hline(yintercept=4, linetype='dashed', color='#333333') +
  facet_grid(implant ~ .) +
  xlab('Time (s)') + ylab('Speed (cm/s)') +
  gtheme

figure.ggsave('Figure_S3D-speed-approach-learning.pdf', width=4.5, height=4.5)

beforetest.aligned.pop.traces %>%
  filter(timestamp_from_end <= 2000, timestamp_from_end >= -5000) %>%
  filter(velocity <= 50) %>% # Skip incorrectly high velocity
  group_by(exp_title, implant, animal, date, day_desc, trial_id, timestamp, event_seq_no, timestamp_from_end) %>%
  dplyr::slice(1) %>% 
  group_by(exp_title,  implant, animal, timestamp_from_end) %>%
  dplyr::summarise(velocity.mean = mean(velocity),
                  velocity.sem = sem(velocity)) %>%
  ggplot(aes(x=timestamp_from_end/1000, group=animal)) +
  geom_ribbon(aes(ymin=velocity.mean-velocity.sem,
                 ymax=pmin(20,velocity.mean+velocity.sem)),
             alpha=0.3) +
  geom_line(aes(y=velocity.mean), color='#666666', size=0.2) +
  geom_vline(xintercept=0, linetype='dashed', color='#333333') +
  geom_hline(yintercept=4, linetype='dashed', color='#333333') +
  facet_grid(implant ~ .) +
  xlab('Time (s)') + ylab('Speed (cm/s)') +
  ylim(c(0,20)) +
  gtheme
figure.ggsave('Figure_S3D-speed-reward-zone.pdf', width=4.5, height=4.5)


# Beforetest velocity at reward proximal vs distal locations
beforetest.velocity.proximal = beforetest.aligned.pop.traces %>%
  filter(velocity <= 50) %>% # Skip incorrectly high velocity
  #filter(timestamp_from_end >= -5000, timestamp_from_start >= 0, timestamp_from_end < 0) %>%
  mutate(min.rew.dist = pmin(dist_reward0, dist_reward1),
         is.proximal = min.rew.dist <= goal.cell.max.dist) %>%
  group_by(implant, animal, day_desc, is.proximal) %>%
  dplyr::summarise(velocity.mean = mean(velocity),
                   velocity.sem = sem(velocity), .groups='drop') 

beforetest.velocity.proximal %>%
  ggplot(aes(x=implant, y=velocity.mean)) +
  geom_jitter(shape=1, width=0.2, height=0, color='#666666') +
  stat_summary(fun.data='mean_cl_boot', geom='crossbar', color='#333333') +
  facet_grid(is.proximal ~ .) +
  ylim(c(0,15)) +
  xlab('') + ylab('Mean speed (cm/s)') +
  gtheme
figure.ggsave('Figure_S3D-speed-proximal-vs-distal.pdf', width=3.8, height=4.5)

lmer.test.print(beforetest.velocity.proximal,
                var=velocity.mean,
                randef.str = '(1 | animal)',
                fixed.effects = implant * is.proximal,
                diagnostics.groupvar = is.proximal)

models = create.bayes.lm.pair(beforetest.velocity.proximal,
                              formula.full = velocity.mean ~ 1 + is.proximal * implant + animal,
                              formula.null = velocity.mean ~ 1 + is.proximal + implant + animal,
                              whichRandom = 'animal')
models$full / models$null
calc.pair.95CI(models$full,
               pair.vars = c("is.proximal:implant-dCA1", "is.proximal:implant-vCA1"),
               show.percent.change = FALSE)
