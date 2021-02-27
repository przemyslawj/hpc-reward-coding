library(dplyr)
library(data.table)
library(datatrace)
library(stringr)

source('behav_analysis.R')
source('traces_input.R')
source('locations.R')
source('plotting_params.R')
source('correlations.R')
source('fixed_effects.R')

timebin.dur.msec = 200
plot.sequences = FALSE

prepare.binned.traces = function(caimg_result_dir, timebin.dur.msec=200) {
  data.traces = read.data.trace(caimg_result_dir)
  # running speed avg > 4 cm/s in 0.5 s window
  data.traces = add.running.col(data.traces, 3.3, 10)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  setorder(data.traces, exp_title, trial_id, cell_id, timestamp)
  detect.events(data.traces, deconv.threshold=0.2)

  data.traces = gauss.smooth.df.var(data.traces, filter.len=20, sigma=5.0)

  data.traces[, moved_distance := cumsum(c(0, norm2(diff(x), diff(y)))), 
              by = .(animal, date, exp_title, trial_id, trial, cell_id)]
  
  binned.traces = timebin.traces(data.traces[x > 0 & y > 0 & exp_title %in% c('trial', 'beforetest'), ],
                                 timebin.dur.msec=timebin.dur.msec) 
  binned.traces[, trial_time_bin := time_bin - min(time_bin), by=.(exp_title, trial)]
  binned.traces[, `:=` (zscored_deconv_trace = zscore(deconv_trace), 
                        zscored_trace = zscore(trace),
                        zscored_smooth_deconv_trace = zscore(smoothed_deconv_trace)),
                by=.(exp_title, cell_id)]
  
  binned.traces$atReward = as.integer(binned.traces$atReward0 >= 0.5) + 
    2 * as.integer(binned.traces$atReward1 >= 0.5)
  
  return(binned.traces)
}


get.behav.traces = function(binned.traces) {
    binned.traces[cell_id == binned.traces$cell_id[1], 
                 .(date, animal, trial, trial_id, exp_title, 
                   timestamp, time_bin, trial_time_bin, 
                   x, y, angle, 
                   atReward, atReward0, atReward1, 
                   is_headdip, is_running, velocity, smooth_velocity, 
                   dist_reward0, dist_reward1, moved_distance)]
}


get.mobility.bouts = function(behav.traces, timebin.dur.msec=200, mobility.expr=is_running,
                              min.bout.len.msec=3000) {
  mobility.expr = enexpr(mobility.expr)
  group_by(behav.traces, date, animal, exp_title, trial, trial_id) %>%
    dplyr::summarise(mobility.bouts.tibble(timestamp, !!mobility.expr, atReward, velocity, 
                                           pmin(dist_reward0, dist_reward1),
                                           window.len = ceiling(500/timebin.dur.msec),
                                           min.bout.len=min.bout.len.msec/timebin.dur.msec),
                     .groups='drop') %>%
    as.data.table()
}


get.reward.approaches = function(behav.traces, timebin.dur.msec=200, min.bout.len.msec=3000) {
  min.bout.len = min.bout.len.msec / timebin.dur.msec
  group_by(behav.traces, date, animal, exp_title, trial, trial_id) %>%
    dplyr::summarise(reward.approaches.tibble(timestamp, dist_reward0, dist_reward1, is_running, velocity,
                                              window.len = 1000 / timebin.dur.msec,
                                              min.bout.len=min.bout.len), .groups='drop') %>%
    arrange(date, animal, exp_title, trial, trial_id, timestamp_start) %>%
    as.data.table()
}


get.event.aligned.pop.traces = function(binned.traces, event.tibble, cell.active.thr=0.5) {
  pop.traces = binned.traces[, .(pop.fr=sum(nevents, na.rm=TRUE) / .N, 
                                 nevents=sum(nevents, na.rm=TRUE),
                                 ncells = .N,
                                 cells.active.pct = sum(zscored_smooth_deconv_trace >= cell.active.thr) / .N * 100,
                                 trace.mean=mean(trace, na.rm=TRUE),
                                 zscored_trace.mean=mean(zscored_trace, na.rm=TRUE),
                                 deconv_trace.mean=mean(smoothed_deconv_trace, na.rm=TRUE),
                                 zscored_deconv_trace.mean=mean(zscored_deconv_trace, na.rm=TRUE),
                                 zscored_smooth_deconv_trace.mean=mean(zscored_smooth_deconv_trace, na.rm=TRUE)
                                 ),
                             by=c('animal', 'date', 'exp_title', 'trial_id', 'trial', 
                                  'timestamp', 'abs_timestamp', 'time_bin', 'x', 'y', 
                                  'velocity', 'arrivedAtReward', 
                                  'dist_reward0', 'dist_reward1',
                                  'atReward', 'atReward0', 'atReward1', 
                                  'is_headdip', 'is_running')] 
  event.aligned.pop.traces = pop.traces[, center.timestamps.around.events(.SD, event.tibble, exp_title[1], 
                                                                          padding.after.event=2000/timebin.dur.msec), 
                                        by=.(exp_title, trial_id)][aligned_event_id >= 0,]
  event.aligned.pop.traces$event_seq_no = as.factor(
    paste(event.aligned.pop.traces$trial_id, event.aligned.pop.traces$aligned_event_id, sep='_')) %>%
    as.integer()
  return(event.aligned.pop.traces)
}


max.omit.na = function(x) {
  max(x, na.rm=TRUE)
}

order.cells.by.max.val = function(trial.traces, xvar = 'time_bin',
                                  value.varname='zscored_smooth_deconv_trace', active.val.threshold=1) {
  trial.traces$cell_id = as.factor(trial.traces$cell_id)
  #trial_trace.matrix = reshape2::acast(trial.traces, time_bin ~ cell_id, value.var=value.varname)
  trial_trace.matrix = reshape2::acast(trial.traces, as.formula(paste(xvar, '~ cell_id')), value.var=value.varname)
  cell.max.time.i = apply(trial_trace.matrix, 2, which.max)
  cell.max.val = apply(trial_trace.matrix, 2, max.omit.na)
  cell.order.df = data.frame(cell_id=as.integer(colnames(trial_trace.matrix)), 
                             max.time.i=cell.max.time.i,
                             cell.max.val=cell.max.val) %>%                
    dplyr::mutate(is.active=cell.max.val > active.val.threshold) %>%
    dplyr::arrange(desc(is.active), cell.max.time.i) %>%                                                       
    dplyr::mutate(cell_order = as.integer(dplyr::row_number())) %>%
    as.data.table()
  
  return(cell.order.df)
}

plot.activity.raster = function(df, fill.var=zscored_smooth_deconv_trace, 
                                xvar=timestamp/1000, 
                                tile.width=timebin.dur.msec/1000) {
  fill.var = enquo(fill.var)
  xvar = enexpr(xvar)
  ggplot(df) +
    geom_tile(aes(x=!!xvar, 
                  y=as.integer(cell_order), 
                  fill=!!fill.var),
              width=tile.width) +
    xlab('Time (s)') +
    ylab('Cell') +
    gtheme +
    theme(legend.position = 'none') +
    scale_fill_gradient(low='white', high = 'black') +
    scale_color_manual(values=c('#0098ffff', 'white')) +
    scale_y_reverse()
}


plot.approaches = function(aligned.trial.traces, trial_id, return.plots=FALSE, 
                           cell.active.thr=0.5, save.plots=FALSE, include.summary.plots=TRUE) {
  g.approaches = list()
  i = 1
  for (event_id in sort(unique(aligned.trial.traces$aligned_event_id))) {
    approach.traces = aligned.trial.traces[aligned_event_id==event_id & timestamp_from_start >= 0 & timestamp_from_end <= 3000,]
    behav.approach.trace = approach.traces[cell_id == approach.traces$cell_id[1],]
    approach.pop.traces = approach.traces[, .(pop.fr=sum(nevents, na.rm=TRUE) / .N, 
                                              nevents=sum(nevents, na.rm=TRUE),
                                              cells.active.pct = sum(zscored_smooth_deconv_trace >= cell.active.thr) / .N * 100,
                                              trace.mean=mean(trace, na.rm=TRUE),
                                              zscored_trace.mean=mean(zscored_trace, na.rm=TRUE),
                                              deconv_trace.mean=mean(smoothed_deconv_trace, na.rm=TRUE),
                                              zscored_deconv_trace.mean=mean(zscored_deconv_trace, na.rm=TRUE),
                                              zscored_smooth_deconv_trace.mean=mean(zscored_smooth_deconv_trace, na.rm=TRUE)),
                                          by=c('timestamp_from_end')]
    
    cell.order.df = order.cells.by.max.val(approach.traces, active.val.threshold=1.0)
    df.joined = approach.traces[cell.order.df, on=.(cell_id)]
    g.activity = plot.activity.raster(df.joined, fill.var=zscored_smooth_deconv_trace, 
                                      xvar=timestamp_from_end/1000,
                                      tile.width=timebin.dur.msec/1000)
    if (include.summary.plots) {
      g.behav = ggplot(behav.approach.trace) +
        geom_line(aes(x=timestamp_from_end/1000, y=dist_approached_rew)) +
        xlab('Time (s)') + ylab('Reward dist') + gtheme
      g.pop.trace = ggplot(approach.pop.traces) + geom_line(aes(x=timestamp_from_end/1000, y=zscored_smooth_deconv_trace.mean)) + 
        xlab('Time (s)') + ylab('df/f') + gtheme
      g.active.cells = ggplot(approach.pop.traces) + geom_line(aes(x=timestamp_from_end/1000, y=cells.active.pct)) +
        xlab('Time (s)') + ylab('Active %') + gtheme
      g.approaches[[i]] = cowplot::plot_grid(g.activity, g.behav, g.pop.trace, g.active.cells, ncol=1, rel_heights = c(2,1,1,1))
    } else {
      g.approaches[[i]] = g.activity
    }
    i = i + 1
  }
  
  napproaches = length(g.approaches)
  if (save.plots & napproaches > 0) {
    nrows = ceil(napproaches/3)
    fname = paste0(trial_id, '.png')
    ggsave(fname, 
           cowplot::plot_grid(plotlist=g.approaches, ncol = 3), 
           path='/home/prez/tmp/cheeseboard/approaches',
           units='cm', height=nrows*11, width=15)
  }
  if (return.plots) {
    return(g.approaches)
  }
  return(napproaches)
}

plot.approaches.dist = function(binned.dist.traces, trial_id, dist.bin.width=5) {
  g.approaches = list()
  i = 1
  for (event_id in unique(binned.dist.traces$aligned_event_id)) {
    approach.traces = binned.dist.traces[aligned_event_id==event_id]
    cell.order.df = order.cells.by.max.val(approach.traces, xvar='bin.dist_approached_rew', active.val.threshold=1.0)
    df.joined = approach.traces[cell.order.df, on=.(cell_id)]
    g.approaches[[i]] = plot.activity.raster(df.joined,
                                             fill.var=zscored_smooth_deconv_trace, 
                                             xvar=dist_approached_rew,
                                             tile.width=dist.bin.width) + 
      xlab('Distance to reward')
    i = i + 1
  }
  
  napproaches = length(g.approaches)
  if (napproaches > 0) {
    nrows = ceil(napproaches/3)
    fname = paste0(trial_id, '.png')
    ggsave(fname, 
           cowplot::plot_grid(plotlist=g.approaches, ncol = 3), 
           path='/home/prez/tmp/cheeseboard/approaches_dist',
           units='cm', height=nrows*7, width=15)
  }
  return(g.approaches)
}


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
seq.cell.timestamp.var.list = list()
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
  seq.cell.timestamp.var.list[[i]] = day.seq.cell.timestmap.var
  
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
  
  # # Activity within population sequence vs at other times
  # rew.aligned.traces = binned.traces[, center.timestamps.around.events(.SD, reward.arriving.bouts, exp_title[1]), 
  #                                    by=.(exp_title, trial_id, cell_id)][aligned_event_id >= 0,]
  # rew.aligned.traces$is.event = rew.aligned.traces$nevents
  # # TimeBin with 600 msec window -- issue when activity distributed between two bins
  # seq.max.sec = 5
  # window.dur.msec = 1000
  # rew.aligned.traces = rew.aligned.traces[timestamp_from_end <= 0.0 & timestamp_from_end >= -seq.max.sec * 1000.0,
  #                                         `:=` (roll.nevents = zoo::rollmean(as.double(nevents), 
  #                                                                            round(window.dur.msec / timebin.dur.msec / 2), 
  #                                                                            fill=c(0,0,0), align='center')),
  #                                         by=.(animal, date, exp_title, trial_id, trial, cell_id, aligned_event_id)][
  #                                           , timestamp_to_max := timestamp_from_start - timestamp_from_start[which.max(roll.nevents)],
  #                                           by=.(animal, date, exp_title, trial_id, trial, cell_id, aligned_event_id)]
  # 
  # rew.aligned.traces[timestamp_from_end <= 0.0 & timestamp_from_end >= -seq.max.sec * 1000.0, 
  #                    .(nevents_inside = as.integer(sum(ifelse(abs(timestamp_to_max) <= window.dur.msec / 2, nevents, 0))),
  #                      nevents = sum(nevents),
  #                      bout_dur = max(timestamp) - min(timestamp)),
  #                    by=.(exp_title, animal, date, trial_id, cell_id, aligned_event_id)][
  #                      bout_dur >= 4000 & nevents_inside > 0,
  #                      perc_inside := nevents_inside / nevents] ->
  #   cell.fr.outside.seq
  # 
  # all.cells.fr.outside.seq.list[[i]] = cell.fr.outside.seq
  
#   running.trace.matrix = reshape2::acast(binned.traces[exp_title == 'trial' & is_running > 0], 
#                                          time_bin ~ cell_id, value.var='zscored_trace')
#   immobility.trace.matrix = reshape2::acast(binned.traces[exp_title == 'trial' & is_running == 0], 
#                                             time_bin ~ cell_id, value.var='zscored_trace')
#   if (dim(running.trace.matrix)[2] > 1) {
#     day.cors.df = bind_rows(
#       dplyr::mutate(create.cor.values.df(cor(running.trace.matrix)), running=TRUE),
#       dplyr::mutate(create.cor.values.df(cor(immobility.trace.matrix)), running=FALSE)) %>% 
#       mutate(animal = binned.traces$animal[1],
#              date = binned.traces$date[1])
#     all.cors.df.list[i] = day.cors.df
#   }
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
seq.cell.timestamp.var = join.meta.df(rbindlist(seq.cell.timestamp.var.list))
#all.cells.fr.outside.seq = join.meta.df(rbindlist(all.cells.fr.outside.seq.list))

save.image(file="data/rew_aligned_pop_activity10.RData")

early.learning.days = c('learning1 day#1')
late.learning.days = c('learning1 day#5', 'learning2 day#2', 'learning3 day#2', 'learning4 day#2')

reward.aligned.pop.traces.earlylate = reward.aligned.pop.traces %>%
  filter(exp_title == 'trial', aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
                is.late.learning = day_desc %in% late.learning.days) %>%
  dplyr::filter(is.early.learning | is.late.learning) %>%
  dplyr::mutate(bout=ifelse(is.early.learning, 'early', 'late'))

reward.aligned.pop.traces.earlylate.by.pc = reward.aligned.pop.traces.by.pc %>%
  filter(exp_title == 'trial', aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
                is.late.learning = day_desc %in% late.learning.days) %>%
  dplyr::filter(is.early.learning | is.late.learning) %>%
  dplyr::mutate(bout=ifelse(is.early.learning, 'early', 'late'))

# No ramping of activity before mobility, possibly synchronous events present not visible in the mean
immobility.aligned.pop.traces %>%
  filter(exp_title == 'trial', 
         aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
                is.late.learning = day_desc %in% late.learning.days) %>%
  dplyr::filter(is.early.learning | is.late.learning) %>%
  group_by(implant, is.early.learning, timestamp_from_end) %>%
  dplyr::summarise(
    zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    zscored_trace.sem=sem(zscored_trace.mean),
    zscored_trace.mean=mean(zscored_trace.mean)) %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=is.early.learning)) +
  geom_ribbon(aes(ymin=zscored_trace.mean-zscored_trace.sem,
                  ymax=zscored_trace.mean+zscored_trace.sem,
                  fill=is.early.learning), alpha=0.5) +
  scale_fill_manual(values=three.colours) +
  facet_wrap(. ~ implant) +
  geom_vline(xintercept = 0, linetype='dashed') +
  labs(title='Population activity before movement onset') +
  xlab('Time (s)') +
  xlim(c(-5, 2)) +
  ylim(c(-0.15, 0.2)) +
  ylab('z-scored deconv signal') +
  gtheme

nonreward.aligned.pop.traces.filtered = filter(nonreward.aligned.pop.traces, exp_title == 'trial', 
       aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  dplyr::mutate(bout='non-reward')


## Beforetest trials
find.peak.vals = function(vals, timestamp_from_end) {
  peak.index = which.max(vals)
  tibble(max.val=vals[peak.index], timestamp_from_end=timestamp_from_end[peak.index])
}

beforetest.aligned.pop.traces %>%
  filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000, timestamp_from_end <= 2000) %>%
  group_by(implant, animal, date, aligned_event_id) %>%
  dplyr::summarise(find.peak.vals(#zscored_smooth_deconv_trace.mean, 
                                  cells.active.pct,
                                  timestamp_from_end)) -> 
  time.to.peak

# time.to.peak %>%
#   #filter(animal=='F-BL', date=='2019-09-08') %>%
#   ggplot(aes(x=timestamp_from_end/1000, y=max.val, text=paste(animal, date, aligned_event_id))) +
#   geom_point() +
#   facet_grid(. ~ implant) +
#   gtheme -> g
# plotly::ggplotly(g, tooltip='text')
# 
# time.to.peak %>%
#   ggplot(aes(x=timestamp_from_end/1000)) +
#   geom_histogram(binwidth=1) +
#   facet_wrap(. ~ implant + animal) +
#   gtheme

# # Mean pop activity before reward and during approach (split into two 4-sec bins)
# beforetest.pop.zone.summary = beforetest.aligned.pop.traces %>%
#   filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000, timestamp_from_end <= 2000) %>%
#   mutate(reward.approach.zone = ifelse(timestamp_from_end < -2000, 'before', 'inside')) %>%
#   mutate(comparison_id = paste(animal, date, aligned_event_id)) %>%
#   group_by(implant, animal, date, aligned_event_id, comparison_id, reward.approach.zone) %>%
#   dplyr::summarise(zscored.activity.mean = mean(zscored_smooth_deconv_trace.mean), 
#                    zscored.activity.max = max(zscored_smooth_deconv_trace.mean),
#                    active.cells.mean = mean(cells.active.pct),
#                    active.cells.max = max(cells.active.pct))
# 
# beforetest.pop.zone.summary %>%
#   filter(implant == 'dCA1') %>%
#   ggplot() +
#   geom_line(aes(x=reward.approach.zone, y=active.cells.max, group=comparison_id)) +
#   facet_wrap(. ~ implant + animal) +
#   gtheme
# 
# pop.zone.peak.model = lmerTest::lmer(active.cells.max ~ reward.approach.zone + (1 + date | animal) + (1 | comparison_id),
#                                      beforetest.pop.zone.summary,
#                                      REML=TRUE, subset=(implant=='dCA1'))
# plot.model.diagnostics(pop.zone.peak.model,
#                        subset(beforetest.pop.zone.summary, (implant=='dCA1'))$animal,
#                        subset(beforetest.pop.zone.summary, (implant=='dCA1'))$reward.approach.zone)
# summary(pop.zone.peak.model)
# anova(pop.zone.peak.model, refit=FALSE, ddf='Satterthwaite')

beforetest.aligned.pop.summary.by.pc = beforetest.aligned.pop.traces.by.pc %>%
  filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000, timestamp_from_end <= 2000) %>%
  group_by(implant, is.pc, timestamp_from_end) %>%
  dplyr::summarise(
    zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    zscored_trace.sem=sem(zscored_smooth_deconv_trace.mean),
    zscored_trace.mean=mean(zscored_smooth_deconv_trace.mean),
    cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
    cells.active.pct.sem=sem(cells.active.pct), .groups='drop') %>% 
  as.data.table()


beforetest.aligned.pop.summary.by.pc %>%
  ggplot(aes(x=timestamp_from_end / 1000)) +
  #geom_ribbon(aes(ymin=zscored_deconv_trace.mean-zscored_deconv_trace.sem, 
  #                ymax=zscored_deconv_trace.mean+zscored_deconv_trace.sem,
  #                group=is.pc, fill=is.pc)) +
  geom_ribbon(aes(ymin=cells.active.pct.mean-cells.active.pct.sem,
                  ymax=cells.active.pct.mean+cells.active.pct.sem,
                  group=is.pc, fill=is.pc), alpha=0.75) +
  geom_vline(xintercept=0, linetype='dashed') +
  facet_wrap(. ~ implant, scales='free') +
  scale_fill_manual(values=side.two.coulours) +
  labs(title='Approach to reward location (beforetest)') +
  ylab('Cells active (%)') +
  xlab('Time (s)') +
  gtheme +
  theme(legend.position='top')
ggsave('~/tmp/cheeseboard/pop_activity/beforetest_reward_approach.svg', 
       height=6.5, width=8, units='cm')

# Stats on bouts arriving at reward in beforetest
timebin.width = 1000
beforetest.aligned.pop.traces.by.pc.binned = beforetest.aligned.pop.traces.by.pc %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial_id, day_desc, is.pc, timestamp_from_end_bin, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
beforetest.aligned.pop.traces.by.pc.binned$proximal.timestamp = as.factor(beforetest.aligned.pop.traces.by.pc.binned$proximal.timestamp)
beforetest.aligned.pop.traces.by.pc.binned$is.pc = as.factor(beforetest.aligned.pop.traces.by.pc.binned$is.pc)

beforetest.aligned.pop.traces.by.pc.binned.per.trial = beforetest.aligned.pop.traces.by.pc %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial_id, day_desc, is.pc, timestamp_from_end_bin) %>%
  dplyr::summarise_all(mean) %>% 
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
beforetest.aligned.pop.traces.by.pc.binned.per.trial$proximal.timestamp = 
  as.factor(beforetest.aligned.pop.traces.by.pc.binned.per.trial$proximal.timestamp)
beforetest.aligned.pop.traces.by.pc.binned.per.trial$is.pc = 
  as.factor(beforetest.aligned.pop.traces.by.pc.binned.per.trial$is.pc)

implant.loc = 'dCA1'
beforetest.aligned.pop.traces.by.pc.binned %>%
  filter(implant == implant.loc) %>%
  mutate(bout_id=paste(trial_id, aligned_event_id)) %>%
  ggplot(aes(x=proximal.timestamp, y=cells.active.pct, group=trial_id)) +
  geom_jitter(size=0.05, shape=1, alpha=0.35, width=0.2, height=0, color='#444444') +
  geom_line(data=subset(beforetest.aligned.pop.traces.by.pc.binned.per.trial, implant==implant.loc),
            size=0.25, color='#333333') +
  facet_wrap(implant ~ is.pc) +
  ylim(c(0, 45)) +
  gtheme +
  xlab('Time to reward location') + scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  ylab('Active cells %')
ggsave(paste0('fig2-approach-stats-beforetest-',implant.loc, '.pdf'),
       path = '/home/prez/tmp/cheeseboard/',
       device=cairo_pdf, units='cm', width=4.5, height=4.6)

m.beforetest = lmer.test.print(subset(beforetest.aligned.pop.traces.by.pc.binned, implant == implant.loc),
                               cells.active.pct,
                               fixed.effects = is.pc * proximal.timestamp,
                               #fixed.effects = proximal.timestamp, 
                               #randef.str = '(1 | animal)',
                               diagnostics.groupvar=is.pc)
pairwise.post.hoc(m.beforetest, factor.interaction = c('is.pc:proximal.timestamp'))

models = create.bayes.lm.pair(subset(beforetest.aligned.pop.traces.by.pc.binned, implant == implant.loc & is.pc=='TRUE'),
                              formula.full = cells.active.pct ~ 1 + proximal.timestamp + animal,
                              formula.null = cells.active.pct ~ 1 + animal,
                              whichRandom = 'animal')
models$full / models$null
calc.pair.95CI(models$full, pair.vars = c('proximal.timestamp-0', 'proximal.timestamp-1'))
group_by(beforetest.aligned.pop.traces.by.pc.binned, implant, is.pc, proximal.timestamp) %>%
  dplyr::summarise(mean(cells.active.pct), sem(cells.active.pct), n())

# Activity averaged over distance - not a good comparison as locations at reward have higher occupancy, and mean activity lower
# dist.bin.width = 10
# dist.pop.traces = stimbin.traces(beforetest.aligned.pop.traces, 
#                                  dist_approached_rew, 100 / dist.bin.width, max.width = 100) %>% as.data.table()
# binned.dist.pop.traces = dist.pop.traces[timestamp_from_start >= 0 & timestamp_from_end <= 2000 & timestamp_from_end >= -5000, 
#                                         .(zscored_smooth_deconv_trace = mean(zscored_smooth_deconv_trace.mean),
#                                           zscored_deconv_trace = mean(zscored_deconv_trace.mean),
#                                           zscored_trace = mean(zscored_trace.mean),
#                                           dist_approached_rew = dist.bin.width/2 + (bin.dist_approached_rew - 1) * dist.bin.width, 
#                                           n=.N),
#                                           by=.(animal, date, aligned_event_id, bin.dist_approached_rew)]
#                                           #by=.(animal, date, bin.dist_approached_rew)]
#                                           #by=.(animal, bin.dist_approached_rew)]
# 
# ggplot(binned.dist.pop.traces, aes(x=dist_approached_rew)) +
#   geom_ribbon(aes(ymin=zscored_smooth_deconv_trace.mean-zscored_smooth_deconv_trace.sem,
#                   ymax=zscored_smooth_deconv_trace.mean+zscored_smooth_deconv_trace.sem)) +
#   #geom_line(aes(y=zscored_smooth_deconv_trace.mean, group=aligned_event_id)) +
#   facet_wrap(. ~ animal + date, scales='free') +
#   geom_hline(yintercept = 0, linetype='dashed') +
#   xlim(c(0, 50)) +
#   labs(title='Population activity when approaching rewarded location (beforetest)') +
#   ylab('z-scored deconv signal') +
#   xlab('Distance (%)') +
#   gtheme
  
reward.pop.traces.aligned.comparison = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>% 
  group_by(implant, timestamp_from_end, bout) %>%
  dplyr::summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
                   zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
                   zscored_trace.sem=sem(zscored_trace.mean),
                   zscored_trace.mean=mean(zscored_trace.mean),
                   cells.active.pct.mean=mean(cells.active.pct, na.rm=TRUE),
                   cells.active.pct.sem=sem(cells.active.pct), .groups='drop') 

reward.pop.traces.aligned.comparison %>%
  ggplot(aes(x=timestamp_from_end / 1000, group=bout)) +
  geom_ribbon(aes(ymin=zscored_deconv_trace.mean-zscored_deconv_trace.sem, 
                  ymax=zscored_deconv_trace.mean+zscored_deconv_trace.sem,
                  fill=bout), alpha=0.6) +
  facet_wrap(. ~ implant, scales='free') +
  scale_fill_manual(values=three.colours) +
  scale_color_manual(values=three.colours) +
  #scale_alpha_manual(values=c(0.25, 0.5, 0.25)) +
  scale_alpha_manual(values=c(0.6, 0.6, 0.6)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  labs(title='Learning-changed activity reward approach') +
  xlab('Time (s)') +
  xlim(c(-5, 2)) +
  #ylim(c(-0.15, 0.4)) +
  ylab('deconvolved dF/F (z-score)') +
  gtheme +
  theme(legend.position = 'top')
ggsave('~/tmp/cheeseboard/pop_activity/reward_approach.svg', 
       height=6.5, width=8, units='cm')


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
ggsave('~/tmp/cheeseboard/pop_activity/learnt_reward_approach.svg', 
       height=6.5, width=8, units='cm')

timebin.width = 1000
reward.aligned.pop.traces.late.by.pc.binned.per.day = filter(reward.aligned.pop.traces.earlylate.by.pc, is.late.learning) %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, day_desc, is.pc, proximal.timestamp, timestamp_from_end_bin) %>%
  dplyr::summarise_all(mean) %>% 
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2,
                aday=paste(animal, day_desc, sep='_'))
reward.aligned.pop.traces.late.by.pc.binned.per.day$is.pc = 
  as.factor(reward.aligned.pop.traces.late.by.pc.binned.per.day$is.pc)
reward.aligned.pop.traces.late.by.pc.binned.per.day$proximal.timestamp = 
  as.factor(reward.aligned.pop.traces.late.by.pc.binned.per.day$proximal.timestamp)

reward.aligned.pop.traces.late.by.pc.binned = filter(reward.aligned.pop.traces.earlylate.by.pc, is.late.learning) %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial, trial_id, day_desc, is.pc, 
                  proximal.timestamp, timestamp_from_end_bin, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>% 
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
reward.aligned.pop.traces.late.by.pc.binned$is.pc = as.factor(reward.aligned.pop.traces.late.by.pc.binned$is.pc)
reward.aligned.pop.traces.late.by.pc.binned$proximal.timestamp = as.factor(reward.aligned.pop.traces.late.by.pc.binned$proximal.timestamp)

implant.loc = 'vCA1'
reward.aligned.pop.traces.late.by.pc.binned %>%
  filter(implant == implant.loc) %>%
  ggplot(aes(x=proximal.timestamp, y=cells.active.pct)) +
  geom_jitter(size=0.05, shape=1, alpha=0.35, width=0.25, height=0, color='#444444') +
  geom_line(data=subset(reward.aligned.pop.traces.late.by.pc.binned.per.day, implant==implant.loc),
            mapping=aes(group=aday),
            size=0.25, color='#333333') +
  facet_wrap(implant ~ is.pc) +
  ylim(c(0, 55)) +
  gtheme +
  xlab('Time to reward location') + 
  scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  ylab('Active cells %')
ggsave(paste0('fig2-approach-stats-latelearning-',implant.loc, '.pdf'),
       path = '/home/prez/tmp/cheeseboard/',
       device=cairo_pdf, units='cm', width=4.5, height=4.6)

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

# Stats on bouts at reward
timebin.width = 1000
reward.aligned.pop.traces.earlylate.binned = bind_rows(
    reward.aligned.pop.traces.earlylate,
    filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, trial, trial_id, day_desc, 
                  proximal.timestamp, bout, aligned_event_id) %>%
  dplyr::summarise_all(mean) %>% 
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2)
reward.aligned.pop.traces.earlylate.binned$bout = as.factor(reward.aligned.pop.traces.earlylate.binned$bout)
reward.aligned.pop.traces.earlylate.binned$proximal.timestamp = as.factor(reward.aligned.pop.traces.earlylate.binned$proximal.timestamp)

reward.aligned.pop.traces.earlylate.binned.per.day = bind_rows(
  reward.aligned.pop.traces.earlylate,
  filter(nonreward.aligned.pop.traces.filtered, day_desc %in% c(early.learning.days, late.learning.days))) %>%
  dplyr::mutate(distal.timestamp = (timestamp_from_end <= -4000) & (timestamp_from_end > -5000),
                proximal.timestamp = (timestamp_from_end > -1000) & (timestamp_from_end <= 0)) %>%
  filter(distal.timestamp | proximal.timestamp) %>%
  dplyr::mutate(timestamp_from_end_bin = ceiling(timestamp_from_end / timebin.width)) %>%
  dplyr::group_by(implant, date, animal, exp, exp_title, day_desc, proximal.timestamp, bout) %>%
  dplyr::summarise_all(mean) %>% 
  dplyr::mutate(timestamp_from_end = timestamp_from_end_bin * timebin.width - timebin.width/2,
                aday=paste(animal, day_desc, sep='_'))
reward.aligned.pop.traces.earlylate.binned.per.day$bout = 
  as.factor(reward.aligned.pop.traces.earlylate.binned.per.day$bout)
reward.aligned.pop.traces.earlylate.binned.per.day$proximal.timestamp = 
  as.factor(reward.aligned.pop.traces.earlylate.binned.per.day$proximal.timestamp)

implant.loc = 'dCA1'
reward.aligned.pop.traces.earlylate.binned %>%
  filter(implant == implant.loc) %>%
  ggplot(aes(x=proximal.timestamp, y=zscored_deconv_trace.mean)) +
  geom_jitter(size=0.05, shape=1, alpha=0.3, width=0.3, height=0, color='#444444') +
  geom_line(data=subset(reward.aligned.pop.traces.earlylate.binned.per.day, implant==implant.loc),
            mapping=aes(group=aday),
            size=0.25, color='#333333') +
  facet_wrap(implant ~ bout) +
  ylim(c(-0.1, 0.85)) +
  #ylim(c(-0.3, 0.65)) +
  gtheme +
  xlab('Time to reward location') + 
  scale_x_discrete(labels=c('4-5 s', '0-1 s')) +
  #scale_y_log10() +
  ylab('zscored dF')
ggsave(paste0('fig2-approach-stats-activity-',implant.loc, '.pdf'),
       path = '/home/prez/tmp/cheeseboard/',
       device=cairo_pdf, units='cm', width=5.6, height=4.6)

m.earlylate = lmer.test.print(
  subset(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc),
  log(1.0 + zscored_deconv_trace.mean),
  fixed.effects = bout * proximal.timestamp,
  diagnostics.groupvar=bout)
group_by(subset(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc), bout, proximal.timestamp) %>%                      
  dplyr::summarise(mean(zscored_deconv_trace.mean), sem(zscored_deconv_trace.mean), n())
pairwise.post.hoc(m.earlylate, factor.interaction = c('bout:proximal.timestamp'))

models = create.bayes.lm.pair(
  filter(reward.aligned.pop.traces.earlylate.binned, implant==implant.loc, bout == 'non-reward') %>% 
    mutate(val = log(1.0 + zscored_deconv_trace.mean)),
  formula.full = zscored_deconv_trace.mean ~ 1 + proximal.timestamp + animal,
  formula.null = zscored_deconv_trace.mean ~ 1 + animal,
  whichRandom = 'animal', iterations = 50000)
models$full / models$null
calc.pair.95CI(models$full, pair.vars = c('proximal.timestamp-FALSE', 'proximal.timestamp-TRUE'),
               ytransform = function(x) {exp(x) - 1 }, show.percent.change = FALSE)

day.mean.on.arrival = reward.aligned.pop.traces %>%
  filter(aligned_event_id >= 0, timestamp_from_start >= 0, timestamp_from_end >= -5000) %>%
  group_by(day_desc, exp_title, implant, animal, timestamp_from_end) %>%
  summarise(zscored_deconv_trace.sem=sem(zscored_smooth_deconv_trace.mean),
            zscored_deconv_trace.mean=mean(zscored_smooth_deconv_trace.mean),
            zscored_trace.sem=sem(zscored_trace.mean),
            zscored_trace.mean=mean(zscored_trace.mean))

day.mean.on.arrival %>%
  filter(day_desc %in% c('learning2 day#2', 'learning3 day#2', 'learning4 day#2')) %>%
  #filter(day_desc %in% c('learning1 day#1', 'learning2 day#1')) %>%
  ggplot() +
  geom_line(aes(x=timestamp_from_end / 1000, y=zscored_deconv_trace.mean, group=day_desc)) +
  facet_wrap(animal + implant ~ ., scales='free') +
  gtheme

### Plot data for a single date
#reward.aligned.pop.traces[, zscored.trace := pmin(10, zscore(smooth_deconv_trace)),
#                         by=.(date, animal, exp_title, trial_id)]
reward.aligned.pop.traces %>%
  filter(day_desc == 'learning2 day#2', implant=='dCA1') %>%
  filter(aligned_event_id >= 0, timestamp_from_start > 0) %>%
  dplyr::mutate(timestamp_from_end_sec = timestamp_from_end / 1000.0) %>% 
  ggplot() +
  geom_tile(aes(x=timestamp_from_end_sec, 
                y=event_seq_no, 
                fill=zscored_deconv_trace.mean,
                width=timebin.dur.msec / 1000, height=1)) +
  #geom_tile(aes(x=round(bin.moved_distance  * max.moved.dist / dist.nbins), y=-2, color=atReward)) +
  xlab('Time (s)') +
  ylab('Trial') +
  gtheme +
  #theme(legend.position = 'none') +
  facet_wrap(implant + animal ~ .) +
  scale_fill_gradient(low='white', high = 'black') +
  scale_color_manual(values=c('#0098ffff', 'white')) +
  xlim(c(-5,0)) +
  scale_y_reverse()

# # Lower correlations can be explained by lower FR
# all.cors.df = join.meta.df(all.cors.df)
# all.cors.df %>%
#   dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
#                 is.late.learning = day_desc %in% late.learning.days,
#                 animal_learning = paste(animal, as.character(is.early.learning), sep='_')) %>%
#   dplyr::filter(is.early.learning | is.late.learning) %>%
#   ggplot() + 
#   stat_ecdf(mapping=aes(x=vals, colour=is.early.learning, group=animal_learning), geom='step') +
#   scale_y_continuous(breaks=c(0.0, 0.5, 1.0)) +
#   scale_x_continuous(breaks=c(-0.5, 0.0, 0.5), limits = c(-0.6,0.6)) +
#   #scale_color_manual(name='', values=is.early) +
#   facet_grid(implant ~ running) +
#   xlab('Correlation') + ylab('Probability')

# What causes the increased firing? more cells firing or higher activity

## Sequential firing

# Variance in sequential firing higher for vCA1 cells
seq.cell.timestamp.var %>% 
  #filter(str_starts(day_desc, 'learning1')) %>%
  filter(day_desc == 'learning2 day#2') %>%
  ggplot() +
  geom_boxplot(aes(x=day_desc, y=sqrt(timestamp.var))) +
  #geom_point(aes(x=sqrt(timestamp.var), y=zscored_smooth_deconv_trace.peak)) +
  facet_wrap(. ~ implant) +
  gtheme

#  Variance of sequence in place cells vs other cells lower in dCA1, their peak higher
left_join(seq.cell.timestamp.var, place.cell.db, by=c('date', 'day_desc', 'animal', 'implant', 'cell_id')) %>%
  #filter(day_desc == 'learning2 day#2') %>%
  ggplot() +
  geom_boxplot(aes(x=is.pc, y=sqrt(timestamp.var))) +
  #geom_boxplot(aes(x=is.pc, y=zscored_smooth_deconv_trace.peak)) +
  facet_wrap(. ~ implant) +
  gtheme


# # Mean nevents inside seq per cell during the day
# all.cells.fr.outside.seq = as.data.table(all.cells.fr.outside.seq)
# 
# # Not exactly correct
# total.day.cells = all.cells.fr.outside.seq[, cell_id, by=.(date, animal, cell_id)][
#   , .(ncells=.N), by=.(animal, date)][, .(max.ncells=max(ncells)), by=.(animal)]
# all.cells.fr.outside.seq[, .(pct.inside = sum(nevents_inside) / sum(nevents_inside + nevents_outside) * 100,
#                              total_events = sum(nevents_inside) + sum(nevents_outside),
#                              fr = (sum(nevents_inside) + sum(nevents_outside)) / bout_dur,
#                              ncells=.N),
#                          by=.(exp_title, implant, animal, date, day_desc, trial_id, aligned_event_id)] [
#                            ,
#                            .(pct.inside.mean = mean(pct.inside),
#                              pct.inside.sem = sem(pct.inside),
#                              fr.per.bout.mean = mean(fr),
#                              fr.per.bout.sem = sem(fr),
#                              ncells.mean = mean(ncells),
#                              ncells.sem = sem(ncells),
#                              nbouts=.N),
#                            by=.(exp_title, implant, animal, date, day_desc)
#                          ] ->
# cell.mean.events.inside.seq
# cell.mean.events.inside.seq = cell.mean.events.inside.seq[total.day.cells, on=.(animal)]
# cell.mean.events.inside.seq[,
#   `:=` (pct.cells.mean = ncells.mean/max.ncells * 100,
#         pct.cells.sem = ncells.sem/max.ncells * 100)]
# 
# summary.mean.events.inside.seq = cell.mean.events.inside.seq %>%
#   dplyr::mutate(is.early.learning = day_desc %in% early.learning.days,
#                 is.late.learning = day_desc %in% late.learning.days) %>%
#   #dplyr::filter(is.early.learning | is.late.learning) %>%
#   dplyr::mutate(bout=ifelse(is.early.learning, 'early', 'late'))
# 
# summary.mean.events.inside.seq %>%
#   filter(exp_title == 'trial') %>%
#   ggplot(aes(x=day_desc, group=animal)) +
#   geom_line(aes(y=ncells.mean, color=animal)) +
#   #geom_line(aes(y=pct.inside.mean, color=animal)) +
#   #geom_ribbon(aes(ymin=pct.inside.mean-pct.inside.sem, ymax=pct.inside.mean+pct.inside.sem), 
#   #            alpha=0.1) +
#   facet_grid(implant ~ ., scales = 'free') +
#   #ylim(c(0.6,1.0)) +
#   gtheme


all.rew.arriving.bouts %>%
  group_by(implant, animal, day_desc) %>%
  dplyr::summarise(bout_sec = mean(timestamp_end - timestamp_start) / 1000) %>%
  ggplot() +
  geom_line(aes(x=day_desc, y=bout_sec, group=animal)) +
  facet_grid(implant ~ .) +
  gtheme


