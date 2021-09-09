library(dplyr)
library(data.table)
library(datatrace)
library(stringr)

source('behav_analysis.R')
source('traces_input.R')
source('locations.R')
source('plotting_params.R')
source('fixed_effects.R')

timebin.dur.msec = 200
plot.sequences = FALSE

prepare.binned.traces = function(caimg_result_dir, 
                                 timebin.dur.msec=200,
                                 # running speed avg > 4 cm/s in 0.5 s window
                                 mean.run.velocity=3.3) {
  data.traces = read.data.trace(caimg_result_dir)
  data.traces = add.running.col(data.traces, mean.run.velocity, 10)
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
