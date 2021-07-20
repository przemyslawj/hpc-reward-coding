library(pracma)
library(purrr)
library(zoo)


mobility.bouts.tibble = function(timestamp, mobility.vals, atReward, velocity, dist_reward,
                                 window.len=3, min.bout.len=2,
                                 mobility.threshold=0.5) {
  mobility.vals[is.na(mobility.vals)] = 0
  if (window.len > 1) {
    movavg_running = zoo::rollmean(mobility.vals, window.len, fill=c(0,mobility.threshold,0))
  } else {
    movavg_running = mobility.vals
  }
  index.running.change = diff(c(0, movavg_running >= mobility.threshold))
  index.starts = which(index.running.change > 0)
  index.ends = which(index.running.change < 0) - 1
  if (length(index.starts) > length(index.ends)) {
    index.ends = c(index.ends, length(mobility.vals))
  }
  rew.at.end = map_int(index.ends, ~ max(atReward[max(1, .x - window.len) : 
                                                    min(.x + window.len, length(atReward))]) %>% 
                         as.integer)
  mean.velocity = map2_dbl(index.starts, index.ends, ~ mean(velocity[.x:.y]))
  durs = index.ends - index.starts
  valid.bouts.index = which(durs >= min.bout.len)
  nbouts = length(valid.bouts.index)
  event_ids = 1:nbouts
  if (nbouts == 0) {
    event_ids = integer(0)
  }
  tibble(index_start=index.starts[valid.bouts.index],
         index_end=index.ends[valid.bouts.index],
         timestamp_start=timestamp[index_start],
         timestamp_end=timestamp[index_end],
         mean_velocity=mean.velocity[valid.bouts.index],
         rew_at_end=rew.at.end[valid.bouts.index],
         rew_dist_at_end=dist_reward[index.ends[valid.bouts.index]],
         event_id=event_ids)
}


reward.approaches.tibble = function(timestamp, dist.reward0, dist.reward1, running.vals, velocity, 
                                    window.len=3, min.bout.len=3, min.dist.run=10,
                                    min.approach.dist=20) {
  if (window.len > 1) {
    movavg.dist.reward0 = zoo::rollmean(dist.reward0, window.len, fill=c(0,50,0))
    movavg.dist.reward1 = zoo::rollmean(dist.reward1, window.len, fill=c(0,50,0))
    movavg.running.vals = zoo::rollmean(running.vals, window.len, fill=c(0,0,0))
  } else {
    movavg.dist.reward0 = dist.reward0
    movavg.dist.reward1 = dist.reward1
    movavg.running.vals = running.vals
  }
  min.dist.reward = pmin(movavg.dist.reward0, movavg.dist.reward1)
  rew.index = map2_dbl(dist.reward0, dist.reward1, ~ ifelse(.x < .y, 0, 1))
  movavg.dist.reward0[which(movavg.running.vals < 0.5)] = 100
  movavg.dist.reward1[which(movavg.running.vals < 0.5)] = 100
  
  peaks0.mat = pracma::findpeaks(-movavg.dist.reward0, nups=window.len, ndowns=0, minpeakheight=-min.approach.dist, 
                                 minpeakdistance=min.bout.len)
  peaks1.mat = pracma::findpeaks(-movavg.dist.reward1, nups=window.len, ndowns=0, minpeakheight=-min.approach.dist, 
                                 minpeakdistance=min.bout.len)
  index.ends = c(peaks0.mat[,2], peaks1.mat[,2])
  if (is.null(index.ends)) {
    index.ends = c(length(timestamp))
  }
  index.running.change = diff(c(0, movavg.running.vals >= 0.5))
  index.running.starts = as.integer(c(which(index.running.change > 0), index.ends + 1))
  index.running.starts = sort(index.running.starts)
  index.starts = map_int(index.ends, ~ index.running.starts[max(which(index.running.starts < .x))])
  index.ends = index.ends[!is.na(index.starts)]
  index.starts = index.starts[!is.na(index.starts)]
  
  durs = index.ends - index.starts + 1
  rew.dist.run = map2_dbl(index.starts, index.ends, ~ sum(abs(movavg.dist.reward0[.y] - movavg.dist.reward0[.x]) +
                                                            abs(movavg.dist.reward1[.y] - movavg.dist.reward1[.x])))
  valid.bouts.index = which(durs >= min.bout.len & rew.dist.run >= min.dist.run)
  mean.velocity = map2_dbl(index.starts, index.ends, ~ mean(velocity[.x:.y]))
  nbouts = length(valid.bouts.index)
  event_ids = 1:nbouts
  if (nbouts == 0) {
    event_ids = integer(0)
  }
  
  tibble(index_start=index.starts[valid.bouts.index],
         index_end=index.ends[valid.bouts.index],
         timestamp_start=timestamp[index_start],
         timestamp_end=timestamp[index_end],
         mean_velocity=mean.velocity[valid.bouts.index],
         rew_dist_at_end=map_dbl(index_end, ~ ifelse(rew.index[.x] == 0, dist.reward0[.x], dist.reward1[.x])),
         rew_dist_run=rew.dist.run[valid.bouts.index],
         rew_at_end=rew.index[index_end] + 1,
         event_id=event_ids)
}

center.timestamps.around.events = function(cell.data.trace, events.tibble, exp_title_name,
                                           padding.after.event=0) {
  events.tibble = as.data.table(events.tibble)
  trial.events = events.tibble[trial==cell.data.trace$trial[1] & exp_title==exp_title_name,]
  if ('cell_id' %in% colnames(trial.events) & length(unique(trial.events$cell_id)) > 0) {
    trial.events = trial.events[cell_id == cell.data.trace$cell_id[1]]
  }

  cell.data.trace$aligned_event_id = -1
  cell.data.trace$timestamp_from_start = -1
  cell.data.trace$timestamp_from_end = -1
  cell.data.trace$rew_at_end = -1
  cell.data.trace$dist_approached_rew = -1
  if (nrow(trial.events) == 0) {
    return(cell.data.trace)
  }
  
  prev.index = 0
  for (i in 1:nrow(trial.events)) {
    timestamp.before.end = which(cell.data.trace$timestamp <= trial.events$timestamp_end)
    index.end = cell.data.trace$timestamp[timestamp.before.end[length(timestamp.before.end)] ]
    event.index.range = (prev.index+1) : min(index.end + padding.after.event, nrow(cell.data.trace))
    if (event.index.range[1] > nrow(cell.data.trace)) {
      break
    }
    cell.data.trace$aligned_event_id[event.index.range] = trial.events$event_id[i]
    cell.data.trace$timestamp_from_start[event.index.range] = 
        cell.data.trace$timestamp[event.index.range] - trial.events$timestamp_start[i]
    cell.data.trace$timestamp_from_end[event.index.range] = 
        cell.data.trace$timestamp[event.index.range] - trial.events$timestamp_end[i]
    cell.data.trace$rew_at_end[event.index.range] = trial.events$rew_at_end[i]
    prev.index = max(event.index.range)
  }
  cell.data.trace[, dist_approached_rew := ifelse(rew_at_end==1, dist_reward0, ifelse(rew_at_end == 2, dist_reward1, -1))]
  cell.data.trace
}

