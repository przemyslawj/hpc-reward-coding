get.arrived.indecies = function(at.reward) {
  reward.arrived = which(diff(at.reward)==1)
  reward.left = which(diff(at.reward)==-1)
  
  if (length(at.reward) == 0) {
    return(c())
  }
  if (at.reward[1]) {
    reward.arrived = c(1, reward.arrived)
  }
  if (length(reward.left) < length(reward.arrived)) {
    reward.left = c(reward.left, length(at.reward))
  }
  if (length(reward.left) != length(reward.arrived)) {
    print('Warning: incorrect lengths')
  }
  
  return(list(arrived=reward.arrived, left=reward.left))
}

fill.gaps = function(pos, timestamps, max.gap.ms=1000) {
  # Replaces 0s in pos vector of 01s if the gaps with 0s are shorter than max.gaps.ms 
  reward.indecies = get.arrived.indecies(pos)
  reward.arrived = reward.indecies[['arrived']]
  reward.left = reward.indecies[['left']]
  
  if (length(reward.arrived) > 1) {
    for (i in 1:(length(reward.left) - 1)) {
      left_dur = timestamps[reward.arrived[i+1]] -  timestamps[reward.left[i]]
      if (left_dur < 1000) {
        pos[reward.left[i] : reward.arrived[i+1]] = 1
      }
    }
  }
  
  return(pos)
}

get.crossing.durations = function(inside_roi, timestamps) {
  
  zone.indecies = get.arrived.indecies(inside_roi)
  zone.arrived = zone.indecies[['arrived']]
  zone.left = zone.indecies[['left']]
  crossing_dur = c()
  
  if (length(zone.left) > 0) {
    crossing_dur = rep(0, length(zone.left))
    for (i in 1:length(zone.left)) {
      crossing_dur[i] = timestamps[zone.left[i]] -  timestamps[zone.arrived[i]]
    }
  }
  
  return(crossing_dur)
}

is.at.reward = function(velocities, rew.dist, timestamps, running.thresh = 1.5, rew.dist.thresh=10) {
  min.duration.at.reward.ms = 2000
  at.reward = rep(0, length(velocities))
  res = rep(0, length(velocities) + 1)
  at.reward = velocities < running.thresh &
      rew.dist < rew.dist.thresh
  
  # trick to handle recordings finished just when arrived to reward
  at.reward = c(at.reward, 1)
  timestamps = c(timestamps, timestamps[length(timestamps)] + min.duration.at.reward.ms)
  
  at.reward = fill.gaps(at.reward, timestamps)
  
  reward.indecies = get.arrived.indecies(at.reward)
  reward.arrived = reward.indecies[['arrived']]
  reward.left = reward.indecies[['left']]
  
  # Discard too short stops at reward
  if (length(reward.left) > 0) {
    at.reward.duration.ms = timestamps[reward.left] - timestamps[reward.arrived]
    keep_indecies = which(at.reward.duration.ms > min.duration.at.reward.ms)
    for (i in keep_indecies) {
      res[reward.arrived[i]:reward.left[i]] = 1
    }
  }
  
  return(res[1:length(res) - 1])
}

# Returns logical vector with TRUE for timestamps which are on the direct approach to reward0.
# The approach is capped at max duration.
is.reward.approach = function(atReward0, atReward1, timestamps, max.approach.dur.sec=3) {
  arrived.rew0 = min(which(atReward0 > 0))
  arrived.rew1 = min(which(atReward1 > 0))
  approach.rew0.start.msec = ifelse(timestamps[arrived.rew0] < timestamps[arrived.rew1],
                                    max(0, timestamps[arrived.rew0] - max.approach.dur.sec * 1000),
                                    max(timestamps[arrived.rew1], timestamps[arrived.rew0] - max.approach.dur.sec * 1000))
  (timestamps >= approach.rew0.start.msec) & (timestamps <= timestamps[arrived.rew0])
}

reward.approach.timestamp = function(atReward, timestamps) {
  arrived.rew = min(which(atReward > 0))
  arrived.rew.msec = timestamps[arrived.rew]
  timestamps - arrived.rew.msec
}

next.goal.reward = function(is_at_rew0, is_at_rew1) {
  trial_end = length(is_at_rew0)
  next.goal = rep(-1, trial_end)
  arrived.rew0 = min(which(data2$is_at_reward0 > 0))
  arrived.rew1 = min(which(data2$is_at_reward1 > 0))
  if (arrived.rew0 < arrived.rew1) {
    next.goal[1 : arrived.rew0] = 0
    next.goal[arrived.rew0+1 : arrived.rew1] = 1
  } else {
    next.goal[1:arrived.rew1] = 1
    next.goal[arrived.rew1+1 : arrived.rew0] = 0
  }
  
  return(next.goal) 
}

# for (trial_id_idx in levels(data$trial_id)) {
#   d = filter(data, trial_id ==trial_id_idx, cell==1)
#   is.at.reward(d$velocity, d$dist_reward0, d$timestamp)
# }
# 
# data2 = data %>%
#   group_by(trial_id, cell) %>%
#   mutate(is_at_reward1 = is.at.reward(velocity, dist_reward1, timestamp))%>%
#   mutate(is_at_reward0 = is.at.reward(velocity, dist_reward0, timestamp)) 
# data3 = data2 %>% 
#   group_by(trial_id) %>%
#   mutate(goal_reward = next.goal.reward(is_at_reward0, is_at_reward1))
