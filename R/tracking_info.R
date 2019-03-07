get.arrived.indecies = function(at.reward) {
  reward.arrived = which(diff(at.reward)==1)
  reward.left = which(diff(at.reward)==-1)
  
  if (at.reward[1]) {
    reward.arrived = c(1, reward.arrived)
  }
  if (length(reward.left) < length(reward.arrived)) {
    reward.left = c(reward.left, length(at.reward))
  }
  if (length(reward.left) != length(reward.arrived)) {
    print('Warning: incorrect lengths')
    print(trial_id[1])
  }
  
  return(list(arrived=reward.arrived, left=reward.left))
}

is.at.reward = function(velocities, rew.dist, timestamps, trial_id, cell, running.thresh = 1, rew.dist.thresh=10) {
  min.duration.at.reward.ms = 3000
  at.reward = rep(0, length(velocities))
  res = rep(0, length(velocities))
  at.reward = velocities < running.thresh &
      rew.dist < rew.dist.thresh
  
  reward.indecies = get.arrived.indecies(at.reward)
  reward.arrived = reward.indecies[['arrived']]
  reward.left = reward.indecies[['left']]
  
  if (length(reward.arrived) > 1) {
    for (i in 1:(length(reward.left) - 1)) {
      left_dur = timestamps[reward.arrived[i+1]] -  timestamps[reward.left[i]]
      if (left_dur < 1000) {
        at.reward[reward.left[i] : reward.arrived[i+1]] = 1
      }
    }
  }
  
  reward.indecies = get.arrived.indecies(at.reward)
  reward.arrived = reward.indecies[['arrived']]
  reward.left = reward.indecies[['left']]
  
  if (length(reward.left) > 0) {
    at.reward.duration.ms = timestamps[reward.left] - timestamps[reward.arrived]
    keep_indecies = which(at.reward.duration.ms > min.duration.at.reward.ms)
    for (i in keep_indecies) {
      res[reward.arrived[i]:reward.left[i]] = 1
    }
  }
  
  return(res)
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

for (trial_id_idx in levels(data$trial_id)) {
  d = filter(data, trial_id ==trial_id_idx, cell==1)
  is.at.reward(d$velocity, d$dist_reward0, d$timestamp, d$trial_id, d$cell)
}

data2 = data %>%
  group_by(trial_id, cell) %>%
  mutate(is_at_reward1 = is.at.reward(velocity, dist_reward1, timestamp, trial_id, cell))%>%
  mutate(is_at_reward0 = is.at.reward(velocity, dist_reward0, timestamp, trial_id, cell)) 
data3 = data2 %>% 
  group_by(trial_id) %>%
  mutate(goal_reward = next.goal.reward(is_at_reward0, is_at_reward1))
