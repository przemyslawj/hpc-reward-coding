library(dplyr)
library(readr)

source('crossings.R')
source('distances.R')
source('locations.R')
source('tracking_files.R')
source('tracking_info.R')

rew_zone_radius = 10
frame_rate = 24

get_time_arrived = function(at.reward.vec, timestamps) {
  indecies = which(at.reward.vec == 1)
  if (length(indecies) == 0) {
    return(-1)
  }
  return(timestamps[min(indecies)] / 1000)
}

get_stats = function(tracking_df, animal_locations.df) {
  valid_pos_df = tracking_df %>%
    filter(x > -1 & y > -1)
  
  dist_df = valid_pos_df %>%
    mutate(dist_trans = vec_dist(smooth_trans_x, smooth_trans_y),
           velocity = get_velocity(dist_trans, timestamp),
           at_rew0 = is.at.reward(velocity, dist_reward0, timestamp),
           at_rew1 = is.at.reward(velocity, dist_reward1, timestamp))
  
  mvelocity = mean(dist_df$velocity[which(dist_df$velocity > 1.0)])
  # Add crossings statistics for previous reward locations
  prev_locs.df = filter(animal_locations.df, previous_loc == TRUE, Valence == 'Positive') %>%
      arrange(position_no)
  loc_summary = list(ncross_prev_rew=NA, ncross_future_roi=NA)
  if (nrow(prev_locs.df) > 0) {
    ncross_prev_rew0 = crossings4positions(dist_df$smooth_trans_x, dist_df$smooth_trans_y,
                                           prev_locs.df$trans_x[1], prev_locs.df$trans_y[1])
    ncross_prev_rew1 = crossings4positions(dist_df$smooth_trans_x, dist_df$smooth_trans_y,
                                           prev_locs.df$trans_x[2], prev_locs.df$trans_y[2])
    loc_summary['ncross_prev_rew'] = ncross_prev_rew0 + ncross_prev_rew1
  }
  
  # Add crossings statistics for future aversive locations
  future_locs.df = filter(animal_locations.df, future_loc == TRUE, Valence == 'Negative') %>%
      arrange(position_no)
  if (nrow(future_locs.df) > 0) {
    loc_summary['ncross_future_roi'] = crossings4positions(dist_df$smooth_trans_x, dist_df$smooth_trans_y,
                                                           future_locs.df$trans_x[1], future_locs.df$trans_y[1]) 
   }
  
  total_dist = sum(dist_df$dist_trans)
  total_frames = dist_df$frame[length(dist_df$frame)] - dist_df$frame[1]
  
  rew_dwell_pct = (filter(valid_pos_df, dist_reward0 < rew_zone_radius) %>% nrow() +
                   filter(valid_pos_df, dist_reward1 < rew_zone_radius) %>% nrow()) * 
                  100.0 / total_frames 
  
  inside_roi = fill.gaps(dist_df$inside_roi, dist_df$timestamp, max.gap.ms=100)
  roi_durs = get.crossing.durations(inside_roi, dist_df$timestamp)
  roi_velocity = dist_df$velocity[which(inside_roi==1)]
  
  frame_30s = frame_rate * 30
  frame_60s = frame_rate * 60
  frame_90s = frame_rate * 90
  frame_120s = frame_rate * 120
  
  main_summary = list(
              total_dist=total_dist, 
              dist_30s=sum(dist_df$dist_trans[1:frame_30s],na.rm=TRUE),
              dist_60s=sum(dist_df$dist_trans[1:frame_60s],na.rm=TRUE),
              total_frames=total_frames, 
              crossings_n=get_crossings_n(valid_pos_df$dist_reward0, frame_rate=frame_rate) +
                          get_crossings_n(valid_pos_df$dist_reward1, frame_rate=frame_rate),
              crossings_0_30s=get_crossings_n(valid_pos_df$dist_reward0[1:frame_30s]) +
                              get_crossings_n(valid_pos_df$dist_reward1[1:frame_30s]),
              crossings_30_60s=get_crossings_n(valid_pos_df$dist_reward0[frame_30s:frame_60s]) +
                              get_crossings_n(valid_pos_df$dist_reward1[frame_30s:frame_60s]),
              crossings_60_90s=get_crossings_n(valid_pos_df$dist_reward0[frame_60s:frame_90s]) +
                              get_crossings_n(valid_pos_df$dist_reward1[frame_60s:frame_90s]),
              crossings_90_120s=get_crossings_n(valid_pos_df$dist_reward0[frame_90s:frame_120s]) +
                              get_crossings_n(valid_pos_df$dist_reward1[frame_90s:frame_120s]),
              mvelocity=mvelocity,
              nroi_crossings=length(roi_durs),
              roi_duration=sum(roi_durs),
              mroi_velocity=mean(roi_velocity),
              rew_dwell_pct=rew_dwell_pct,
              start_x=valid_pos_df$trans_x[1], 
              start_y=valid_pos_df$trans_y[1],
              arrived_rew0=get_time_arrived(dist_df$at_rew0, dist_df$timestamp),
              arrived_rew1=get_time_arrived(dist_df$at_rew1, dist_df$timestamp))
  return(c(main_summary, loc_summary))
}


output_df = data.frame(date=character(), 
                       animal=character(), 
                       trial_n=numeric(), 
                       dist=numeric(), 
                       total_frames=numeric(), 
                       is_test=factor(),
                       start_x=numeric(),
                       start_y=numeric(),
                       rew1_x=numeric(),
                       rew1_y=numeric(),
                       rew2_x=numeric(),
                       rew2_y=numeric(),
                       arrived_rew0=numeric(),
                       arrived_rew1=numeric(),
                       time_finished=numeric(),
                       crossings_n=numeric(),
                       mvelocity=numeric(),
                       ncross_prev_rew=numeric(),
                       ncross_future_roi=numeric(),
                       nroi_crossings=numeric(),
                       roi_duration=numeric(),
                       mroi_velocity=numeric(),
                       rew_dwell_pct=numeric())

root_dat_dir = '/mnt/DATA/Prez/cheeseboard/2019-02-learning'
#root_dat_dir = '~/neurodata/cheeseboard/2019-02-learning'

locations.df = read_locations(root_dat_dir) %>% 
  add_prev_locations() %>%
  add_future_neg_locations()

files.df = get_tracking_files(root_dat_dir)

for (i in 1:nrow(files.df)) {
  tracking_df = read_csv(files.df$filepath[i], col_types = cols())
  tracking_df$inside_roi = as.logical(tracking_df$inside_roi)
  print(files.df$filepath[i])
      
  animal_pos.df = filter(locations.df, Animal == files.df$animal[i], 
                         date == files.df$date[i], is_test == files.df$is_test[i])
  
  res = get_stats(tracking_df, animal_pos.df)
  current_reward_pos = filter(animal_pos.df, Valence == 'Positive', 
                              current_loc == TRUE, is_test == files.df$is_test[i])
  time_finished_sec = -1
  if (res[['arrived_rew0']] >= 0 && res[['arrived_rew1']] >= 0) {
    time_finished_sec = max(res[['arrived_rew0']], res[['arrived_rew1']])
  }
  
  new_row = data.frame(date=files.df$date[i],
                       animal_id=files.df$animal[i],
                       trial_n=files.df$trial[i],
                       dist=res['total_dist'],
                       dist_30s=res['dist_30s'],
                       dist_60s=res['dist_60s'],
                       total_frames=res['total_frames'],
                       is_test=files.df$is_test[i],
                       location_set=current_reward_pos$location_set[1],
                       start_x=res['start_x'],
                       start_y=res['start_y'],
                       rew1_x=current_reward_pos$trans_x[1],
                       rew1_y=current_reward_pos$trans_y[1],
                       rew2_x=current_reward_pos$trans_x[2],
                       rew2_y=current_reward_pos$trans_y[2],
                       arrived_rew0=res[['arrived_rew0']],
                       arrived_rew1=res[['arrived_rew1']],
                       time_finished=time_finished_sec,
                       crossings_n=res['crossings_n'],
                       mvelocity=res['mvelocity'],
                       crossings_0_30s=res['crossings_0_30s'],
                       crossings_30_60s=res['crossings_30_60s'],
                       crossings_60_90s=res['crossings_60_90s'],
                       crossings_90_120s=res['crossings_90_120s'],
                       ncross_prev_rew = res$ncross_prev_rew[1],
                       ncross_future_roi = res$ncross_future_roi[1],
                       nroi_crossings=res['nroi_crossings'],
                       roi_duration=res['roi_duration'],
                       mroi_velocity=res['mroi_velocity'],
                       rew_dwell_pct=res['rew_dwell_pct'])
  output_df = rbind(output_df, new_row) 
}

output_df$animal_id = as.factor(output_df$animal_id)
output_df$date = as.Date(output_df$date)
output_df$trial_n = as.integer(output_df$trial_n)

# Add absolute trial id/number
fst_loc_date = as.Date('2019-02-20')
snd_loc_date = as.Date('2018-03-04')
output_df = output_df %>% 
    mutate(learning_day1 = as.integer(date - fst_loc_date) + 1,
           learning_day2 = as.integer(date - snd_loc_date) + 1,
           learning_day = if_else(location_set == 1, learning_day1, learning_day2))  %>%
    mutate(trial_id = learning_day1 * 10 + trial_n) %>%
    select(-learning_day2) 
output_df$learning_day = as.integer(output_df$learning_day)
output_df$learning_day1 = as.integer(output_df$learning_day1)
output_df$trial_id = as.integer(output_df$trial_id)

write_csv(output_df, paste0(root_dat_dir, '/trial_stats.csv'))

