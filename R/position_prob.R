library(dplyr)
library(readr)

source('tracking_files.R')

root_dat_dir = '/mnt/DATA/Prez/cheeseboard/2018-07-habituation'
files.df = get_tracking_files(root_dat_dir)

map_rows = 20
map_cols = 20

trial_occupancy = function(tracking.df) {
  map = matrix(rep(0, map_rows*map_cols), nrow=map_rows, ncol=map_cols)
  frames.total = 0
  for (i in 1:nrow(tracking.df)) {
    x = tracking.df$smooth_trans_x[i]
    y = tracking.df$smooth_trans_y[i]
    if ( x >= 0 & y >= 0 & x <= 100 & y <= 100) {
      x = ceiling(x / (100 / map_rows))
      y = ceiling(y / (100 / map_cols))
      map[x, y] = map[x, y] + 1
      frames.total = frames.total + 1
    }
  }
  return(list(map=map, frames=frames.total))
}

occupancy_maps = list()
total_occupancy_frames = list()
for (i in 1:nrow(files.df)) {
  tracking.df = read_csv(files.df$filepath[i], col_types = cols())
  trial_map = trial_occupancy(tracking.df)
  animal_id = as.character(files.df$animal[i])
  animal_occupancy = occupancy_maps[[animal_id]]
  animal_occupancy_frames = total_occupancy_frames[[animal_id]]
  if (is.null(animal_occupancy)) {
    animal_occupancy = matrix(rep(0, map_rows*map_cols), nrow=map_rows, ncol=map_cols)
    animal_occupancy_frames = 0
  }
  occupancy_maps[[animal_id]] = animal_occupancy + trial_map[['map']]
  total_occupancy_frames[[animal_id]] = animal_occupancy_frames + trial_map[['frames']]
}

probability_maps = list()
for (animal_id in names(occupancy_maps)) {
  probability_maps[[animal_id]] = occupancy_maps[[animal_id]] / total_occupancy_frames[[animal_id]]
}
