require(dplyr)

getPlaceField = function(trial.df) {
  
  minX = 0
  minY = 0
  maxX = max(trial.df$smooth_trans_x) - minX
  maxY = max(trial.df$smooth_trans_y) - minY
  binSize = 10
  
  trial.df = trial.df %>%
    mutate(binned_posx = round(smooth_trans_x / binSize)) %>%
    mutate(binned_posy = round(smooth_trans_y / binSize)) 
    
  #trace = trial.df$trace - min(trial.df$trace)
  trace = trial.df$nevents
  
  totalActivityMap = matrix(0, nrow=ceiling(maxX / binSize), ncol=ceiling(maxY / binSize))
  occupancyMap = matrix(-1, nrow=ceiling(maxX / binSize), ncol=ceiling(maxY / binSize))
  
  prevTimestamp = 0
  for (i in 1:length(trace)) {
      x = max(1, trial.df$binned_posx[i])
      y = max(1, trial.df$binned_posy[i])
  
      if (occupancyMap[x, y] < 0) {
          occupancyMap[x, y] = 0
      }
      occupancyMap[x, y] = occupancyMap[x, y] + trial.df$timestamp[i] - prevTimestamp
      prevTimestamp = trial.df$timestamp[i]
      totalActivityMap[x, y] = totalActivityMap[x, y] + trace[i]
  }
  
  field = totalActivityMap / occupancyMap
  return(list(occupancy=occupancyMap, field=field))
}