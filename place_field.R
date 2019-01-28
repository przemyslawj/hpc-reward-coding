require(dplyr)

get.entropy = function(freqs) {
  -sum(freqs * log2(freqs))
}

to.probs = function(x) { 
  x / sum(x)
}

getPlaceField = function(trial.df) {
  
  binSizeX = 17
  binSizeY = 10
  minX = 0
  minY = 0
  maxX = max(trial.df$smooth_trans_x) - minX 
  maxY = max(trial.df$smooth_trans_y) - minY
  
  trial.df = trial.df %>%
    mutate(binned_posx = round(smooth_trans_x / binSizeX)) %>%
    mutate(binned_posy = round(smooth_trans_y / binSizeY))  %>%
    arrange(timestamp)
    
  #trace = trial.df$trace - min(trial.df$trace)
  trace = trial.df$nevents
  
  totalActivityMap = matrix(0, nrow=ceiling(maxX / binSizeX), ncol=ceiling(maxY / binSizeY))
  occupancyMap = matrix(-1, nrow=ceiling(maxX / binSizeX), ncol=ceiling(maxY / binSizeY))
  
  for (i in 1:length(trace)) {
      x = max(1, trial.df$binned_posx[i])
      y = max(1, trial.df$binned_posy[i])
  
      if (occupancyMap[x, y] < 0) {
          occupancyMap[x, y] = 0
      }
      occupancyMap[x, y] = occupancyMap[x, y] + 1
      totalActivityMap[x, y] = totalActivityMap[x, y] + trial.df$nevents[i]
  }
  
  PCI = 0;
  trial.dur = nrow(trial.df)
  events.total = sum(trial.df$nevents)
  mfr = events.total / trial.dur;
  for (y in 1:dim(totalActivityMap)[2]) {
    for (x in 1:dim(totalActivityMap)[1]) {
      occupancyProb = max(0, occupancyMap[x, y]) / trial.dur
      fr = totalActivityMap[x, y] / occupancyMap[x, y]
      if (occupancyProb > 0 && fr > 0) {
        PCI = PCI + occupancyProb * fr / mfr * log2(fr / mfr)
      }
    }
  }
  Ent = to.probs(as.vector(totalActivityMap) + 0.001) %>% get.entropy
  
  field = totalActivityMap / occupancyMap
  return(list(occupancy=occupancyMap, field=field, pci=PCI, entropy=Ent))
}