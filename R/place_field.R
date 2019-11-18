require(dplyr)
#library(mmand) 
library(smoothie) # for gaussian smoothing

library(data.table)

library(reshape2) # for melt

get.entropy = function(freqs) {
  -sum(freqs * log2(freqs))
}

to.probs = function(x) { 
  x / sum(x)
}

getPlaceField = function(trial.df, trace.col) {
  binSizeX = 3
  binSizeY = 3
  minX = 0
  minY = 0
  maxX = max(trial.df$smooth_trans_x, na.rm=TRUE) - minX 
  maxY = max(trial.df$smooth_trans_y, na.rm=TRUE) - minY
  
  #trial.df = trial.df %>%
  #  mutate(binned_posx = round(smooth_trans_x / binSizeX)) %>%
  #  mutate(binned_posy = round(smooth_trans_y / binSizeY))  %>%
  #  arrange(timestamp)
  trial.df = data.table(trial.df)
    
  trace = trial.df[, ..trace.col][[1]]
  
  totalActivityMap = matrix(0, nrow=ceiling(maxX / binSizeX), ncol=ceiling(maxY / binSizeY))
  occupancyMap = matrix(0, nrow=ceiling(maxX / binSizeX), ncol=ceiling(maxY / binSizeY))
  
  for (i in 1:length(trace)) {
    x = max(1, round(trial.df$smooth_trans_x[i] / binSizeX))
    y = max(1, round(trial.df$smooth_trans_y[i] / binSizeY))
    #x = max(1, trial.df$binned_posx[i])
    #y = max(1, trial.df$binned_posy[i])
    occupancyMap[x, y] = occupancyMap[x, y] + 1
    totalActivityMap[x, y] = totalActivityMap[x, y] + trace[i]
  }
  
  SI = 0
  frame_rate = 20
  fr = totalActivityMap / occupancyMap 
  fr_offset = min(fr, na.rm=TRUE) - 10e-10
  fr = fr - fr_offset
  
  mfr = mean(trace) - fr_offset;
  field = matrix(NA, nrow=nrow(totalActivityMap), ncol=ncol(totalActivityMap))
  for (y in 1:dim(totalActivityMap)[2]) {
    for (x in 1:dim(totalActivityMap)[1]) {
      occupancyProb = max(0, occupancyMap[x, y], na.rm=TRUE) / nrow(trial.df)
      if (occupancyProb > 0) {
        SI = SI + occupancyProb * fr[x,y] / mfr * log2(fr[x,y] / mfr)
        field[x,y] = totalActivityMap[x,y] / occupancyMap[x,y]
      }
    }
  }
  Ent = to.probs(as.vector(totalActivityMap) + 0.001) %>% get.entropy
  
  field[is.na(field)] = mean(field, na.rm=TRUE)
  return(list(occupancy=occupancyMap, 
              activity=totalActivityMap, 
              field=field, 
              spatial.information=SI, 
              spatial.information.perspike=SI/mfr,
              entropy=Ent))
}

# Plot smoothed values from matrix representation
plot.pf = function(M, occupancyM, min.zscore=-3, max.zscore = 3) {
  #M1 = gaussianSmooth(M, 2)
  M1=gauss2dsmooth(M,lambda=2, nx=15, ny=15)
  df1 = reshape2::melt(M1) %>%
    mutate(value = ifelse(value < min.zscore, min.zscore, value)) %>%
    mutate(value = ifelse(value > max.zscore, max.zscore, value)) 
  
  smoothedOccupancy = gauss2dsmooth(occupancyM,lambda=1, nx=3, ny=3)
  df_org = reshape2::melt(smoothedOccupancy) %>%
    filter(value > 1)
  #df_org = reshape2::melt(M)  %>%
    #filter(value > -0.02)
  df2 = left_join(df_org, df1, by=c('Var1'='Var1', 'Var2'='Var2'), suffix=c('.occupancy', '.conv')) %>%
    #filter(value.conv > -0.015)
    filter(value.conv > -10.015)
    
  
  jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  df2 %>%
    ggplot(aes(x=Var1, y=Var2)) +
    geom_raster(aes(fill=value.conv), interpolate=FALSE) +
    scale_fill_gradientn(colours=jet.colours(7),
                         limits=c(min.zscore, max.zscore)) +
    theme_void() 
}