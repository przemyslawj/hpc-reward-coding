library(dplyr)
library(smoothie) # for gaussian smoothing
library(data.table)
library(reshape2) # for melt


get.entropy = function(freqs) {
  -sum(freqs * log2(freqs))
}

to.probs = function(x) { 
  x / sum(x)
}

# For performance use C++ implementation instead: getCppPlaceField
getPlaceField = function(trial.df, trace.col) {
  binSizeX = 3
  binSizeY = 3
  minX = 0
  minY = 0
  maxX = max(trial.df$smooth_trans_x, na.rm=TRUE) - minX 
  maxY = max(trial.df$smooth_trans_y, na.rm=TRUE) - minY
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

norm2 = function(x, y) {
  sqrt(x**2 + y**2)
}

# Create df with smoothed values from matrix representation
create.pf.df = function(M, occupancyM, min.zscore=0, max.zscore=2, min.occupancy.sec=1, frame.rate=20, max.y=34) {
  sigma = 2
  M1=gauss2dsmooth(M,lambda=sigma, nx=11, ny=11)
  df1 = reshape2::melt(M1) %>%
    mutate(value = ifelse(value < min.zscore, min.zscore, value)) %>%
    mutate(value = ifelse(value > max.zscore, max.zscore, value)) 
  
  #smoothedOccupancy = gauss2dsmooth(occupancyM,lambda=1, nx=3, ny=3)
  min.occupancy = min.occupancy.sec * frame.rate
  min.smoothed.occupancy = min.occupancy * 1/(2*pi*sigma^2)
  smoothedOccupancy = gauss2dsmooth(occupancyM,lambda=sigma, nx=11, ny=11)
  mid.pt = mean(1:max.y)
  df_org = reshape2::melt(smoothedOccupancy) %>%
    filter(value >= min.smoothed.occupancy) %>%
    filter(norm2(Var1 - mid.pt, Var2 - mid.pt) <= mid.pt)
  df2 = left_join(df_org, df1, by=c('Var1'='Var1', 'Var2'='Var2'), suffix=c('.occupancy', '.conv'))
    #filter(value.conv > -10.0)
  
  return(df2)
}
  
plot.pf = function(df, min.zscore=0, max.zscore=2, max.y=34) {
  jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  ggplot() +
    geom_raster(data=df, aes(x=Var1, y=max.y-Var2, fill=value.conv), interpolate=FALSE) +
    scale_fill_gradientn(colours=jet.colours(7)) +
                         #limits=c(min.zscore, max.zscore)) +
    theme_void() 
}


field.cor = function(field1, field2, make.cor.plot=FALSE, max.xy=32) {
  result = list()
  joined.fields = inner_join(field1, field2, by=c('Var1', 'Var2')) %>%
    # Avoid calculating correlation on the edges: the blue highly correlated and increases the overall correlation
    filter(Var1 < max.xy, Var2 < max.xy, Var1 > 2, Var2 > 2)  
  result$cor = cor(joined.fields$value.conv.x, joined.fields$value.conv.y)
  
  if (make.cor.plot) {
    m.conv.x = mean(joined.fields$value.conv.x)
    sd.conv.x = sd(joined.fields$value.conv.x)
    m.conv.y = mean(joined.fields$value.conv.y)
    sd.conv.y = sd(joined.fields$value.conv.y)
    
    joined.fields = mutate(
      joined.fields,
      value.conv = (value.conv.x - m.conv.x) * (value.conv.y - m.conv.y) /
        (nrow(joined.fields) - 1) / sd.conv.x / sd.conv.y)
    
    g = ggplot(joined.fields) +
      geom_raster(aes(x=Var1, y=34-Var2, fill=value.conv), interpolate=FALSE) +
      scale_fill_gradient2(low = 'blue', mid = 'white', high='red', midpoint = 0.0) +
      theme_void()
    result$g = g
  }
  
  return(result)
}