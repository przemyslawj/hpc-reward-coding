library(dplyr)
library(smoothie) # for gaussian smoothing
library(data.table)
library(reshape2) # for melt
Rcpp::sourceCpp('place_field.cpp')

source('utils.R')

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
create.pf.df = function(M, occupancyM, max.xy, min.occupancy.sec=1, frame.rate=20) {
  sigma = 1.4
  M1=gauss2dsmooth(M,lambda=sigma, nx=11, ny=11)
  df1 = reshape2::melt(M1) 
  
  min.occupancy = min.occupancy.sec * frame.rate
  min.smoothed.occupancy = min.occupancy * 1/(2*pi*sigma^2)
  smoothedOccupancy = gauss2dsmooth(occupancyM,lambda=sigma, nx=11, ny=11)
  mid.pt = mean(1:max.xy)
  df_org = reshape2::melt(smoothedOccupancy) %>%
    filter(value >= min.smoothed.occupancy) %>%
    filter(norm2(Var1 - mid.pt, Var2 - mid.pt) <= mid.pt)
  df2 = left_join(df_org, df1, by=c('Var1'='Var1', 'Var2'='Var2'), suffix=c('.occupancy', '.conv'))
  
  return(df2)
}
  
plot.pf = function(df, max.xy) {
  jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  ggplot() +
    geom_raster(data=df, aes(x=Var1, y=max.xy-Var2, fill=value.conv), interpolate=FALSE) +
    scale_fill_gradientn(colours=jet.colours(7)) +
    xlim(c(0, max.xy)) + ylim(c(0, max.xy)) +
    theme_void() 
}


field.cor = function(field1, field2, max.xy, make.cor.plot=FALSE) {
  result = list()
  joined.fields = inner_join(field1, field2, by=c('Var1', 'Var2')) %>%
    # Avoid calculating correlation on the edges: the blue highly correlated and increases the overall correlation
    filter(Var1 < max.xy - 1, Var2 < max.xy - 1, Var1 > 2, Var2 > 2)  
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
      geom_raster(aes(x=Var1, y=max.xy-Var2, fill=value.conv), interpolate=FALSE) +
      scale_fill_gradient2(low = 'blue', mid = 'white', high='red', midpoint = 0.0) +
      xlim(c(0, max.xy)) + ylim(c(0, max.xy)) +
      theme_void()
    result$g = g
  }
  
  return(result)
}


cell.spatial.info = function(cell.df, generate.plots=FALSE, nshuffles=0,
                             trace.col='trace',
                             timebin.size=1,
                             frame.hz=20) {
  nbins.xy = getNBinsXY()
  cell.events = cell.df[nevents > 0,]
  trace.vals = cell.df[[trace.col]]
  trace.vals = trace.vals - min(trace.vals)

  if (length(trace.vals) == 0) {
    return(list(cell_info=data.frame(),
                field=matrix(),
                occupancy=matrix(),
                g=NULL))
  }
  cell_name = cell.df$cell_id[1]
  trial_ends = get.trial.ends(cell.df$timestamp)
  #trace.quantiles = quantile(trace.vals, c(0.85, 0.95, 1.0), na.rm=TRUE) %>% unname + 0.01
  #trace.quantiles = c(0.5, 1000.0)
  trace.quantiles = quantile(trace.vals, c(0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0), na.rm=TRUE)  + 0.01
  pf = with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, trace.vals,
                                      c(trace.quantiles[["20%"]], trace.quantiles[["50%"]],
                                        trace.quantiles[["80%"]], trace.quantiles[["90%"]], 
                                        trace.quantiles[["95%"]], trace.quantiles[["100%"]]), 
                                      trial_ends, nshuffles,
                                      2 * frame.hz,
                                      timebin.size))

  si.signif.thresh = quantile(pf$shuffle.si, 0.95, na.rm=TRUE)[[1]]
  si.signif = pf$spatial.information >= si.signif.thresh
  mi.signif.thresh = quantile(pf$shuffle.mi, 0.95, na.rm=TRUE)[[1]]
  mi.signif = pf$spatial.information >= mi.signif.thresh

  cell_info = list(cell_id=cell_name,
                   spatial.information=pf$spatial.information,
                   si.signif.thresh=si.signif.thresh,
                   signif.si=si.signif,
                   mi.signif.thresh=mi.signif.thresh,
                   signif.mi=mi.signif,
                   #field.centre.x=pf$field.centre[1],
                   #field.centre.y=pf$field.centre[2],
                   #field.max.x=pf$field.max.xy[1],
                   #field.max.y=pf$field.max.xy[2],
                   #field.max=pf$field.max,
                   field.size.50=pf$field.size.50,
                   field.size.25=pf$field.size.25,
                   spatial.information.perspike=pf$spatial.information.perspike,
                   mfr=pf$mfr,
                   mutual.info=pf$mutual.info,
                   mutual.info.bias=pf$mutual.info.bias,
                   space.sampling.factor = pf$space.sampling.factor,
                   sparsity = pf$sparsity,
                   quantile20=trace.quantiles[["20%"]],
                   quantile99=trace.quantiles[["99%"]],
                   nevents=nrow(cell.events))

  # find field max value and pos in the smoothed values
  pf.df = create.pf.df(pf$field, pf$occupancy, max.xy=nbins.xy)
  max.row = pf.df[which.max(pf.df$value.conv),]
  cell_info$field.max = max.row$value.conv
  cell_info$field.mean = mean(pf.df$value.conv)
  cell_info$field.max.x = max.row$Var1 / nbins.xy * 100
  cell_info$field.max.y = max.row$Var2/ nbins.xy * 100
  
  g.placefield=NA
  if (generate.plots) {
    cell_event_rate = nrow(cell.events) / nrow(cell.df) * frame.hz

    g.placefield = plot.pf(pf.df, max.xy=getNBinsXY()) +
      labs(title=paste0('Cell ', cell_name, ' MER = ', format(cell_event_rate, digits=2), ' Hz',
                        '\nSI = ', format(pf$spatial.information, digits=2),
                        ifelse(si.signif, '*', ''),
                        ' SI per spike = ', format(pf$spatial.information.perspike, digits=2),
                        '\nMI - bias = ', format(pf$mutual.info - pf$mutual.info.bias, digits=3),
                        ifelse(mi.signif, '*', ''),
                        '\nsparsity = ', format(pf$sparsity, digits=2)),
           fill='Deconv trace') +
      theme(text = element_text(size=4),
            plot.title = element_text(size=4))


  }

  return(list(cell_info=cell_info,
             field=pf$field,
             occupancy=pf$occupancy,
             g=g.placefield))
}
