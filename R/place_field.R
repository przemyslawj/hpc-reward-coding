library(dplyr)
library(smoothie) # for gaussian smoothing
library(data.table)
library(reshape2) # for melt
Rcpp::sourceCpp('place_field.cpp')

source('utils.R')


norm2 = function(x, y) {
  sqrt(x**2 + y**2)
}


to.matrix = function(M, max.xy) {
  if (!is.matrix(M)) {
    M = matrix(M, nrow=max.xy, ncol=max.xy, byrow=TRUE) 
  }
  return(M)
}

# Create df with smoothed values from matrix representation
create.pf.df = function(M, occupancyM, max.xy, min.occupancy.sec=1, frame.rate=20, sigma = 1.4) {

  M = to.matrix(M, max.xy)
  occupancyM = to.matrix(occupancyM, max.xy)
  M1 = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  df1 = reshape2::melt(M1) 
  
  min.occupancy = min.occupancy.sec * frame.rate
  min.smoothed.occupancy = min.occupancy * 1/(2*pi*sigma^2)
  smoothedOccupancy = gauss2dsmooth(occupancyM, lambda=sigma, nx=11, ny=11)
  mid.pt = mean(1:max.xy)
  df_org = reshape2::melt(smoothedOccupancy) %>%
    dplyr::filter(value >= min.smoothed.occupancy) %>%
    dplyr::filter(norm2(Var1 - mid.pt, Var2 - mid.pt) <= mid.pt)
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
                             bin.hz=5,
                             trace.var='trace',
                             min.occupancy.sec=1) {
  nbins.xy = getNBinsXY()
  nevents = nrow(cell.df[nevents > 0,])
  trace.vals = cell.df[[trace.var]]
  trace.vals = trace.vals - min(trace.vals)

  if (length(trace.vals) == 0) {
    return(list(cell_info=list(),
                field=matrix(),
                occupancy=matrix(),
                g=NULL))
  }
  cell_name = cell.df$cell_id[1]
  trial_ends = get.trial.ends(cell.df$time_bin)
  #trace.quantiles = quantile(trace.vals, c(0.85, 0.95, 1.0), na.rm=TRUE) %>% unname + 0.01
  #trace.quantiles = c(0.5, 1000.0)
  trace.quantiles = quantile(trace.vals, c(0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0), na.rm=TRUE)  + 0.01
  max.traceval = max(trace.vals)
  trace.bins = c(trace.quantiles[["20%"]], trace.quantiles[["50%"]],
                 trace.quantiles[["80%"]], trace.quantiles[["90%"]], 
                 trace.quantiles[["95%"]], trace.quantiles[["100%"]])
  #trace.bins = c(max(0.1, trace.quantiles[["50%"]]),
  #               trace.quantiles[["80%"]],
  #               trace.quantiles[["95%"]],
  #               trace.quantiles[["99%"]],
  #               trace.quantiles[["100%"]])
  pf = getCppPlaceField(to_1dim(cell.df$bin.x, cell.df$bin.y), 
                        trace.vals,
                        trace.bins,
                        trial_ends, 
                        nshuffles,
                        2 * bin.hz,
                        min.occupancy.sec * bin.hz) # min occupancy

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
                   nevents=nevents)

  # find field max value and pos in the smoothed values
  pf.df = create.pf.df(pf$field, pf$occupancy, max.xy=nbins.xy, frame.rate=bin.hz)
  if (nrow(pf.df) > 0) {
    max.row = pf.df[which.max(pf.df$value.conv),]
  } else { # no bin with high enough occupancy
    max.row=list(value.conv=0, Var1=-1, Var2=-1)
  }
  cell_info$field.max = max.row$value.conv
  cell_info$field.mean = mean(pf.df$value.conv)
  cell_info$field.max.x = max.row$Var1 / nbins.xy * 100
  cell_info$field.max.y = max.row$Var2/ nbins.xy * 100
  
  g.placefield=NA
  if (generate.plots) {
    cell_event_rate = nevents / nrow(cell.df) * bin.hz

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
             field=to.matrix(pf$field, nbins.xy),
             occupancy=to.matrix(pf$occupancy, nbins.xy),
             g=g.placefield))
}
