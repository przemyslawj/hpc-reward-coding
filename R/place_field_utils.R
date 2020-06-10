library(data.table)
library(datatrace)
library(raster)
library(smoothie)


calc.spatial.info = function(binned.traces,
                             plot.dir='/tmp/pf_stability/',
                             generate.plots=FALSE,
                             nshuffles=0,
                             trace.var='trace',
                             timebin.dur.msec=100,
                             min.occupancy.sec=1,
                             nbins=30,
                             gaussian.var=2) {
  binned.traces = data.table(binned.traces)
  cells = unique(binned.traces$cell_id) %>% sort
  pci.df = data.frame()
  fields = list()
  occupancies = list()

  for (cell_name in cells) {
    cell.df = binned.traces[cell_id == cell_name ,]
    pf = cell.spatial.info(cell.df, nbins, nbins, generate.plots, nshuffles,
                           trace.var=trace.var,
                           bin.hz=1000/timebin.dur.msec,
                           shuffle.shift.sec=5,
                           min.occupancy.sec=min.occupancy.sec,
                           kernel.size = 9,
                           gaussian.var = gaussian.var)
    if (length(pf$cell_info) > 0) {
      fields[[format(cell_name)]] = pf$field
      occupancies[[format(cell_name)]] = pf$occupancy
      pci.df = bind_rows(pci.df, pf$cell_info)

      if (!is.na(plot.dir) && generate.plots && !is.na(pf$g)) {
        if (!file.exists(plot.dir)) {
          dir.create(plot.dir, recursive=TRUE)
        }
        ggsave(paste0(plot.dir, 'place_field_cell_', cell_name, '.jpg'), pf$g,
               width=3.5, height=3.0, units='cm', dpi=300)
      }
    }
  }

  return(list(df=pci.df, field=fields, occupancy=occupancies))
}


fieldPeaks = function(M, minpeaksize=3, minpeakdistance=4, min.peakheight=0.33, sigma=1.4) {
  if (sigma > 0) {
    M1 = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  } else {
    M1 = M
  }
  minpeakval = max(M1, na.rm=TRUE) * min.peakheight
  M.thr = M1
  M.thr[M1 < minpeakval] = NA
  r = raster(M.thr)

  # Align by 0.5 so the cells are index by integers
  extent(r) = extent(c(0, dim(M)[1], 0, dim(M)[2]) + 0.5)

  kernel.matrix = matrix(1, nrow=2*minpeakdistance + 1, ncol=2 * minpeakdistance + 1)
  localmax = focal(r,
                   fun=function(X) { ifelse(all(is.na(X)), NA, max(X, na.rm=TRUE)) },
                   w=kernel.matrix,
                   pad=TRUE, padValue=NA)

  peaksizes = focal(r,
                    fun=function(X) { sum(X > minpeakval * 0.25, na.rm=TRUE) },
                    w=kernel.matrix,
                    pad=TRUE, padValue=NA)
  r2 = r==localmax

  maxima.index = Which(r2==1, cells=TRUE)
  keep.peaks = which(peaksizes[maxima.index] >= minpeaksize)
  maxRowCol = rowColFromCell(r2, maxima.index[keep.peaks])

  return(maxRowCol)
}

findPeaksAtRewards = function(M, rewLocXY, minpeaksize=3, min.peakheight=0.20, max.rewdistance=2) {
  sigma = 1.4
  M1 = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  minpeakval = max(M1, na.rm=TRUE) * min.peakheight
  M.thr = M1
  M.thr[M1 < minpeakval] = NA
  r = raster(M.thr)

  rew.M = matrix(0, nrow=dim(M)[1], ncol=dim(M)[2])
  rew.mask = raster(rew.M)

  # Align by 0.5 so the cells are index by integers
  extent(rew.mask) = extent(c(0, dim(M)[1], 0, dim(M)[2]) + 0.5)
  extent(r) = extent(c(0, dim(M)[1], 0, dim(M)[2]) + 0.5)

  rew.cells = cellFromRowCol(r, rewLocXY$x, rewLocXY$y)
  rew.mask[rew.cells] = 1

  kernel.matrix = matrix(1, nrow=2*minpeakdistance + 1, ncol=2 * minpeakdistance + 1)
  rew.kernel = matrix(1, nrow=2*max.rewdistance + 1, ncol=2*max.rewdistance + 1)

  #create rings around the rewards
  rew.mask.calc = focal(rew.mask,
                        fun=max,
                        w=rew.kernel,
                        pad=TRUE,
                        padValue=0)
  r.masked = mask(r, rew.mask.calc, maskvalue=0)

  localmax = focal(r.masked,
                   fun=function(X) { ifelse(all(is.na(X)), NA, max(X, na.rm=TRUE)) },
                   w=kernel.matrix,
                   pad=TRUE, padValue=NA)


  # peak size f
  peaksizes = focal(r,
                    fun=function(X) { sum(X >= minpeakval * 0.25, na.rm=TRUE) },
                    w=kernel.matrix,
                    pad=TRUE, padValue=NA) %>%
              focal(fun=max, w=kernel.matrix, pad=TRUE, padValue=0)
  r2 = r.masked==localmax

  maxima.index = Which(r2==1, cells=TRUE)
  keep.peaks = which(peaksizes[maxima.index] >= minpeaksize)
  maxRowCol = rowColFromCell(r2, maxima.index[keep.peaks])

  return(maxRowCol)
}


cell.field.peaks.info = function(field,
                                 day.rew.df,
                                 bin.size=100/nbins,
                                 max.rew.dist.thr) {
  peaks.rowcol = fieldPeaks(field, min.peakheight=0.5, sigma=1.4)
  peaks.rowcol = peaks.rowcol * bin.size
  
  if (nrow(peaks.rowcol) > 0) {
    min.rew.dist.res = calc.min.rew.dist(day.rew.df, peaks.rowcol[,'row'], peaks.rowcol[,'col'])
  } else {
    min.rew.dist.res = list(rew.dist=100, location.ordinal=-1)
  }
  peak.atrew.index = which(min.rew.dist.res$rew.dist <= max.rew.dist.thr)
  closer.reward.index = which.min(min.rew.dist.res$rew.dist)
  
  loc.ordinals.with.peaks = unique(min.rew.dist.res$location.ordinal[peak.atrew.index])
  return(list(
    maxpeakval = max(field, na.rm=TRUE),
    npeaks = nrow(peaks.rowcol),
    npeaks.at.rew = length(peak.atrew.index),
    rew.peaks.count = length(loc.ordinals.with.peaks),
    current.rew.peaks.count = sum(day.rew.df[location_ordinal %in% loc.ordinals.with.peaks, current_loc], na.rm = TRUE),
    min.rew.dist = min.rew.dist.res$rew.dist[closer.reward.index],
    closer.rew.ordinal = min.rew.dist.res$location.ordinal[closer.reward.index]
  ))
}

calc.field.peaks.info = function(day,
                                 animal_name,
                                 field.list,
                                 rew.df,
                                 max.rew.dist.thr) {
  if (nrow(rew.df) == 0) {
    #print(paste0('No rewards to compare for animal=', animal_name, ' date=', day))
    return(data.frame())
  }
  bin.size = 100 / nbins
  cell_names = names(field.list)
  if (length(cell_names) == 0) {
    #print(paste0('Empty list of fields for animal=', animal_name, ' date=', day))
    return(data.frame())
  }
  
  map_dfr(cell_names, ~ {
    field.info = cell.field.peaks.info(field.list[[.x]], rew.df, bin.size, max.rew.dist.thr)
    meta = list(animal=animal_name, date=day, cell_id=as.integer(.x))
    return(append(meta, field.info))
  })
}


day.normalize.fields = function(pf.df) {
  pf.df %>%
    group_by(animal, date, day_desc, cell_id) %>%
    dplyr::mutate(value.field = value.field - min(value.field, na.rm=TRUE),
                  value.field = value.field / max(value.field, na.rm=TRUE))
}

geom_rewards = function(rewards.df, subject=NULL, day=NULL, nbins=20) {
  dayreward.df = rewards.df
  if (!is.null(subject)) {
    dayreward.df = filter(rewards.df, animal==subject)
  }
  if (!is.null(day)) {
    dayreward.df = filter(dayreward.df, date==day)
  }
  if (!('current_loc' %in% names(dayreward.df))) {
    dayreward.df$current_loc = TRUE
  }
  dayreward.df = dplyr::mutate(
    dayreward.df,
    rew_name = ifelse(current_loc, 'current_rew', 'prev_rew'))
  rew.colours = c('prev_rew'='gray60', 'current_rew'='white', 'TRUE'='red', 'FALSE'='blue')
  #rew.shape = 0
  rew.shape = 2
  rew.size = 2
  bin.size = 100 / nbins
  list(geom_point(data=dayreward.df, 
                  aes(x=trans_x / bin.size, y=(100 - trans_y) / bin.size, color=rew_name), 
                  shape=rew.shape,
                  size=rew.size,
                  stroke=2),
       scale_color_manual(values=rew.colours))
}
