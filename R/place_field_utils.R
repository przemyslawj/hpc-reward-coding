library(data.table)
library(datatrace)
library(raster)
library(smoothie)

source('utils.R')


calc.spatial.info = function(binned.traces,
                             plot.dir='/tmp/pf_stability/',
                             generate.plots=FALSE,
                             nshuffles=0,
                             shuffle.shift.sec = 5,
                             # Spatial Info assumes poisson / gamma distribution of trace values
                             trace.var='deconv_trace',
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
                           shuffle.shift.sec=shuffle.shift.sec,
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


fieldCOM = function(M, field.xy, min.fieldval) {
  M.thr = M
  M.thr[M.thr < min.fieldval] = NA
  mid.pt = (dim(M)[1] + 1) / 2 
  
  # Queue for Breadth First Search on the neighbouring cells
  cell_queue = list(field.xy)
  cell_queue_index = 1
  
  com_weighted_sum = list(x=0.0, y=0.0)
  weights_sum = 0.0
  peaksize = 0
  while (cell_queue_index <= length(cell_queue)) {
    field.x = cell_queue[[cell_queue_index]]['row']
    field.y = cell_queue[[cell_queue_index]]['col']
    if (field.x >= 1 && field.x <= dim(M)[1] && 
          field.y >= 1 && field.y <= dim(M)[2] &&
          #norm2(field.x - mid.pt, field.y - mid.pt) <= mid.pt + 1.5 &&
          !is.na(M.thr[field.x, field.y])) {
      val = M.thr[field.x, field.y]
      com_weighted_sum$x = com_weighted_sum$x + val * field.x
      com_weighted_sum$y = com_weighted_sum$y + val * field.y
      weights_sum = weights_sum + val
      
      M.thr[field.x, field.y] = NA
      peaksize = peaksize + 1
      cell_queue_len = length(cell_queue)
      for (shift_x in c(-1,0,1)) {
        for (shift_y in c(-1,0,1)) {
          if (shift_x != 0 || shift_y != 0) {
            cell_queue[[length(cell_queue) + 1]] = c(field.x + shift_x, field.y + shift_y)
          }
        }
      }
    }
    cell_queue_index = cell_queue_index + 1
  }
  
  return(list(
    com=c(com_weighted_sum$x / weights_sum, com_weighted_sum$y / weights_sum),
    peaksize=peaksize))
}

fieldPeaks = function(M, 
                      minpeaksize=2, 
                      maxpeaksize=100,
                      minpeakdistance=5, 
                      min.peakheight=0.5, 
                      sigma=0,
                      min.fieldheight=0.75) {
  if (sigma > 0) {
    M1 = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  } else {
    M1 = M
  }
  minpeakval = max(1e-10, max(M1, na.rm=TRUE)) * min.peakheight
  #min.fieldval = max(M1, na.rm=TRUE) * min.fieldheight
  M.thr = M1
  M.thr[M1 < minpeakval] = NA
  r = raster(M.thr)

  # Align by 0.5 so the cells are indexed by integers
  extent(r) = extent(c(0, dim(M)[1], 0, dim(M)[2]) + 0.5)

  kernel.matrix = matrix(1, nrow=2*minpeakdistance + 1, ncol=2 * minpeakdistance + 1)
  localmax = focal(r,
                   fun=function(X) { ifelse(all(is.na(X)), NA, max(X, na.rm=TRUE)) },
                   w=kernel.matrix,
                   pad=TRUE, padValue=NA)
  r2 = r==localmax

  maxima.index = Which(r2==1, cells=TRUE)
  maxRowCol = rowColFromCell(r2, maxima.index)
  
  peaksizes = rep(0, nrow(maxRowCol))
  res = matrix(nrow=nrow(maxRowCol), ncol=2)
  colnames(res) = c('row', 'col')
  
  if (nrow(maxRowCol) == 0) {
    return(res)
  }

  for (i in 1 : nrow(maxRowCol)) {
    min.fieldval = M.thr[maxRowCol[i,1], maxRowCol[i,2]] * min.fieldheight
    com = fieldCOM(M, maxRowCol[i,], min.fieldval)
    
    res[i, ] = com$com
    peaksizes[i] = com$peaksize
  }
  
  return(res[peaksizes >= minpeaksize & peaksizes <= maxpeaksize, , drop=FALSE])
}

cell.field.peaks.info = function(field,
                                 day.rew.df,
                                 bin.size=100/nbins,
                                 max.rew.dist.thr) {
  peaks.rowcol = fieldPeaks(field, min.peakheight=0.5, sigma=0)
  peaks.rowcol = peaks.rowcol * bin.size
  
  if (nrow(peaks.rowcol) > 0) {
    min.rew.dist.res = calc.min.rew.dist(day.rew.df, peaks.rowcol[,'row'], peaks.rowcol[,'col'])
  } else {
    min.rew.dist.res = list(rew.dist=100, location.ordinal=-1)
  }
  peak.atrew.index = which(min.rew.dist.res$rew.dist <= max.rew.dist.thr)
  closer.reward.index = which.min(min.rew.dist.res$rew.dist)
  if (length(min.rew.dist.res$rew.dist[closer.reward.index]) == 0) {
    print('Error')
  }
  
  loc.ordinals.with.peaks = unique(min.rew.dist.res$location.ordinal[peak.atrew.index])
  return(list(
    maxpeakval = max(field, na.rm=TRUE),
    npeaks = nrow(peaks.rowcol),
    npeaks.at.rew = length(peak.atrew.index),
    rew.peaks.count = length(loc.ordinals.with.peaks),
    current.rew.peaks.count = sum(day.rew.df[location_ordinal %in% loc.ordinals.with.peaks, current_loc], na.rm = TRUE),
    min.rew.dist = min.rew.dist.res$rew.dist[closer.reward.index],
    closer.rew.angle = min.rew.dist.res$rew.angle[closer.reward.index],
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

geom_rewards = function(rewards.df, subject=NULL, day=NULL, nbins=20,
                        rew.colours = c('prev_rew'='gray60', 'current_rew'='white', 'TRUE'='red', 'FALSE'='blue')) {
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
