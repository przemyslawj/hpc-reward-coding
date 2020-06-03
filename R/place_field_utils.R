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
