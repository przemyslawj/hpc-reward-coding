source('plotting_params.R')

library(dplyr)
library(datatrace)
library(R.matlab)
library(imager)
ca_img_result_dir = '/mnt/DATA/Prez/cheeseboard-down/down_2/2020-10/habituation/2020_10_08/caiman/M-BR/'
#ca_img_result_dir = '/mnt/DATA/Prez/cheeseboard-down/down_2/2019-07/habituation/2019-07-18/caiman/D-BR/'
#ca_img_result_dir = '/mnt/DATA/Prez/cheeseboard-down/down_2/2020-01/habituation/2020-01-28/caiman/G-BR/'

mat_file = paste0(ca_img_result_dir, 'ms.mat')
ca.mat = readMat(mat_file)
meanFrame = matrix(ca.mat$ms[,,1]$meanFrame, 
                   nrow = ca.mat$ms[,,1]$height, ncol=ca.mat$ms[,,1]$width) 

filtered_mat_file = file.path(ca_img_result_dir, 'filtered', 'ms.mat')
ca.mat = readMat(filtered_mat_file)
dat.mat = ca.mat$ms[,,1]

filtered.SFPs = dat.mat$SFPs
for (i in 1:dim(dat.mat$SFPs)[3]) {
  M = filtered.SFPs[,,i]
  thr.val = max(M) * 0.5
  M.thr = M
  M.thr[M<thr.val] = 0
  filtered.SFPs[,,i] = M.thr
}

grays = rgb(red = 0:255/255, blue = 0:255/255, green = 0:255/255)

M.sum = apply(filtered.SFPs, c(1,2), max) #/ dim(filtered.SFPs)[3]
image(log(1-M.sum), col=grays)
as.cimg(1-M.sum) %>% plot()

selected.cell_ids = seq(1,100,10)
#selected.cell_ids = seq(1,71,7) 
M.minus.selected = apply(filtered.SFPs[,,-selected.cell_ids], c(1,2), sum)
M.selected = apply(filtered.SFPs[,,selected.cell_ids], c(1,2), sum)

#im = rgb(M,M,M)
#dim(im) = dim(M)
M.cimg = as.cimg(c(M.minus.selected, M.minus.selected, M.sum + 3*M.selected), x=nrow(M.sum), y=ncol(M.sum), cc=3) %>%
  imager::resize_tripleXY() %>%
  imager::imsharpen(3) %>%
  imager::isoblur(1.5)
M.meanframe = as.cimg(meanFrame / max(meanFrame)) %>%
  imager::resize_tripleXY() %>%
  identity()
#imager::imsharpen(3) %>%
#imager::isoblur(1.5)
plot(M.meanframe, rescale=FALSE)

plot(M.cimg)
imdraw(M.meanframe, M.cimg, opacity=0.85) %>% 
  plot()


data.traces = read.data.trace(file.path(ca_img_result_dir, 'filtered'))

g1 = data.traces[cell %in% selected.cell_ids & exp_title=='trial' & trial==2,] %>%
  ggplot(aes(x=timestamp/1000, y=trace, color=cell_id)) +
  geom_line(size=0.5) +
  gtheme +
  facet_grid(cell_id ~ .) +
  scale_color_continuous(low='#81b1e4', high=low2high.colours[2]) +
  xlim(c(200,400))

g2 = data.traces[cell %in% selected.cell_ids & exp_title=='trial' & trial==2,] %>%
  ggplot(aes(x=timestamp/1000, y=deconv_trace, color=cell_id)) +
  geom_line(size=0.5) +
  gtheme +
  facet_grid(cell_id ~ .) +
  scale_color_continuous(low='#81b1e4', high=low2high.colours[2]) +
  xlim(c(200,400))

cowplot::plot_grid(g1, g2, nrow=1)
ggsave('/home/prez/tmp/cheeseboard/examples_traces.svg', units='cm',
       width=18, height=8)
