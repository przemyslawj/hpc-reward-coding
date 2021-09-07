library(purrr)
library(stringr)

source('plotting_params.R')

v4.px2um = 4.31/2
v3.px2um = 2.84/2
v4.animals = c('M-BL', 'N-BL', 'N-BR')

reg_stat_files_dir = '/mnt/DATA/Prez/cheeseboard-down/down_2/registration_stats'
reg_stat_files = list.files(reg_stat_files_dir, pattern='*.csv')
reg_stat_df = map_dfr(reg_stat_files, ~ {
  df = read.csv(file.path(reg_stat_files_dir, .x))
  df$animal = stringr::str_replace(.x, '.csv', '')
  px2um = v3.px2um
  if (df$animal[1] %in% v4.animals) {
    px2um = v4.px2um
  }
  df$reg_dist_um = df$reg_dist_px * px2um
  df
})

reg_stat_df %>% 
  left_join(mouse.meta.df) %>%
  ggplot(aes(x=reg_dist_px, y=matched_cors)) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) + scale_fill_continuous(low=low2high.colours[1], high=low2high.colours[2]) +
  #geom_point(size=0.2, alpha=0.2, shape=1) +
  #facet_wrap(animal ~ ., scales = 'free') +
  #facet_grid(implant ~ .) +
  ylim(c(0, 1.0)) + xlim(c(0, 6)) +
  gtheme +
  ylab('Spatial correlation') + xlab('Centroid distance (μm)')
ggsave('cell_registration_matched_stat.pdf', device=cairo_pdf,
       path='/home/prez/tmp/cheeseboard', units='cm', height=3.3, width=5.0)


reg_stat_df$reg_dist_um %>% quantile(c(0.25, 0.5, 0.75))
reg_stat_df$matched_cors %>% quantile(c(0.25, 0.5, 0.75))
