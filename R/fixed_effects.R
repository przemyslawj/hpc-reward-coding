library(dplyr)
library(ggplot2)
library(lmerTest)


plot.model.diagnostics = function(m, animals, exp_names) {
  df = data.frame(res=unname(residuals(m)), 
                  fitted=unname(fitted(m)), 
                  animal=animals, 
                  exp_name=exp_names) 
  
  ggplot(df) +
    geom_point(aes(x=res, y=fitted, colour=animal, shape=exp_name)) +
    theme(legend.position = 'none') -> g1
  
  X = qqnorm(residuals(m), plot.it=FALSE)
  data.frame(x=X$x, y=X$y, animal=animals, exp_name=exp_names) %>%
    ggplot() +
    geom_point(aes(x=x, y=y, colour=animal, shape=exp_name)) +
    stat_qq_line(mapping=aes(sample=res), data=df) +
    xlab('Theoretical quantiles') +
    ylab('Sample quantiles') +
    theme(legend.position = 'right') -> g2
  
  plot_grid(g1,g2)
}

create.animal.summary = function(df, var) {
  var = enquo(var)
  group_by(df, implant, animal, exp) %>% 
    dplyr::summarise(var.mean=mean(!!var, na.rm=TRUE),
                     var.sem=sem(!!var))
}

plot.val.change.between.exp = function(df, var) {
  var = enquo(var)
  summary.df = create.animal.summary(df, !!var)
  summary.df %>%
    ggplot() +
    #geom_point(aes(x=exp, y=var.mean, color=animal)) +
    geom_line(aes(x=exp, y=var.mean, group=animal)) +
    geom_ribbon(aes(x=exp, 
                    ymin=var.mean - var.sem, 
                    ymax=var.mean + var.sem,
                    group=animal), 
                alpha=0.1) +
    facet_grid(. ~ implant, scales = 'free') +
    xlab('') +
    gtheme
}
