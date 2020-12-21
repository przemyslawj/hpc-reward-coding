library(dplyr)
library(ggplot2)
library(rlang)
library(lmerTest)


plot.model.diagnostics = function(m, animals, effect_vars) {
  df = data.frame(res=unname(residuals(m)), 
                  fitted=unname(fitted(m)), 
                  animal=animals, 
                  effect_var=effect_vars) 
  
  ggplot(df) +
    geom_point(aes(x=res, y=fitted, colour=animal, shape=effect_var)) +
    theme(legend.position = 'none') -> g1
  
  X = qqnorm(residuals(m), plot.it=FALSE)
  data.frame(x=X$x, y=X$y, animal=animals, effect_var=effect_vars) %>%
    ggplot() +
    geom_point(aes(x=x, y=y, colour=animal, shape=effect_var)) +
    stat_qq_line(mapping=aes(sample=res), data=df) +
    xlab('Theoretical quantiles') +
    ylab('Sample quantiles') +
    theme(legend.position = 'right') -> g2
  
  plot_grid(g1,g2)
}

create.animal.summary = function(df, var, ...) {
  var = enquo(var)
  group.vars = enquos(...)
  group_by(df, implant, animal, !!!group.vars) %>% 
    dplyr::summarise(var.mean=mean(!!var, na.rm=TRUE),
                     var.sem=sem(!!var))
}

plot.val.change.between.groups = function(summary.df, var, group.var) {
  var = enquo(var)
  group.var = enquo(group.var)
  summary.df %>%
    ggplot() +
    geom_point(aes(x=!!group.var, y=var.mean, color=animal)) +
    geom_line(aes(x=!!group.var, y=var.mean, group=animal)) +
    geom_ribbon(aes(x=!!group.var, 
                    ymin=var.mean - var.sem, 
                    ymax=var.mean + var.sem,
                    group=animal), 
                alpha=0.1) +
    facet_grid(. ~ implant, scales = 'free') +
    xlab('') +
    gtheme
}


lmer.test.print = function(df, 
                           var = log(1.0 + zscored_deconv_trace.mean), 
                           fixed.effects = is.early.learning,
                           randef.str = '(1 + day_desc | animal)',
                           diagnostics.groupvar = is.early.learning,
                           print.mean.stats = FALSE) {
  var = enexpr(var)
  fixed.effects = enexpr(fixed.effects)
  diagnostics.groupvar = enquo(diagnostics.groupvar)
  
  form = glue::glue('{varname} ~ {fixed.effects} + {random.effects}',           
                    varname=quo_name(var),                                           
                    fixed.effects=quo_name(fixed.effects),
                    random.effects=randef.str)
  model = lmerTest::lmer(form,
                         data=df,
                         REML=TRUE)
  
  fixed_varvals = mutate(df, x = !!diagnostics.groupvar) %>% pull(x)
  g = plot.model.diagnostics(model, df$animal, fixed_varvals)
  print(g)
  
  if (print.mean.stats) {                                                       
    print('Mean stats')                                                         
    group_by(df, !!fixed.effects) %>%                      
      dplyr::summarise(mean(!!var), sem(!!var), n()) %>% print                  
  } 
  print(summary(model))
  print(anova(model, refit=FALSE, ddf='Satterthwaite'))
  return(model)
}

pairwise.post.hoc = function(m.full, factor.interaction=c('implant:is.early.learning'), p.adjust.method='holm') {
  x = lsmeansLT(m.full, pairwise=TRUE, which=factor.interaction, ddf='Satterthwaite')
  x$`Pr(>|t|).adj` = p.adjust(x$`Pr(>|t|)`, method=p.adjust.method)
  x                                                                             
}
