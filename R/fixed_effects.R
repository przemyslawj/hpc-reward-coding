library(dplyr)
library(ggplot2)
library(rlang)
library(lmerTest)
library(BayesFactor)


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

create.summary = function(df, var, ...) {
  var = enquo(var)
  group.vars = enquos(...)
  group_by(df, !!!group.vars) %>% 
    dplyr::summarise(var.mean=mean(!!var, na.rm=TRUE),
                     var.median=median(!!var, na.rm=TRUE),
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
  
  df = dplyr::mutate(df, yvar = !!var)
  form = glue::glue('yvar ~ {fixed.effects} + {random.effects}',           
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

get.effect.size.t.test = function(t.test.res, estimated.val, p.val.col='Pr(>|t|).adj') {
  print(t.test.res[, c('df', 't value', p.val.col)])
  sprintf('%s, estimate change by %.0f pct (%.0f; %.0f), ', 
          row.names(t.test.res),
          t.test.res$Estimate / estimated.val * 100,
          t.test.res$lower / estimated.val * 100,
          t.test.res$upper / estimated.val * 100)
}

#############################################
### Bayes factor and credibility intervals  #
#############################################
# Bayes factor from Bayesian Information Criterion, formula (10) from http://www.ejwagenmakers.com/2007/pValueProblems.pdf
show.bayes.factor = function(model, null.model) {
  BF_BIC = exp((BIC(null.model) - BIC(model))/2)
  #BF_BIC = exp((BIC(model) - BIC(null.model))/2)
  print(BF_BIC)
}

create.bayes.lm.pair = function(df, 
                                formula.full = val ~ 1 + implant + aday,
                                formula.null = val ~ 1 + aday,
                                whichRandom='aday',
                                rscaleRandom='medium',
                                rscaleFixed='wide',
                                iterations=1000) {
  m.full = BayesFactor::lmBF(formula.full, data=df,
                             whichRandom = whichRandom, 
                             rscaleRandom = rscaleRandom, 
                             rscaleFixed = rscaleFixed,
                             iterations = iterations)
  m.null = BayesFactor::lmBF(formula.null, data=df,
                             whichRandom = whichRandom, 
                             rscaleRandom = rscaleRandom, 
                             rscaleFixed = rscaleFixed,
                             iterations = iterations)
  list(full=m.full, null=m.null)
}


calc.pair.95CI = function(full.model, 
                          ytransform=function(x) {x}, 
                          show.percent.change=TRUE,
                          pair.vars = c('implant-dCA1', 'implant-vCA1')) {
  samples = BayesFactor::posterior(full.model, iterations = 10000, columnFilter="^aday$")
  samples.dca1 = ytransform(samples[,'mu'] + samples[, pair.vars[1]])
  samples.vca1 = ytransform(samples[,'mu'] + samples[, pair.vars[2]])
  implant.diff.samples = samples.vca1 - samples.dca1
  out.format = 'Mean diff %.3f CI=[%.3f, %.3f]'
  if (show.percent.change) {
    out.format = 'Mean diff %.f pct CI=[%.f, %.f]'
    implant.diff.samples = samples.vca1 / samples.dca1 * 100 - 100
  } 
  conf.inter = quantile(implant.diff.samples, c(0.025, 0.975))
  sprintf(out.format,
          mean(implant.diff.samples),
          conf.inter[1], conf.inter[2]) %>% print
}
