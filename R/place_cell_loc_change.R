px2pt = 1/(ggplot2::.pt*72.27/96)

plot.cluster.assignments = function(clust.df, clust.distance=50, ...) {
  xvars = enquos(...)
  xvar_names = lapply(xvars, quo_name)
  
  max.ncells = unique(clust.df$cell_id) %>% max
  clust.df$animal = as.factor(clust.df$animal)
  clust.df = dplyr::mutate(clust.df,
                           cell_id=(as.numeric(animal) - 1) * max.ncells + cell_id)
  
  my.clust.df = dplyr::select(clust.df, cell_id, exp, my.clust) %>%
    spread(exp, my.clust)
  if (! (all(levels(clust.df$exp) %in% colnames(my.clust.df)))) {
    return(ggplot())
  }
  
  get.pattern.desc = function(clust_char) {
    # if (clust_char %in% c('22222')) {
    #   return('stable_place_cell')
    # }
    # 
    # if (clust_char %in% c('21222', '21122', '21221', '21121')) {
    #   return('reward_silenced_cell')
    # }
    
    if (clust_char %in% c('12111', '12211')) {
      return('1-reward_cell')
    }
    # 
    # if (clust_char %in% c('11211', '11221', '11222', '11121', '11122')) {
    #   return('active_after_transloc')
    # }

    return('other')
  }
  
  ncells = unique(clust.df$cell_id) %>% length
  my.clust.df = dplyr::arrange(my.clust.df, !!!xvars) 
  my.clust.df$ordinal = 1:nrow(my.clust.df)
  
  assignment.groups = my.clust.df %>% 
    dplyr::group_by(!!!xvars) %>%
    dplyr::summarise(min.ordinal = min(ordinal),
                     max.ordinal = max(ordinal),
                     pct = round((max.ordinal - min.ordinal + 1)/ ncells * 100)) %>%
    dplyr::mutate(clust_char = paste0(!!!xvars)) %>%
    dplyr::mutate(group_desc = map_chr(clust_char, get.pattern.desc))
  
  assignment.groups$group_ordinal = 1:nrow(assignment.groups)
  
  assignment.points = gather(assignment.groups, key='exp', value='clust', 
                             unlist(xvar_names))
  
  assignment.points = assignment.points %>%
    arrange(clust, group_ordinal) %>%
    ddply(.(exp), mutate,
          cum.pct=cumsum(pct))
  
  assignment.points$group_ordinal = as.factor(assignment.points$group_ordinal)
  assignment.points$exp = factor(assignment.points$exp, levels=unlist(xvar_names))
  assignments.polypoints = assignment.points %>%
    mutate(cum.pct.start = cum.pct - pct) %>%
    gather(cum.pct, cum.pct.start, key='src', value='y') %>%
    arrange(clust_char, exp, src)
  
  assignments.polypoints$exp = as.numeric(assignments.polypoints$exp)
  assignments.polypoints = bind_rows(
    dplyr::mutate(assignments.polypoints, exp=exp-0.05),
    dplyr::mutate(assignments.polypoints, exp=exp+0.05)
  )
  ordered.polypoints = bind_rows(
    filter(assignments.polypoints, src=='cum.pct.start') %>% arrange(exp),
    filter(assignments.polypoints, src=='cum.pct') %>% arrange(desc(exp)))
  
  yintercept = max((filter(ordered.polypoints, clust==1))$y) + clust.distance/2
  ordered.polypoints %>%
    ggplot() +
    geom_hline(aes(yintercept = yintercept), linetype='dashed') +
    geom_polygon(aes(x=exp, 
                     y=ifelse(clust==1, y, y + clust.distance), 
                     #y=cum.pct,
                     fill=group_desc,
                     group=group_ordinal), alpha=0.7) +
    scale_fill_grey(start=0.0, end=0.6) +
    scale_x_continuous(breaks=seq_along(xvar_names), labels = xvar_names) +
    theme(axis.line.x = element_blank(),
          line = element_line(size = px2pt * 0.5),
          legend.position = 'top') +
    geom_text(x=1, y=10, label=paste0('ncells=',ncells)) +
    xlab('') + ylab('') 
  
}


selected.peaks.df = bind_rows(
  #filter(field.peaks.habituation, day_desc=='habituation day#1') %>% 
  #  dplyr::mutate(exp='hday1',
  #                active=rew.peaks.count > 0),
  filter(field.peaks.habituation, day_desc=='habituation day#3') %>% 
    dplyr::mutate(exp='hday3',
                  active=rew.peaks.count > 0),
  filter(field.peaks.df, day_desc=='learning1 day#1') %>% 
    dplyr::mutate(exp='l1d1',
                  active=rew.peaks.count > 0),
  filter(field.peaks.df, day_desc=='learning1 day#5') %>% 
    dplyr::mutate(exp='l1d5',
                  active=rew.peaks.count > 0),
  #filter(field.peaks.beforetest, day_desc=='learning2 day#1') %>% 
  #  dplyr::mutate(exp='l2d1_test',
  #                active=rew.peaks.count > 0),
  filter(field.peaks.df, day_desc=='learning2 day#1') %>% 
    dplyr::mutate(exp='l2d1',
                  active=rew.peaks.count > 0),
  #filter(field.peaks.df, day_desc=='learning2 day#2') %>% 
  #  dplyr::mutate(exp='l2d2',
  #                active=rew.peaks.count > 0),
  filter(field.peaks.df, day_desc=='learning3 day#2') %>% 
    dplyr::mutate(exp='l3d2',
                  active=rew.peaks.count > 0)
) %>%
  dplyr::select(animal, day_desc, cell_id, exp, active, peak2rew.mindist)

ndays = selected.peaks.df$exp %>% unique %>% length
selected.peaks.df$active = as.integer(selected.peaks.df$active)
peakatrew.clusters = reshape2::dcast(selected.peaks.df, 
                  animal + cell_id ~ exp, 
                  value.var='active') 
peakatrew.clusters$days.present = apply(peakatrew.clusters[,3:ncol(peakatrew.clusters)], 1, 
                                        FUN=function(vals) {sum(!is.na(vals))})
peakatrew.clusters = peakatrew.clusters %>%
  dplyr::filter(days.present >= ndays) %>% 
  dplyr::select(-days.present) %>%
  replace_na(list(hday1=0, 
                  hday3=0,
                  l1d1=0,
                  l1d5=0, 
                  l2d1_test=0, 
                  l2d1=0, 
                  l2d2=0, 
                  l3d2=0)) %>% 
  reshape2::melt(id.vars=c('animal', 'cell_id'), variable.name='exp', value.name='active')

peakatrew.clusters$my.clust = as.integer(peakatrew.clusters$active) + 1

peakatrew.clusters = left_join(peakatrew.clusters, mouse.meta.df)
plot.cluster.assignments(filter(peakatrew.clusters, 
                                animal=='F-TL'),
                                #implant=='vCA1'),
                         clust.distance = 100,
                         #hday1, 
                         hday3, 
                         l1d1,
                         l1d5, 
                         #l2d1_test, 
                         l2d1, 
                         #l2d2
                         l3d2
                         )

peakatrew.mindist.wide = reshape2::dcast(field.peaks.df, 
                                         animal + cell_id ~ day_desc, 
                                         value.var='peak2rew.mindist') 
peakatrew.mindist.wide$days.present = apply(peakatrew.mindist.wide[,3:ncol(peakatrew.mindist.wide)], 1, 
                                        FUN=function(vals) {sum(!is.na(vals))})
peakatrew.mindist = filter(peakatrew.mindist.wide, days.present >= 6) %>%
  reshape2::melt(peakatrew.mindist.wide, id.vars=c('animal', 'cell_id'),
                 measure.vars=3:(ncol(peakatrew.mindist.wide)-1), variable.name='day_desc', value.name='peak2rew.mindist') 

peakatrew.mindist = left_join(peakatrew.mindist, mouse.meta.df)
peakatrew.mindist$animal = as.factor(peakatrew.mindist$animal)
peakatrew.mindist %>%
#field.peaks.df$animal = as.factor(field.peaks.df$animal)
#field.peaks.df %>%
  #filter(implant=='dCA1') %>%
  filter(animal=='D-BR') %>%
  dplyr::mutate(unique_cell_id=(as.numeric(animal) - 1)*1000 + cell_id) %>%
  ggplot(aes(x=day_desc, y=peak2rew.mindist, group=unique_cell_id)) +
  geom_point(alpha=0.1) +
  #geom_smooth(formula=y ~ poly(x, 2), method='lm', se=FALSE, alpha=0.2) +
  geom_path(alpha=0.2) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('') + ylab('Zone distance (%)')




day.ordinal.rew.changed = 9
field.peaks.df$day_desc = as.factor(field.peaks.df$day_desc)
field.peaks.df = field.peaks.df %>%
  dplyr::mutate(day_ordinal = as.numeric(day_desc),
                day_ordinal = day_ordinal + ifelse(day_ordinal >= day.ordinal.rew.changed, 1, 0),
                next_day_ordinal = day_ordinal + 1)

l2d1.beforetest = filter(field.peaks.beforetest, day_desc=='learning2 day#1') %>%
  dplyr::mutate(day_ordinal=day.ordinal.rew.changed,
                day_desc = 'learning1 day#6 test',
                next_day_ordinal = day_ordinal + 1) %>%
  left_join(mouse.meta.df) %>%
  left_join(dplyr::select(pc.test.df, cell_id, animal, date, signif.si))
joined.test.trial.peaks = bind_rows(field.peaks.df, l2d1.beforetest)

joined.nextday.peaks =  joined.test.trial.peaks %>%
  left_join(joined.test.trial.peaks, 
            by=c('implant'='implant', 'animal'='animal', 'cell_id'='cell_id', 'next_day_ordinal'='day_ordinal'),
            suffix=c('.fst', '.snd')) %>%
  filter(!is.na(peak2rew.mindist.snd))

# How many times does a cell cross from close to far
cell.activity.changes = joined.nextday.peaks %>%
  dplyr::mutate(changed.rew.activity=rew.peaks.count.fst != rew.peaks.count.snd ) %>%
  dplyr::group_by(implant, animal, cell_id) %>%
  dplyr::summarise(nday.pairs=n(),
                   ndays=unique(c(day_desc.fst, day_desc.snd)) %>% length,
                   nchanged.rew.activity=sum(changed.rew.activity),
                   nactive=unique(c(ifelse(rew.peaks.count.fst > 0, day_desc.fst, '-1'),
                                    ifelse(rew.peaks.count.snd > 0, day_desc.snd, '-1'))) %>% 
                            setdiff('-1') %>% length
                   ) 

stable.animal.cells = cell.activity.changes %>%
  filter(animal=='K-BR') %>%
  filter(nchanged.rew.activity / nday.pairs < 0.3)
  #filter(nday.pairs >= 6) %>%
  #filter(nchanged.rew.activity <= 2) 
  
perc2dist = 1.8

joined.nextday.peaks %>%
  #filter(animal=='K-BR') %>%
  filter(implant=='vCA1') %>%
  dplyr::mutate(stable_cell=(cell_id %in% stable.animal.cells$cell_id)) %>%
  filter(signif.si.fst) %>%
  #filter(implant=='vCA1') %>%
  filter(day_desc.snd != 'learning3 day#3') %>%
  ggplot(aes(x=day_desc.fst, y=peak2rew.mindist.fst * perc2dist)) +
  geom_segment(aes(xend=day_desc.snd, yend=peak2rew.mindist.snd * perc2dist),#, color=stable_cell),
               alpha=0.2) +
  #facet_grid(signif.si.fst ~ .) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('') + ylab('Distance (cm)')


get.stage.from.day.desc = function(day_desc) {
  if (stringr::str_starts(day_desc, 'habituation')) {
    return('H1')
  }
  if (stringr::str_starts(day_desc, 'learning1 day#6 test')) {
    return('L1_test')
  }
  if (stringr::str_starts(day_desc, 'learning1')) {
    return('L1')
  }
  if (stringr::str_starts(day_desc, 'learning2')) {
    return('L2')
  }
  if (stringr::str_starts(day_desc, 'learning3')) {
    return('L3')
  }
}

# Group days by stages and show active change schematic on groups
stage.peaks.df = joined.test.trial.peaks %>%
  dplyr::mutate(exp=map_chr(day_desc, get.stage.from.day.desc)) %>%
  group_by(animal, cell_id, exp) %>%
  dplyr::summarise(med.rew.dist=median(peak2rew.mindist, na.rm=TRUE),
                   ndays.stage=n()) %>%
  dplyr::mutate(active = med.rew.dist <= 20)

stage.peaks.wide = reshape2::dcast(stage.peaks.df, 
                                   animal + cell_id ~ exp, 
                                   value.var='active') 
stage.peaks.wide$nstages.present = apply(stage.peaks.wide[,3:ncol(stage.peaks.wide)], 1, 
                                            FUN=function(vals) {sum(!is.na(vals))})

stage.peaks.wide = stage.peaks.wide %>%
  dplyr::filter(nstages.present >= 4) %>% 
  dplyr::select(-nstages.present) %>%
  replace_na(list(H1=0, 
                  L1=0,
                  L1_test=0,
                  L2=0,
                  L3=0)) %>% 
  reshape2::melt(id.vars=c('animal', 'cell_id'), variable.name='exp', value.name='active')

stage.peaks.wide$my.clust = as.integer(stage.peaks.wide$active) + 1

stage.peaks.wide = left_join(stage.peaks.wide, mouse.meta.df)
plot.cluster.assignments(filter(stage.peaks.wide, 
                                #cell_id %in% stable.animal.cells$cell_id,
                                animal=='K-BR'),
                                #implant=='vCA1'),
                         clust.distance = 100,
                         H1,
                         L1,
                         L1_test,
                         L2,
                         L3
) +
  scale_x_continuous(labels=c('habituation','learning 1', 'test trial', 'learning 2', 'learning 3')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#TODOs regenerate the peaks and present the figures..
# name the groups: group_desc
# place_cell, reward cell, reward silenced cell, reward_move activated
# quantify 
# Should use only cells with signif spatial information in the quantification?
