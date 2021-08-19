eval.decoder.other.days.on.traces = function(binned.traces.filtered,
                                             decoder.config,
                                             test.day.indecies,
                                             equal.prior=FALSE,
                                             filter.train.traces.fun=NULL,
                                             filter.test.traces.fun=NULL) {
  traces.dates = char2date(unique(binned.traces.filtered$date)) %>% sort()
  eval.dfs = data.frame()
  for (i in test.day.indecies) {
    train.traces = binned.traces.filtered[date != as.Date(traces.dates[i]),]
    if (!is.null(filter.train.traces.fun)) {
      train.traces = filter.train.traces.fun(train.traces)
    }
    test.traces = binned.traces.filtered[date == as.Date(traces.dates[i]),]
    print('Test traces size before')
    print(dim(test.traces))
    if (!is.null(filter.test.traces.fun)) {
      test.traces = filter.test.traces.fun(test.traces, train.traces)
    }
    if (nrow(test.traces) == 0) {
      print(sprintf('Test traces empty for test date %s, skipping', traces.dates[i]))
      next
    }
    model = train.decoder.model(train.traces, decoder.config)
    if (equal.prior) {
      model$prior[model$prior > 0] = 1.0 / sum(model$prior > 0)
    }
    eval.res = eval.model(model, test.traces, decoder.config)
    eval.res$df$model_date = str_c(format(traces.dates[-i]), collapse=',')
    eval.res$df$model_id = i
    eval.res$df$ncells_present = length(unique(binned.traces.filtered$cell_id))
    eval.dfs = bind_rows(eval.dfs, eval.res$df)
  }
  return(eval.dfs)
}

crosseval.decoder.on.paired.days = function(binned.traces.filtered,
                                            decoder.config,
                                            test.day.indecies, # excluded from the cross-validation
                                            equal.prior=FALSE,
                                            filter.train.traces.fun=NULL,
                                            filter.test.traces.fun=NULL,
                                            min.samples=5) {
  traces.dates = char2date(unique(binned.traces.filtered$date)) %>% sort()
  eval.dfs = data.frame()
  for (i in test.day.indecies) {
    traces = binned.traces.filtered[date != as.Date(traces.dates[i]),]
    if (!is.null(filter.train.traces.fun)) {
      traces = filter.train.traces.fun(traces)
    }
    print(sprintf('Traces timestamps after filering %d', unique(traces$timestamp) %>% length()))
    print(sprintf('Timestamps for test %d', unique(filter.test.traces.fun(traces)$timestamp) %>% length()))
    if (nrow(traces) == 0) {
      print(sprintf('Train traces empty for test date %s, skipping', traces.dates[i]))
      next
    }
    
    eval.res = eval.decoder(
      traces,
      nclasses=decoder.config$nclasses,
      predicted.var=!!decoder.config$predicted.var,
      error.fun=decoder.config$error.fun,
      train.fun=decoder.config$train.fun,
      predict.fun=decoder.config$predict.fun,
      predictor.var=!!decoder.config$predictor.var,
      cv=TRUE,
      min.samples=min.samples,
      filter.test.traces.fun=filter.test.traces.fun)

    if (nrow(eval.res) > 0) {
      eval.res$model_date = str_c(format(traces.dates[-i]), collapse=',')
      eval.res$model_id = i
      eval.res$ncells_present = length(unique(traces$cell_id))
      eval.res$animal = traces$animal[1]
      eval.dfs = bind_rows(eval.dfs, eval.res)
    }
  }
  return(eval.dfs)
}


eval.decoder.all.other.days = function(decoder.config, 
                                       test.days.df, 
                                       equal.prior=FALSE,
                                       prepare.vars.fun=NULL,
                                       filter.train.traces.fun=NULL,
                                       filter.test.traces.fun=NULL,
                                       test.day.indecies=c(1,2,3),
                                       min.ncells.for.eval=10,
                                       min.present.times=2,
                                       traces.decode.fun=eval.decoder.other.days.on.traces) {
  ndays = 3  # Number of the first beforetest days to use
  eval.dfs = data.frame()
  for (animal_name in unique(test.days.df$animal)) {
    beforetest_dates = filter(test.days.df, animal == animal_name) %>% pull(date)
    animal.caimg_dirs = map_chr(beforetest_dates, ~ find.caimg.dir(caimg_result_dirs, animal_name, .x))[1:ndays]
    data.traces = map_dfr(animal.caimg_dirs, ~ {
      x = read.data.trace(.x, filter_exp_title = 'beforetest')
      x$date = rep(char2date(x$date[1]), nrow(x))
      x
    })
    binned.traces = prepare.traces(data.traces,
                                   filter.running=TRUE,
                                   timebin.dur.msec = 200,
                                   get.bin.thresholds.fun = get.quantiles.fun(c(0.9, 1.0)),
                                   binned.var='zscored_trace')
    prepared.traces = binned.traces
    if (!is.null(prepare.vars.fun)) {
      prepared.traces = prepare.vars.fun(binned.traces)
    }
    cells.present = prepared.traces[, (1), by=.(date, cell_id)][, .(ntimes=.N), by=.(cell_id)][
      ntimes >= min.present.times, cell_id]
    sprintf('%d cells present on at least one day for animal=%s', length(cells.present), animal_name) %>% print
    if (length(cells.present) < min.ncells.for.eval) {
      sprintf('Few cells present across days, skipping eval for animal=%s', animal_name) %>% print
      next
    }
    
    prepared.traces.filtered = prepared.traces[cell_id %in% cells.present]
    eval.res = traces.decode.fun(prepared.traces.filtered, 
                                 decoder.config,
                                 test.day.indecies,
                                 equal.prior,
                                 filter.train.traces.fun,
                                 filter.test.traces.fun)
    eval.dfs = bind_rows(eval.dfs, eval.res)
  }
  return(eval.dfs)
}


avg.pop.activity.debug = function(binned.traces, 
                                  min.present.times=1,
                                  min.cells.per.group=3) {
  min.times_rew_cell = 1
  trace.dates = unique(binned.traces$date) %>% char2date()
  cell.db.daily =  beforetest.field.dist2current.rew[
    (date %in% trace.dates & animal == binned.traces$animal[1]), 
    .(animal, date, cell_id, signif.si)]
  
  cell.db = beforetest.field.dist2current.rew[
    (date %in% trace.dates & animal == binned.traces$animal[1]), 
    .(animal, date, cell_id, 
      is.rew.cell = min.rew.dist <= min.rew.dist.thresh,
      is.signif.rew.cell = signif.si & (min.rew.dist <= min.rew.dist.thresh) )][
        ,.(times_signif_rew_cell = sum(is.signif.rew.cell),
           times_rew_cell = sum(is.rew.cell),
           npresent = .N,
           times_signif_rew_cell_pct = sum(is.signif.rew.cell) / .N),
        .(animal, cell_id)]
  
  # restrict to cells present on all 3 days
  sprintf('#cells at all three beforetests: %d', nrow(cell.db[npresent >= min.present.times, ]) ) %>% print
  sprintf('#signif reward cells amongst them: %d', 
          nrow(cell.db[npresent >= min.present.times & times_signif_rew_cell >= min.times_rew_cell]) ) %>% print
  sprintf('# reward cells amongst them: %d',
          nrow(cell.db[npresent >= min.present.times & times_rew_cell >= 2]) ) %>% print
  joined.data.traces = merge(binned.traces, cell.db, by=c('animal','cell_id'), all.x = TRUE)
  joined.data.traces[, min.rew.dist := pmin(dist_reward0, dist_reward1)]
  joined.data.traces[, is.rew.cell := (times_signif_rew_cell >= min.times_rew_cell)]
  joined.data.traces[, rew.cell.group := pmin(times_signif_rew_cell, 2)]
  #joined.data.traces[, rew.cell.group := as.integer(times_signif_rew_cell_pct >= 0.0)]
  #joined.data.traces[, rew.cell.group := pmin(pmin(1, times_signif_rew_cell) * times_rew_cell, 2)]
  joined.data.traces[, close2rew := as.integer(min.rew.dist <= min.rew.dist.thresh) + 1]  
  
  filtered.cell.groups = joined.data.traces[, .(ncells.in.group=length(unique(cell_id))), by=list(rew.cell.group)][
    ncells.in.group >= min.cells.per.group, rew.cell.group]
  
  pop.traces.rew.cells = joined.data.traces[(npresent >= min.present.times) & 
                                              (rew.cell.group %in% filtered.cell.groups) &
                                              (rew.cell.group >= 0), 
                                            #.(avg.val = mean(smoothed_deconv_trace)),
                                            .(avg.val = mean(response_bin)-1),
                                            by=list(animal, date, exp_title, trial_id, trial, x, y,
                                                    close2rew, min.rew.dist, timestamp, time_bin, rew.cell.group)]
  #pop.traces.rew.cells$cell_id = cell.db[npresent==3, cell_id][1]
  # workaround to use existing cell ids
  # cell.group.cell.id = joined.data.traces[(npresent == 3), .(fst.cell_id=cell_id[1]), by=list(rew.cell.group)]
  # pop.traces.rew.cells = pop.traces.rew.cells[cell.group.cell.id, on=.(rew.cell.group)]
  # pop.traces.rew.cells$cell_id = pop.traces.rew.cells$fst.cell_id
  pop.traces.rew.cells$cell_id = pop.traces.rew.cells$rew.cell.group
  ngroups = length(unique(pop.traces.rew.cells$rew.cell.group))
  sprintf('# cell groups: %d', ngroups) %>% print
  if (ngroups > 0) {
    g = ggplot(pop.traces.rew.cells, aes(x=timestamp/1000)) +
      geom_tile(aes(fill=as.factor(close2rew), y=0.1), height=0.3) +
      geom_line(aes(y=avg.val)) +
      facet_grid(rew.cell.group ~ trial_id) + theme(legend.position = 'none')
    ggsave(paste0('/home/prez/tmp/cheeseboard/', binned.traces$trial_id[1], '.png'))
  }
  pop.traces.rew.cells
}

caret.svm.predict = function(prior, svmFit, predictor.vals, nclasses=2) {
  testing.wide.df = data.frame(t(predictor.vals))
  colnames(testing.wide.df) = paste0('cell_', names(predictor.vals))
  predicted <- predict(svmFit, testing.wide.df)
  predicted.class.index = as.integer(substr(as.character(predicted), 7, 100))
  list(s=predicted.class.index,
       density=rep(1/nclasses, nclasses),
       prob=1)
}
