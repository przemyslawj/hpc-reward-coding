library(dplyr)
library(data.table)
library(datatrace)
library(permute)


get.max.thr = function(max.fraction=0.25, absolute.min=0.01) {
  function(vals) {  
    max.val = max(c(vals, absolute.min), na.rm = TRUE)
    c( max.val* max.fraction, max.val+1) 
  }
}

get.nvenent.thr = function(min.nevents=0.1) {
  function(vals) { c(min.nevents, max(vals) + 1)}
}

prepare.traces = function(data.traces,
                          binned.var,
                          get.bin.thresholds.fun,
                          filter.running=TRUE,
                          timebin.dur.msec=200) {
  setorder(data.traces, exp_title, trial_id, cell_id, timestamp)
  data.traces = detect.events(data.traces, deconv.threshold=0.1)
  # running speed avg > 4 cm/s in 0.5 s window
  data.traces = add.running.col(data.traces, 3.3, 10)
  data.traces = gauss.smooth.df.var(data.traces, filter.len=20, sigma=4.0)
  data.traces[, `:=` (zscored_deconv_trace = zscore(deconv_trace), 
                      zscored_trace = zscore(trace),
                      zscored_smooth_deconv_trace = zscore(smoothed_deconv_trace)),
              by=.(exp_title, cell_id)]
  
  if (filter.running) {
    data.traces.filtered = data.traces[ x > 0 & y > 0 & is_running, ]
  } else {
    data.traces.filtered = data.traces[x > 0 & y > 0,]
  }
  
  date_str = format(data.traces$date[1])
  animal_name = data.traces$animal[1]
  
  binned.traces = bin.time.space(data.traces.filtered,
                                 nbins.x = xybins,
                                 nbins.y = xybins,
                                 get.bin.thresholds.fun = get.bin.thresholds.fun,
                                 binned.var=binned.var,
                                 timebin.dur.msec=timebin.dur.msec)
  binned.traces = gauss.smooth.df.var(binned.traces, var='nevents', out.var='smoothed_nevents',
                                      filter.len=5, sigma=1.0)
  
  return(binned.traces)
}


discrete.bayes.spatial.decoder = list(
  stim.var=quo(bin.xy),
  value.var=quo(response_bin),
  nstim.bins=xybins^2,
  train.fun = create.discrete.bayes,
  predict.fun = bayesmax,
  error.fun=bin.distance.error(xybins, xybins)
)

mfr.bayes.spatial.decoder = list(
  stim.var=quo(bin.xy),
  value.var=quo(smoothed_deconv_trace),
  nstim.bins=xybins^2,
  train.fun = create.mfr.bayes,
  predict.fun = bayesmax_poisson,
  error.fun=bin.distance.error(xybins, xybins)
)


read.binned.traces = function(caimg_result_dir, 
                              exp_title='trial') {
  data.traces = read.data.trace(caimg_result_dir, filter_exp_title = exp_title)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  prepare.traces(data.traces, 
                 binned.var='zscored_trace',
                 get.bin.thresholds.fun = get.quantiles.fun(c(0.9, 1.0)),
                 timebin.dur.msec = 200)
}

eval.decoder.same.day = function(binned.traces,
                                 decoder.config, 
                                 prepare.vars.fun=NULL,
                                 filter.test.traces.fun=NULL,
                                 min.samples=10,
                                 ncells=30, # If ncells == 0, take all cells
                                 nshuffles=20) {
  
  if (nrow(binned.traces) == 0) {
    return(data.frame())
  } 
  
  print('Started evaluating the model')
  if (!is.null(prepare.vars.fun)) {
    binned.traces = prepare.vars.fun(binned.traces)
  }
  cell_ids = binned.traces$cell_id %>% unique
  
  all.shuffles.eval.res = map_dfr(1:nshuffles, function(j) {
    subset.traces = binned.traces
    if (ncells > 0) {
      cell_ids = cell_ids[shuffle(cell_ids)]
      subset.traces = binned.traces[cell_id %in% cell_ids[1:min(ncells, length(cell_ids))],]
    }
    
    eval.df = eval.decoder(subset.traces,
                           nstim.bins=decoder.config$nstim.bins,
                           stim.var=!!decoder.config$stim.var,
                           error.fun=decoder.config$error.fun,
                           train.fun=decoder.config$train.fun,
                           predict.fun=decoder.config$predict.fun,
                           value.var=!!decoder.config$value.var,
                           cv=TRUE,
                           min.samples=min.samples,
                           filter.test.traces.fun=filter.test.traces.fun)
    
    if (nrow(eval.df) > 0) {
      eval.df$animal = binned.traces$animal[1]
      eval.df$date = binned.traces$date[1]
      eval.df$shuffle = j
    }
    eval.df
  })
  
  return(all.shuffles.eval.res)
}


scale.bins2cm = function(decoding.df) {
  decoding.df$error = decoding.df$error * cheeseboard.bin.cm
  if ('random_error' %in% names(decoding.df)) {
    decoding.df$random_error = decoding.df$random_error * cheeseboard.bin.cm
  }
  return(decoding.df)
}

join.meta.dfs = function(df) {
  df = left_join(df, mouse.meta.df, by='animal') %>%
    left_join(trials.meta.df, by=c("animal", "date"))
  df$day_desc = as.factor(df$day_desc)
  return(df)
}


eval.model = function(tested_model, binned.traces, decoder.config) {
  model.visited = which(tested_model$prior > 0)
  binned.traces.filtered = binned.traces[get(quo_name(decoder.config$predicted.var)) %in% model.visited,]
  system.time(eval.res <- eval.testdata2(binned.traces.filtered, 
                                         tested_model, 
                                         predict.fun=decoder.config$predict.fun,
                                         error.fun=decoder.config$error.fun,
                                         predicted.var=!!decoder.config$predicted.var,
                                         nclasses=decoder.config$nclasses,
                                         predictor.var=!!decoder.config$predictor.var))
  if (nrow(eval.res$df) == 0) {
    return(list(df=data.frame()))
  }
  
  random.classifier.res <- eval.testdata2(binned.traces.filtered, 
                                          tested_model, 
                                          predict.fun=random.prior.classifier,
                                          error.fun=decoder.config$error.fun,
                                          predicted.var=!!decoder.config$predicted.var,
                                          nclasses=decoder.config$nclasses,
                                          predictor.var=!!decoder.config$predictor.var)
  
  eval.res$df$random_error = random.classifier.res$df$error
  eval.res$df$test_date = as.Date(binned.traces$date[1])
  eval.res$df$animal = binned.traces$animal[1]
  return(eval.res)
}

train.decoder.model = function(training.traces,
                               decoder.config,
                               filter.training.trace.fun=NULL) {
  if (!is.null(filter.training.trace.fun)) {
    training.traces = filter.training.trace.fun(training.traces)  
  }
  if (nrow(training.traces) == 0) {
    return(NULL)
  }
  model = decoder.config$train.fun(
    training.traces, 
    nclasses=decoder.config$nclasses,
    predictor.var=!!decoder.config$predictor.var,
    predicted.var=!!decoder.config$predicted.var)
  
  model$model_dates = unique(training.traces$date)
  return(model)
}


# Evaluates the decoders against the test data in caimg_result_dir, create and append a model for the current dir
eval.decoder.other.day = function(caimg_result_dir, 
                                  decoder.config,
                                  trials.meta.df,
                                  models,
                                  prepare.vars.fun=NULL,
                                  filter.training.trace.fun=NULL,
                                  filter.test.traces.fun=NULL,
                                  exp_title='trial',
                                  max.days.diff=2,
                                  min.days.diff=1,
                                  min.samples=10,
                                  binned.varname='zscored_trace',
                                  nshuffles=1,
                                  ncells=0) {
  data.traces = read.data.trace(caimg_result_dir, filter_exp_title = exp_title)
  data.traces$date = rep(char2date(data.traces$date[1]), nrow(data.traces))
  date_str = format(char2date(data.traces$date[1]))
  animal_name = data.traces$animal[1]
  
  eval.dfs = data.frame()
  binned.traces = prepare.traces(data.traces,
                                 filter.running=TRUE,
                                 timebin.dur.msec = 200,
                                 #get.bin.thresholds.fun = get.nvenent.thr(0.1),
                                 get.bin.thresholds.fun = get.quantiles.fun(c(0.9, 1.0)),
                                 binned.var=binned.varname)
  
  if (!is.null(prepare.vars.fun)) {
    binned.traces = prepare.vars.fun(binned.traces)
  }
  
  for (model_date in names(models[[animal_name]])) {
    model.date.ordinal = filter(trials.meta.df, animal==animal_name, date==model_date) %>% pull(exp_day_ordinal)
    test.date.ordinal = filter(trials.meta.df, animal==animal_name, date==date_str) %>% pull(exp_day_ordinal)
    exp.days.diff = test.date.ordinal - model.date.ordinal
    if (exp.days.diff <= max.days.diff && exp.days.diff >= min.days.diff) {
      tested_model = models[[animal_name]][[model_date]]
      model.cell_ids = as.integer(colnames(tested_model$likelihood))
      cell_ids = intersect(unique(binned.traces$cell_id), model.cell_ids)
      if (length(cell_ids) < 5) {
        warning('Few overlapping cells (', length(cell_ids), ') for animal=', animal_name, 
                ' between day=', format(model_date), ' and day=', date_str)
        next
      }
      test.traces = binned.traces
      if (!is.null(filter.test.traces.fun)) {
        test.traces = filter.test.traces.fun(test.traces, model_date)
      }
      for (shuffle_i in 1:nshuffles) {
        shuffle_cell_ids = cell_ids
        if (ncells > 0) {
          shuffle_cell_ids = shuffle(cell_ids)[1:min(ncells, length(cell_ids))]
        }
        subset.model = tested_model
        subset.model$likelihood = tested_model$likelihood[, which(model.cell_ids %in% cell_ids),]
        eval.res = eval.model(subset.model, binned.traces[cell_id %in% cell_ids], decoder.config)
        eval.res$df$model_date = char2date(model_date)
        eval.res$df$shuffle_i = shuffle_i
        eval.dfs = bind_rows(eval.dfs, eval.res$df)
      }
    }
  }
  
  filtered.dfs = filter.sampled.during.training(binned.traces, binned.traces[1,], min.samples=min.samples)
  model = train.decoder.model(filtered.dfs$train, decoder.config, filter.training.trace.fun=filter.training.trace.fun)
  models[[animal_name]][[date_str]] = model
  return(list(df=eval.dfs, models=models))
}

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
