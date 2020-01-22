list.subdir = function(subdir) {
  learning_datestr = list.files(subdir)
  sapply(learning_datestr, FUN=function(x) {paste(subdir, x, sep='/')}) %>% unname
}

norm2 = function(x, y) {
  sqrt(x**2 + y**2)
}


get.subject.result.dirs = function(root_dir, subject) {
  dated_dirs = c(list.subdir(paste(root_dir, 'habituation', sep='/')),
                 list.subdir(paste(root_dir, 'learning', sep='/')))
  sapply(dated_dirs, FUN=function(x) {paste0(x, '/caiman/', subject, '/filtered/')})
}


max.consecutive.days = function(dates) {
  ordered.dates = ordered(dates) %>% as.Date
  days.gap = diff(c(as.Date('2018-01-01'), ordered.dates, as.Date('2100-01-01')))
  max(diff(which(days.gap != 1)))
}


get.trial.ends = function(timestamps) {
  c(which(diff(timestamps) < 0), length(timestamps))
}

pairdiff = function(x, empty.val=0) {
  if (length(x) == 0) {
    return(0)
  }
  if (length(x) == 1) {
    x=c(x, empty.val)
  }
  diff(x)[[1]]
}

pairequal = function(x) {
  if (length(x) < 2) {
    return(TRUE)
  }
  # make two NAs return TRUE
  x[is.na(x)] = -1
  return(x[[1]] == x[[2]])
}

read.mouse.meta = function(rootdirs) {
  mouseinfo.filename = 'mouse_info.csv'
  mouse.infos = lapply(rootdirs, function(d) { 
    tryCatch(
      read.csv(file.path(d, mouseinfo.filename)),
      error = function(e) {
        warning(paste('Can not read trials meta for file', d), e)
        data.frame()
      }
    )
  })
  mouse.info.df = do.call('rbind', mouse.infos) %>%
    dplyr::rename('animal'='mouse') 
  return(mouse.info.df)
}

create.partial.df = function(fieldList, occupancyList, animal_name, day, cell_name, group_name) {
  field = fieldList[[animal_name]][[day]][[paste(cell_name)]]
  if (is.null(field)) {
    #print(paste('Field not found for', animal_name, day, cell_name))
    return(data.frame())
  }
  df = create.pf.df(field,
                    occupancyList[[animal_name]][[day]][[paste(cell_name)]])
  if (nrow(df) > 0) {
    df$group = group_name
  }
  return(df)
}

