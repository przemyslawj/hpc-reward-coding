list.subdir = function(subdir) {
  learning_datestr = list.files(subdir)
  sapply(learning_datestr, FUN=function(x) {paste(subdir, x, sep='/')}) %>% unname
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