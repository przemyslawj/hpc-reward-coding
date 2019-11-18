max.consecutive.days = function(dates) {
  ordered.dates = ordered(dates) %>% as.Date
  days.gap = diff(c(as.Date('2018-01-01'), ordered.dates, as.Date('2100-01-01')))
  max(diff(which(days.gap != 1)))
}