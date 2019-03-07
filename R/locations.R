library(dplyr)

is.date <- function(x) !is.na(as.Date(x, '%Y-%m-%d'))

read_locations = function(root.data.dir) {
  merged.df = data.frame()
  cheeseboard.map = read.csv(paste0(root.data.dir, '/cheeseboard_map.csv')) %>%
    select(Row_X, Row_Y, trans_x, trans_y)
  
  all_subfiles = list.files(root.data.dir, full.names=TRUE)
  dated_subdirs = all_subfiles[file.info(all_subfiles)$isdir]
  for (dated_dir in dated_subdirs) {
    date_str = basename(dated_dir)
    
    is_test = FALSE
    if (endsWith(date_str, '_test')) {
      date_str = substring(date_str, 1, nchar(date_str) - nchar('_test'))
      is_test = TRUE
    }
    if (is.date(date_str)) {
      print(dated_dir)
      locations.df = read.csv(paste0(dated_dir, '/locations.csv'), stringsAsFactors=TRUE)
      locations.df$date = rep(date_str, nrow(locations.df))
      locations.df$is_test = rep(is_test, nrow(locations.df))
      locations.df.pos = left_join(locations.df, cheeseboard.map, by=c("Well_row"="Row_X", "Well_col"="Row_Y"))
      merged.df = bind_rows(merged.df, locations.df.pos)
    }
  }
 
  
  result.df = data.frame()
  # Add location_set column 
  merged.df$location_set = rep(0, nrow(merged.df))
  for (animal in levels(merged.df$Animal)) {
    animal.locations = filter(merged.df, Animal == animal) %>% 
      arrange(date, desc(is_test))
    location_set = 0
    prev_locs = c()
    prev_date = '0000-00-00'
    for (i in 1:nrow(animal.locations)) {
      loc = paste0(animal.locations$Well_row[i], 'x', animal.locations$Well_col[i])
      if (!(loc %in% prev_locs)) {
        if (prev_date == animal.locations$date[i]) {
          prev_locs = append(prev_locs, loc)
        } else {
          location_set = location_set + 1
          prev_locs = c(loc)
          prev_date = animal.locations$date[i]
        }
      }
      animal.locations$location_set[i] = location_set
      
    }
    
    result.df = bind_rows(result.df, animal.locations)
  }
  
  return(result.df)
}