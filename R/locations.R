library(dplyr)
library(purrr)
library(stringr)

source('utils.R')
perc2dist = 1.2
goal.cell.max.dist = 20 / perc2dist
rew_zone_radius = goal.cell.max.dist

is.date <- function(x) !is.na(char2date(x))

get.subdirs = function(path) {
  subdirs = list.files(path, full.names=TRUE)
  return(subdirs[file.info(subdirs)$isdir])
}

add_location_set = function(merged.df) {
  result.df = data.frame()
  
  if (nrow(merged.df) == 0) {
    return(merged.df)
  }
  
  # Add location_set column 
  merged.df$location_set = rep(0, nrow(merged.df))
  merged.df$exp_title = stringr::str_replace(merged.df$exp_title, '^test$', 'beforetest')
  merged.df$exp_title = factor(merged.df$exp_title,
                               levels=c('beforetest', 'trial', 'aftertest'))
  for (animal_name in levels(merged.df$animal)) {
    animal.locations = filter(merged.df, .data$animal == animal_name) %>% 
      arrange(date, exp_title)
    location_set = 0
    animal_locs = c()
    prev_locs = c()
    prev_date = '0000-00-00'
    for (i in 1:nrow(animal.locations)) {
      loc = paste0(animal.locations$Well_row[i], 'x', animal.locations$Well_col[i])
      if (!(loc %in% prev_locs) && animal.locations$Valence[i] == 'Positive') {
        if (prev_date == animal.locations$date[i]) {
          prev_locs = append(prev_locs, loc)
        } else {
          location_set = location_set + 1
          prev_locs = c(prev_locs, loc)
          prev_date = animal.locations$date[i]
        }
      }
      if (!loc %in% animal_locs) {
        animal_locs = append(animal_locs, loc)
      }
      
      animal.locations$location_set[i] = location_set
      animal.locations$location_ordinal[i] = which(loc == animal_locs)
      
    }
    
    result.df = bind_rows(result.df, animal.locations)
  }
  
  # Make location set consistent in one day: if one reward changed then update
  # the location set for all rewards
  result.df = group_by(result.df, animal, date, exp_title, is_test) %>%
    dplyr::mutate(location_set = max(location_set)) %>%
    ungroup() %>%
    arrange(animal, date, exp_title)
  return(result.df)
}


add_prev_locations = function(locations.df, prev.loc.set.diff=1) {
  locations.df = locations.df %>%
    dplyr::mutate(prev_loc_set = location_set - prev.loc.set.diff) %>%
    dplyr::group_by(animal, date, Valence, is_test) %>%
    dplyr::mutate(position_no = row_number())
  
  if (!('current_loc' %in% colnames(locations.df))) {
    locations.df$current_loc = TRUE
    locations.df$previous_loc = FALSE
  }
  
  set.locations.df = filter(locations.df, current_loc) %>%
    group_by(animal, is_test, exp_title, location_set, Valence, position_no) %>%
    dplyr::select(-date, -prev_loc_set) %>%
    dplyr::distinct() %>%
    dplyr::mutate(current_loc = FALSE, previous_loc = TRUE)
  
  previous.locations.df = filter(locations.df, current_loc) %>%
    ungroup() %>%
    dplyr::select(date, location_set, prev_loc_set, animal, exp_title, Valence, position_no) %>%
    left_join(set.locations.df, by=c('prev_loc_set'='location_set',
                                     'animal'='animal', 'exp_title'='exp_title',
                                     'Valence'='Valence', 'position_no'='position_no')) %>%
    filter(prev_loc_set > 0, !is.na(Well_row)) 
  
  joined.locations.df = bind_rows(locations.df, previous.locations.df) %>%
    dplyr::distinct(animal, date, is_test, Valence, location_set, location_ordinal, .keep_all=TRUE)
  return (joined.locations.df)
}

add_future_neg_locations = function(locations.df) {
  locations.df = mutate(locations.df, future_loc = FALSE)
  # Negative locations for each set
  set.locations.df = filter(locations.df, Valence == 'Negative', current_loc == TRUE) %>%
    group_by(animal, is_test, location_set, position_no) %>%
    slice(1) %>%
    dplyr::select(-date, -prev_loc_set) %>%
    rename(next_location_set = location_set) %>%
    mutate(future_loc = TRUE, current_loc = FALSE)
  
  next.locations.df = locations.df %>%
    dplyr::select(date, location_set, prev_loc_set, animal, is_test, Valence, position_no) %>%
    left_join(set.locations.df, by=c('location_set'='next_location_set',
                                     'animal'='animal', 'is_test'='is_test',
                                     'position_no'='position_no'),
              suffix=c('.x', '')) %>%
    group_by(animal, date, is_test, location_set) %>% 
    slice(1)
  
  joined.locations.df = bind_rows(locations.df, next.locations.df) %>%
    group_by(animal, date, is_test, location_set, Valence, position_no, previous_loc) %>%
    arrange(future_loc) %>%
    slice(1) %>%
    dplyr::select(-Valence.x)
  return (joined.locations.df)
}


read_locations = function(root.data.dir) {
  cheeseboard.map = read.csv(file.path(root.data.dir, 'cheeseboard_map.csv')) %>%
    dplyr::select(Row_X, Row_Y, trans_x, trans_y)
  
  expdirs = get.subdirs(root.data.dir)
  merged.df = data.frame()
  for (expdir in expdirs) {
    dated_subdirs = get.subdirs(expdir)
    
    for (dated_dir in dated_subdirs) {
      date_str = basename(dated_dir)
      
      exp_titles = get.subdirs(dated_dir)
      for (exp_titledir in exp_titles) {
        is_test = FALSE
        if (endsWith(exp_titledir, 'test')) {
          #date_str = substring(date_str, 1, nchar(date_str) - nchar('_test'))
          is_test = TRUE
        }
        if (is.date(date_str)) {
          fpath = file.path(exp_titledir, 'locations.csv')
          #print(paste('Reading locations from: ', fpath))
          if (file.exists(fpath)) {
            locations.df = read.csv(fpath, stringsAsFactors=TRUE)
            locations.df$date = rep(date_str, nrow(locations.df))
            locations.df$is_test = rep(is_test, nrow(locations.df))
            locations.df$exp_title = basename(exp_titledir)
            locations.df.pos = left_join(locations.df, cheeseboard.map, by=c("Well_row"="Row_X", "Well_col"="Row_Y"))
            merged.df = bind_rows(merged.df, locations.df.pos)
          }
        }
      }
    }
  }
  
  if (nrow(merged.df) == 0) {
    return(merged.df)
  }
  
  merged.df = dplyr::rename(merged.df, animal=Animal)
  merged.df$animal = as.factor(merged.df$animal)
  result = add_location_set(merged.df)
  result$date = char2date(result$date)
  
  dplyr::distinct(result, animal, Valence, date, is_test, location_set, location_ordinal, .keep_all=TRUE)
}


read.trials.meta = function(rootdirs) {
  df = data.frame()
  for (rootdir in rootdirs) {
    habit.days = list.files(file.path(rootdir, 'habituation'))
    habit.dates = char2date(habit.days) %>% sort()
    habit.df.parts = lapply(habit.days, function(.x) {
      caiman.dir = file.path(rootdir, 'habituation', .x, 'caiman')
      if (!file.exists(caiman.dir)) {
        return(data.frame())
      }
      animal_names = list.files(caiman.dir)
      df = data.frame(animal = animal_names)
      df$date = char2date(.x)
      df$day_ordinal = which(habit.dates == char2date(.x))
      df
    })
    habit.df = do.call('rbind', habit.df.parts)
    habit.df$exp = 'habituation'
    habit.df$location_set = 0
    df = bind_rows(df, habit.df)
    
    rewards.df = read_locations(rootdir) 
    
    if (nrow(rewards.df) > 0) {
      rewards.df = rewards.df %>% 
        group_by(date, animal) %>%
        dplyr::arrange(is_test) %>%
        top_n(1) %>%
        ungroup()
      #dplyr::filter(!is_test)
      
      rewards.df = rewards.df %>%
        dplyr::rename('animal' = 'animal') %>%
        dplyr::select('animal', 'date', 'location_set') %>%
        distinct()
      
      .dayorder = function(days) {
        days = char2date(unique(days)) %>% sort
        map_int(days, ~ which(days == .x))
      }
      
      learning.df = rewards.df %>% 
        dplyr::group_by(animal, location_set) %>%
        dplyr::mutate(day_ordinal = .dayorder(date),
                      exp = 'learning')
    }
    
    df = bind_rows(df, learning.df)
  }
  
  df = mutate(df, day_desc = paste0(
    ifelse(exp == 'habituation', exp, paste0(exp, location_set)), 
    ' day#', 
    day_ordinal))
  
  df = df %>% 
    group_by(animal) %>%
    arrange(date) %>%
    dplyr::mutate(exp_day_ordinal = row_number())
  df
}

filter.rews.df = function(locs.df, day, animal_name) {
  locs.df = as.data.table(locs.df)
  locs = locs.df[animal == animal_name, ]
  if (!is.null(day) && ('date' %in% colnames(locs))) {
    locs = locs[date==format(day) ,]  
  }
  
  return(locs)
}
