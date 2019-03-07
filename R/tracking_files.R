library(dplyr)

get_tracking_files = function(root_dat_dir) {
  
  tracking_files = tibble(filepath = character(), 
                          filename = character(),
                          date = character(),
                          animal = factor(),
                          trial = integer(),
                          is_test = logical())
  
  all_subfiles = list.files(root_dat_dir, full.names=TRUE)
  dated_subdirs = all_subfiles[file.info(all_subfiles)$isdir]
  for (dated_dir in dated_subdirs) {
    tracking_dat_dir = paste0(dated_dir, '/movie/tracking/')
    
    for (filename in list.files(tracking_dat_dir, pattern='*.csv$')) {
      fileparts = strsplit(filename,'_')[[1]]
      date_str = basename(dated_dir)
      
      is_test = FALSE
       animal = fileparts[2]
       trial_n=fileparts[4]
      if (endsWith(date_str, '_test')) {
        date_str = substring(date_str, 1, nchar(date_str) - nchar('_test'))
        is_test = TRUE
        if (animal == 'test') {
          animal = fileparts[3]
          trial_n=fileparts[5]
        }
      } 
      tracking_files = add_row(tracking_files, filepath=paste0(tracking_dat_dir, filename),
                               date=date_str, filename=basename(filename), animal=animal, 
                               trial=as.integer(trial_n), is_test=is_test)
    }
  }
  return(tracking_files)
}