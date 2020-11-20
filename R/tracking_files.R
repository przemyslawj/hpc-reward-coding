library(dplyr)

source('utils.R')


get_tracking_files = function(root_dat_dir) {
  print(paste0('Getting tracking files for dir=', root_dat_dir))
  tracking_files = tibble(filepath = character(), 
                          filename = character(),
                          date = character(),
                          animal = factor(),
                          trial = integer(),
                          exp_title = character())
  
  all_subfiles = list.files(root_dat_dir, full.names=TRUE)
  dated_subdirs = all_subfiles[file.info(all_subfiles)$isdir]
  for (dated_dir in dated_subdirs) {
    for (exp_title in list.files(dated_dir)) {
      
      if (exp_title == 'caiman') {
        next
      }
      
      tracking_files_list = list.files(file.path(dated_dir, exp_title), 
                                       pattern='*_positions.csv', 
                                       recursive=TRUE,
                                       full.names = TRUE)

      for (filepath in tracking_files_list) {
        filename = basename(filepath)
        fileparts = strsplit(filename,'_')[[1]]
        date_str = basename(dated_dir)
        
        animal = fileparts[length(fileparts)-3]
        trial_n = fileparts[length(fileparts)-1]
        if (endsWith(date_str, '_test')) {
          date_str = substring(date_str, 1, nchar(date_str) - nchar('_test'))
          if (animal == 'test') {
            animal = fileparts[3]
            trial_n = fileparts[5]
          }
        } 
        tracking_files = add_row(tracking_files, 
                                 filepath=filepath,
                                 date=date_str,
                                 filename=basename(filename), 
                                 animal=animal, 
                                 trial=as.integer(trial_n), 
                                 exp_title=exp_title)
      }
    }
  }
  tracking_files$date = char2date(tracking_files$date)
  return(tracking_files)
}