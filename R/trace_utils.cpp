#include <Rcpp.h>
using namespace Rcpp;

/***
 * Requires DF to be ordered by trial_id, cell_id and timestamp, so the ajacent rows have
 * sequential timestamps.
 */
// [[Rcpp::export]]
LogicalVector isRunning(DataFrame& df,
                        double min_run_velocity,
                        double mean_run_velocity,
                        double window_dur_ms) {

  LogicalVector result(df.nrows());
  
  NumericVector velocity = df["velocity"];
  NumericVector pos_x = df["x"];
  NumericVector pos_y = df["y"];
  IntegerVector timestamp = df["timestamp"];
  CharacterVector trial_id = df["trial_id"];
  
  int i = 0;
  while (i < df.nrows()) {
    int j = i + 1;
    while(j < df.nrows() && 
          velocity[i] >= min_run_velocity && 
          trial_id[i] == trial_id[j] &&
          velocity[j] >= min_run_velocity) {
      j = j + 1;
    }
    j = j - 1;
    
    double distx = pos_x[j] - pos_x[i];
    double disty = pos_y[j] - pos_y[i];
    double dist = sqrt(distx * distx + disty * disty);
    double dur_ms = std::max(timestamp[j] - timestamp[i], 1);
    double vel = dist / dur_ms * 1000;
    bool is_running = (vel >= mean_run_velocity) && (dur_ms >= window_dur_ms);
    for (int k = i; k <= j; ++k) {
      result[k] = is_running;
    }
    i = j + 1;
  }
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
