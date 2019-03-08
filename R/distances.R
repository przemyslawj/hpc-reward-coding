vec_dist = function(x, y) {
  dist_x = c(0, diff(x))
  dist_y = c(0, diff(y))
  return(sqrt(dist_x^2 + dist_y^2))
}

get_velocity = function(vec_dist, timestamps) {
  t_diff_sec = c(1, diff(timestamps)) / 1000
  vec_dist / t_diff_sec
}