vec_dist = function(x, y) {
  dist_x = c(0, diff(x))
  dist_y = c(0, diff(y))
  return(sqrt(dist_x^2 + dist_y^2))
}
