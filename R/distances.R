vec_dist = function(x, y) {
  dist_x = c(0, diff(x))
  dist_y = c(0, diff(y))
  return(sqrt(dist_x^2 + dist_y^2))
}

get_velocity = function(vec_dist, timestamps) {
  t_diff_sec = c(1, diff(timestamps)) / 1000
  vec_dist / t_diff_sec
}

calc.min.rew.dist = function(locs, field.x, field.y) {
  if (nrow(locs) == 0) {
    nrows.out = length(field.x)
    return(list(location.ordinal=rep(NA, nrows.out),
                rew.dist=rep(NA, nrows.out)))
  }
  
  dist.vals = c()
  for (rew_i in 1:nrow(locs)) {
    dist.vals = cbind(
      dist.vals,
      norm2(locs$trans_x[rew_i] - field.x, locs$trans_y[rew_i] - field.y))
  }
  
  min.i = apply(dist.vals, 1, which.min)
  min.dist = apply(dist.vals, 1, min)
  return(list(location.ordinal=locs$location_ordinal[min.i],
              rew.dist=min.dist))
}
