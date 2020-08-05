library(purrr)

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
                rew.dist=rep(NA, nrows.out),
                rew.angle=rep(NA, nrows.out)))
  }
  
  if (length(field.x) != length(field.y)) {
    warning('Vectors with field x and y coordinates should have equal length')
  }
  
  dist.vals = matrix(data=0, nrow=length(field.x), ncol=nrow(locs))
  angle.vals = matrix(data=0, nrow=length(field.x), ncol=nrow(locs))
  for (rew_i in 1:nrow(locs)) {
    z = complex(real=locs$trans_x[rew_i] - field.x, 
                imaginary=locs$trans_y[rew_i] - field.y)
    dist.vals[, rew_i] = Mod(z)
    angle.vals[, rew_i] = (Arg(z) / pi * 180) %% 360
  }
  
  min.i = apply(dist.vals, 1, which.min)
  min.dist = apply(dist.vals, 1, min)
  min.angle = purrr::map_dbl(1:nrow(angle.vals), ~ angle.vals[.x, min.i[.x]])
  
  return(list(location.ordinal=locs$location_ordinal[min.i],
              rew.dist=min.dist,
              rew.angle=min.angle))
}
