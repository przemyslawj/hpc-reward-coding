zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
}

df = function(f) {
  f0 = mean(f)
  (f - f0) / f0
}