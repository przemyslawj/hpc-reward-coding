library(dplyr)

sim.xy.coord = function(nplace.fields=100) {
  set.seed(3456)
  r = runif(nplace.fields)
  set.seed(3457)
  alpha = runif(nplace.fields, 0, 2)
  field.x = (r * cospi(alpha) / 2 + 0.5) * 100
  field.y = (r * sinpi(alpha) / 2 + 0.5) * 100
  
  list(field.x=field.x, field.y=field.y)
}

calc.rew.distances.df = function(loc.df, field.x, field.y) {
  location.sets.df = dplyr::select(loc.df, -date) %>%
    dplyr::distinct()
  
  map_dfr(1:nrow(location.sets.df), ~ data.frame(
    animal=location.sets.df$animal[.x],
    location_set=location.sets.df$location_set[.x],
    min.rew.dist=calc.min.rew.dist(loc.df[.x, ], field.x, field.y)$rew.dist)
  )
}
