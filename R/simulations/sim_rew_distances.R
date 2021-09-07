library(dplyr)
library(data.table)

sim.xy.coord = function(nplace.fields=100, max.val=100, fixed.seed=FALSE) {
  if (fixed.seed) {
    set.seed(3456)
  }
  r = runif(nplace.fields)
  if (fixed.seed) {
    set.seed(3457)
  }
  alpha = runif(nplace.fields, 0, 2)
  field.x = (r * cospi(alpha) / 2 + 0.5) * max.val
  field.y = (r * sinpi(alpha) / 2 + 0.5) * max.val
  
  list(field.x=field.x, field.y=field.y)
}

# Restricts the simulated fields to occupied bins
sim.xy.coord.occupancy_mat = function(is_occupancied_mat, nplace.fields, search.multiplier=4) {
  nbins = dim(is_occupancied_mat)[1]
  random.coors = sim.xy.coord(nplace.fields = search.multiplier * nplace.fields, max.val = nbins)
  field.scalar = 100 / nbins
  kept.random.coors.index = Filter(
    function(i) { is_occupancied_mat[random.coors$field.x[i], random.coors$field.y[i]]},
    1:length(random.coors$field.x))
  if (length(kept.random.coors.index) >= nplace.fields) {
    kept.random.coors = list(field.x=random.coors$field.x[kept.random.coors.index],
                             field.y=random.coors$field.y[kept.random.coors.index])
    return(list(field.x=kept.random.coors$field.x[1:nplace.fields] * field.scalar,
                field.y=kept.random.coors$field.y[1:nplace.fields] * field.scalar))
  }
  
  sim.xy.coord.occupancy_mat(is_occupancied_mat, nplace.fields, search.multiplier * 4)
}

# Create DF with simulated coordinates of fields per animal per day using day occupancy matrices to
# restrict the field coordinates to occupied bins
sim.fields.from.occupied = function(fields.list, nplace.fields=100) {
  result = data.frame()
  for (animal in names(fields.list)) {
    for (date_str in names(fields.list[[animal]])) {
      is_occupancied_mat = !is.na(fields.list[[animal]][[date_str]][[1]])
      random.coors = sim.xy.coord.occupancy_mat(is_occupancied_mat, nplace.fields)
      
      ncoords = length(random.coors$field.x)
      result = bind_rows(result, 
                         data.frame(animal=rep(animal, ncoords),
                                    date=rep(date_str, ncoords),
                                    field.x=random.coors$field.x,
                                    field.y=random.coors$field.y,
                                    stringsAsFactors = FALSE))
    }
  }
  return(result)
}

calc.rew.distances.df = function(loc.df, field.x, field.y) {
  location.sets.df = dplyr::select(loc.df, -date) %>%
    dplyr::distinct()
  
  map_dfr(1:nrow(location.sets.df), ~ data.frame(
    animal=location.sets.df$animal[.x],
    location_set=location.sets.df$location_set[.x],
    min.rew.dist=calc.min.rew.dist(loc.df[.x, ], field.x, field.y)$rew.dist,
    stringsAsFactors=FALSE)
  )
}

calc.rew.distances.for.fields = function(loc.df, field.df) {
  field.df = data.table(field.df)
  
  map_dfr(1:nrow(loc.df), ~ {
    animal_x = loc.df$animal[.x]
    location_set_x = loc.df$location_set[.x]
    date_x = loc.df$date[.x]
    date.fields = field.df[animal == animal_x & date == date_x, ]
    if (nrow(date.fields) == 0) {
      return(data.frame())
    }
    data.frame(
      animal=animal_x,
      location_set=location_set_x,
      date=date_x,
      min.rew.dist=calc.min.rew.dist(loc.df[.x, ], date.fields$field.x, date.fields$field.y)$rew.dist,
      stringsAsFactors=FALSE)
  })
}
