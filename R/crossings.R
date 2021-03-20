source('locations.R')

get_crossings_n = function(zone_dist, zone_radius=rew_zone_radius, frame_rate=24) {
  is_inside = zone_dist < zone_radius
  min_frames_inside = 0.25  * frame_rate
  inside_diff = diff(is_inside) 
  entries = which(inside_diff == 1)
  exits = which(inside_diff == -1)
  if (length(entries) > length(exits)) {
    exits = c(exits, length(is_inside))
  }
  crossing_durations = exits - entries
  crossings_n = sum(crossing_durations > min_frames_inside)
  return(crossings_n)
}

crossings4positions = function(pos_x, pos_y, zone_x, zone_y, zone_radius=rew_zone_radius) {
  zone_dist = sqrt((pos_x - zone_x)^2 + (pos_y - zone_y)^2)
  return(get_crossings_n(zone_dist, zone_radius=zone_radius))
}
