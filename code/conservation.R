conservation <- function(initial_temps, intensity, effect_size) {
  
  # vector of all clutch temps, not organized by female
  raw_clutch_temps <- unlist(initial_temps)
  
  # generate vector of indexes of clutches with temperature changes
  clutches_affected <- sample(x = 1:length(raw_clutch_temps), 
                              size = round(length(raw_clutch_temps) * intensity), 
                              replace = FALSE)
  
  # create vector of temp changes for all clutches
  temp_changes <- rep(0, length(raw_clutch_temps))
  temp_changes[clutches_affected] <- effect_size
  
  # new clutch temperatures with temperature change
  new_clutch_temps <- raw_clutch_temps - temp_changes
  
  # and now put back into list form
  final_temps <- relist(flesh = new_clutch_temps, 
                        skeleton = initial_temps)
  
  return(final_temps)
  
}