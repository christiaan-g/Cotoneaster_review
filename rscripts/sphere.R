sphere <- function(lat, lon, R = 6371000) {
  # Convert to radians
  phi <- lat * pi / 180
  lambda <- lon * pi / 180
  
  # Calculate coordinates
  x <- R * cos(phi) * cos(lambda)
  y <- R * cos(phi) * sin(lambda)
  z <- R * sin(phi)
  
  return(list(x, y, z))
}
