# Ensure counter-clockwise direction
is.ccw <- function(p) {
  tryCatch({owin(poly=p); T}, error=function(e) F)
}

# Check that all polygons were traversed in the right direction
ensure.ccw <- function(polygon.list) {
  lapply(polygon.list, function(p) {
    # Check the first polygon (outer boundary)
    if (!is.ccw(p)) {
      p$x <- rev(p$x)
      p$y <- rev(p$y)
    }
    p
  })
}

# Function to convert polygon to required list format
convert.to.list.format <- function(sf.object) {
  if (class(sf.object)[[1]] == "sf") {
    polygons.list <- list()
    for ( i in 1:nrow(sf.object)) {
      sfc <- st_geometry(sf.object[i,])
      if (class(sfc)[[1]] == "sfc_MULTIPOLYGON") {
        sfc <- st_cast(sfc, "POLYGON")
      }
      for (poly in sfc) {
        polygons.list[[length(polygons.list) + 1]] <- list(x = st_coordinates(poly)[,1], 
                                                           y = st_coordinates(poly)[,2])
      }
    }
  } else if (class(sf.object)[[1]] == "sfc_MULTIPOLYGON") {
    sfc <- st_geometry(sf.object) %>%
      st_cast("POLYGON")
    polygons.list <- lapply(sfc, function(poly) {
      list(x = st_coordinates(poly)[,1], 
           y = st_coordinates(poly)[,2])
    })
  } else if (class(sf.object)[[1]] == "sfc_POLYGON") {
    polygons.list <- lapply(sf.object, function(poly) {
      list(x = st_coordinates(poly)[,1], 
           y = st_coordinates(poly)[,2])
    })
  }
  
  # If the object has only one row, we unlist one level to adhere 
  # to the described format for `ppp` objects
  if (length(polygons.list) == 1) {
    polygons.list <- polygons.list[[1]]
  }
  
  # Ensure counter-clockwise
  polygons.list <- ensure.ccw(polygons.list)
  
  return(polygons.list)
}