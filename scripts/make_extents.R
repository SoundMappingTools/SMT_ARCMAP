# Script to create a square buffer around an input point

# Author: A.C. Keyel
# Copied from Sierrahlpr.R on 12/6/2016
# Last modified: 4/11/2017

# Needs to be run manually
#install.packages("maptools")
#install.packages("sp")
#install.packages("rgdal")

# First, a function to create the square buffer is defined
# Second, the inputs are read from a Python call and the function is run

#' Define function to create a square buffer
#' Create a square polygon around a point to use as a spatial extent
#' 
make_square_buffer = function(X, Y, buffer.size, buffer.file, utm.zone){

  library(maptools)
  library(sp)
  library(rgdal)

  # Convert to numeric format
  X = as.numeric(as.character(X))
  Y = as.numeric(as.character(Y))
  buffer.size = as.numeric(as.character(buffer.size))
  
  xmin = X - buffer.size
  xmax = X + buffer.size
  ymin = Y - buffer.size
  ymax = Y + buffer.size
  pt1 = c(xmin,ymin)
  pt2 = c(xmin, ymax)
  pt3 = c(xmax, ymax)
  pt4 = c(xmax, ymin)
  my.crs = sprintf("+proj=utm +zone=%s+ datum=WGS84", utm.zone)
  coords = matrix(c(pt1,pt2,pt3,pt4,pt1), ncol = 2, byrow = T)
  poly = Polygon(coords, hole = FALSE) #proj4string = my.crs,  
  polys = Polygons(list(poly), ID = 1)
  sp.poly = SpatialPolygons(list(polys), proj4string = CRS(my.crs))
  sp.poly.df = SpatialPolygonsDataFrame(sp.poly, data = data.frame(ID=1))
  # Write polygon to file
  writeOGR(sp.poly.df, buffer.file, layer=buffer.file, driver="ESRI Shapefile", overwrite_layer = TRUE)
  
}

# Run the function for a set of inputs
args = commandArgs(trailingOnly=TRUE)
X = args[1]
Y = args[2]
buffer.size = args[3]
buffer.file = args[4]
utm.zone = args[5]

make_square_buffer(X, Y, buffer.size, buffer.file, utm.zone)
  


# End of File