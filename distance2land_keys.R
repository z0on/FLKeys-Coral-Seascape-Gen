rm(list=ls())
library(RDAforest)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(viridis)
library(sf)
library(spatstat)
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")

mapdata= ne_countries(scale="large",continent="north america")

ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/rasters_keys_jan2025_gebcoDepth.RData")
ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Agaricia_env_ll_jan2025.RData")
env.a=env
latlon.a=latlon
ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Porites_env_ll_jan2025.RData")
env.p=env
latlon.p=latlon

library(maps)

# Load necessary packages
library(sf)
library(tmap)

# Download and read a shapefile of US states (or a coastline dataset)
florida <- st_read("~/Dropbox/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  dplyr::filter(STUSPS == "FL")  # Filter only Florida

# Get the boundaries of the Florida Keys using a bounding box
keys_bbox <- st_bbox(c(xmin = min(XY$lon), xmax = max(XY$lon), ymin = min(XY$lat), ymax = max(XY$lat)), crs = st_crs(florida))

# Crop Florida to only the Florida Keys area
florida_keys <- st_crop(florida, keys_bbox)

# # Plot the result
tmap_mode("view")
# tm_shape(florida_keys) + 
#   tm_polygons(col = "lightblue") + 
#   tm_borders()
# 
tm_shape(florida_keys)

plot(st_geometry(florida_keys))
points_input <- locator(type = "o")  # Collects x, y coordinates

# Convert selected points into a dataframe
coords <- data.frame(x = points_input$x, y = points_input$y)

# Ensure the polygon is closed (first and last points must be the same)
coords <- rbind(coords, coords[1,]) 

# dry tortugas
plot(st_geometry(florida_keys))
points_nput <- locator(type = "o")  # Collects x, y coordinates
coords.dt <- data.frame(x = points_input$x, y = points_input$y)
coords.dt=rbind(coords.dt, coords.dt[1,]) 


# Convert to two-polygon thingy
polygon_sf <- st_sf(
  id = c(1, 2),
  geometry = st_sfc(
    st_polygon(list(as.matrix(coords))),
    st_polygon(list(as.matrix(coords.dt)))
  ),
  crs = 4326
)

polygon_sf <- st_make_valid(polygon_sf)

# Plot selected polygon on top of the shapefile
plot(st_geometry(florida_keys))
plot(polygon_sf, col = rgb(1, 0, 0, 0.5), add = TRUE)  # Semi-transparent red polygon
save(polygon_sf,file="keys_polygon_sf.RData")
load("keys_polygon_sf.RData")

points_sf=st_as_sf(XY, coords = c("lon", "lat"), crs = st_crs(4326))
# Compute distances from each point to all polygons
dist_matrix <- st_distance(points_sf, polygon_sf)
# Find the minimum distance for each point (i.e., shortest distance to any polygon)
rasters.pre$dis2land=rasters.post$dis2land= apply(dist_matrix, 1, min)
plot(lonlat2raster(XY,data.frame(rasters.pre[,c("dis2land")])),bty="n")

points_sf=st_as_sf(latlon.a, coords = c("lon", "lat"), crs = st_crs(4326))
# Compute distances from each point to all polygons
dist_matrix <- st_distance(points_sf, polygon_sf)
env.a$dis2land=apply(dist_matrix, 1, min)

points_sf=st_as_sf(latlon.p, coords = c("lon", "lat"), crs = st_crs(4326))
# Compute distances from each point to all polygons
dist_matrix <- st_distance(points_sf, polygon_sf)
env.p$dis2land=apply(dist_matrix, 1, min)

names(rasters.pre)

save(XY,rasters.pre,rasters.post,file="~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/rasters_keys_feb2025.RData")
env=env.a
latlon=latlon.a
save(env,latlon,file="~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Agaricia_env_feb2025.RData")
env=env.p
latlon=latlon.p
save(env,latlon,file="~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Porites_env_feb2025.RData")
