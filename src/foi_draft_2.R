# Revised FOI script 65
# Developed in conjunction with Daniel Erro, Dillon Murphy, 
# Liam Quach, Rachel Roggenkemper from statistics consulting department
# 4/29/2025

require(terra)
require(sf)
require(svMisc)

data <- st_read("blitzdata.shp")
comphost <- subset(data, 
                   (Origin %in% c("Umbellularia californica",
                                  "Umbellularia  californica",
                                  "Notholithocarpus densiflorus")) &
                     (format(as.Date(Date), "%Y" ) <= "2023") &
                     SymbolID >= 4
)

# Reproject from lon/lat (EPSG:4326) to California Albers (EPSG:3310), 
# units in meters
comphost_proj <- st_transform(comphost, 3310)
pts_vect      <- vect(comphost_proj)    # SpatVector in meters
#   → pts_vect is a SpatVector of point geometries.  
#      You can index it like a matrix: pts_vect[i,] is the i‑th point.

# Create empty raster from that vector’s true extent
r_template <- rast(pts_vect, res = 100)  # auto-uses pts_vect’s bbox & CRS
values(r_template) <- 0
#   → A multi‑cell raster in EPSG:3310, 1 km resolution (set for testing), all values initialized to 0

# Accumulate FOI
alpha    <- 25  # for testing we’re using 5 km; later can change to 25 m
foi_rast <- r_template
for (i in seq_len(nrow(pts_vect))) { # loop over each point index
  # a) Extract the i-th point as a single‑feature SpatVector:
  this_pt <- pts_vect[i, ]
  
  # b) Compute a distance‐raster: for every cell in foi_rast, how far is its center
  #    from this single point? Resulting d_rast has same dimensions as foi_rast,
  #    with values in meters.
  d_rast <- distance(foi_rast, this_pt)
  
  # c) Apply the negative‐exponential kernel to those distances:
  #    exp(-d/alpha) yields a weight between 0 (far away) and 1 (zero distance).
  kernel_rast <- exp(-d_rast / alpha)
  
  # d) Add this point’s kernel to the running total FOI raster:
  #    cell‐wise addition accumulates contributions from each infection point.
  foi_rast <- foi_rast + kernel_rast
  
  # e) Optional progress bar to monitor status of for loop
  progress <- progress(i, 7764)
    Sys.sleep(0.02)
    if (i == 7764) message("Done!") 
}

    # At the end of the loop, foi_rast[cell] = sum_i exp(-dist(cell, p_i)/alpha).

# Mask out zeros
foi_mask <- foi_rast
values(foi_mask)[values(foi_mask) == 0] <- NA

# Write to tif file and change name
writeRaster(foi_rast, "foi_a25_2023.tif", overwrite=TRUE) 
