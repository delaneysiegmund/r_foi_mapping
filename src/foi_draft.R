

library(sf)         # For handling shapefiles
library(raster)     # For raster operations
library(sp)         # For spatial operations
library(geosphere)  # For calculating distances (alternative to rgeos)

#Load the shapefile containing disease incidence points (assumed as 'disease_points.shp')
data <- st_read("C:/Users/delan/OneDrive/Documents/Thesis/R/blitzdata/blitzdata.shp")

#Set parameters of infection data 
#Origin is viable host species, Date filters by year, and SymbolID refers to 
#infection status (which is on a scale from 0-6, 4+ is SOD positive)
###Change year as needed

comphost <- subset(data, (data$Origin == "Umbellularia californica" |
                     data$Origin == "Umbellularia  californica" | 
                     data$Origin == "Notholithocarpus densiflorus") 
                   & format(as.Date(data$Date),"%Y") == 2010
                   & data$SymbolID >= 4)
###

#Set the CRS 
st_crs(comphost) <- 4326  # WGS84 

#Convert the shapefile to a spatial object if it's not already
comphost_sp <- as(comphost, "Spatial")

# Define the parameters for the negative exponential kernel
# For the negative exponential kernel, we define the decay parameter (alpha)
alpha <- 25  #Consistent with literature

#Calculate the distance matrix (from each point to every other point)
#Here we use geosphere's distVincentySphere function for distance calculation in meters
dist_matrix <- distVincentySphere( 
  cbind(comphost_sp@coords[, 1], comphost_sp@coords[, 2])
)

#Calculate the dispersal probability for each point
#Using a negative exponential function
prob_matrix <- exp(-dist_matrix / alpha)

#Creating raster grid 
raster_extent <- extent(-124.0, -118.0, 34.0, 43.0)  # California and Southern Oregon
raster_res <- 100  # 100m resolution
study_area <- raster(raster_extent, res = raster_res)

### END OF WORK 4/10

###THING THAT DIDNT WORK #1:
#I asked Chat GPT for help writing this part, but it didn't run, even after
#checking arguments

for (row in 1:nrow(study_area)) {
  for (col in 1:ncol(study_area)) {
    
    # Get the center coordinates of the current raster cell
    cell_coords <- xyFromCell(comphost_sp)
    
    # Initialize the FOI value for the current cell
    FOI_value <- 0
    
    # Iterate over all points in the shapefile
    for (i in 1:length(comphost_sp)) {
      
      # Get the coordinates of the point
      point_coords <- coordinates(comphost_sp)[i, ]
      
      # Calculate the distance between the raster cell and the point
      distance <- spDistsN1(as.matrix(cell_coords), point_coords, longlat = TRUE)
      
      # Apply the negative exponential kernel to the distance
      FOI_value <- FOI_value + negative_exponential_kernel(distance, lambda = 25)
    }
    
    # Assign the FOI value to the raster cell
    FOI_raster[row, col] <- FOI_value
  }
}  

### THING THAT DIDNT WORK #2 
#Some arguments on this one might be really funky 

#Initialize an empty raster with zeros (or any default value)
values(study_area) <- 0  # Initialize all values to 0

# Loop through each point to calculate the dispersal and populate the raster
for (i in 1:nrow(study_area)) {
# Calculate distances from the ith point to all other points
  distances <- distVincentySphere(
    cbind(comphost_sp@coords[i, 1], comphost_sp@coords[i, 2]),
    cbind(comphost_sp@coords[, 1], comphost_sp@coords[, 2])
  )
    
# Calculate the probability values using the negative exponential kernel
prob_values <- exp(-distances / alpha)
  
#For each cell in the raster, sum the probabilities based on the distance from the point
  for (cell in 1:length(study_area)) {
    # Get the coordinates of the cell
    cell_coords <- xyFromCell(study_area, cell)
    
    # Calculate distance from the current point (i) to the raster cell
    distance_to_cell <- distVincentySphere(
      cbind(comphost_sp@coords[i, 1], comphost_sp@coords[i, 2]),
      cell_coords
    )
    
    # Calculate the probability for this cell
    prob_for_cell <- exp(-distance_to_cell / alpha)
    
    # Get the current values of the raster
    current_values <- values(study_area)
    
    # Add the calculated probability to the raster at the appropriate cell
    current_values[cell] <- current_values[cell] + prob_for_cell
    
    # Set the updated values back to the raster
    values(study_area) <- current_values
  }
}

#This is other processing code to 
#Normalize the raster (optional, to get probabilities between 0 and 1)
spread_raster <- FOI_raster / max(values(FOI_raster), na.rm = TRUE)

#Plot the resulting raster to visualize the probability of disease spread
plot(spread_raster, main = "FOI 2023")

writeRaster(spread_raster, "C:/Users/delan/OneDrive/Documents/Thesis/R/blitzdata/blitzdata.shp/FOI2023.tif", format = "GTiff")
