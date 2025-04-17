############# - Convert lidar footprint to polygons
############# - This can save time for selecting tile_id, and reducing number of
############# - files to be processed
#------------------------------------------------------------------------------#

#### - 1. Set global option -####
#------------------------------------------------------------------------------#
options(continue = "...", digits = 5, width = 80, timeout = 30, prompt = "> ")

#------------------------------------------- remove everything from environment
rm(list = ls(all.names = TRUE))
gc()

#### - 2. Load library - ####
#------------------------------------------------------------------------------#
library(lidR)
library(sf)

#### - 3. Read data - ####
#------------------------------------------------------------------------------#
st_abbr <- c("TX_2011", "TX_2018")

state = st_abbr[2]


# create directory if it does not exist
# if(!dir.exists(paste0("../00_lidTilesBB/", state))){
#   dir.create(path = paste0("../00_lidTilesBB/", state))
# }

lidar_path <- "./TX_2018/"


tile_list <- list.files("./TX_2018/",
                        pattern = ".laz$",
                        full.names = FALSE
)
# lid_state = read.delim(file = "./0_file_download_link_2018.txt",sep = "\t",header = FALSE)
# lid_state$sub_proj<- substr(x = lid_state$V1,start = 1,stop = nchar(lid_state$V1)-46)
# lid_state$laz_file<- substr(x = lid_state$V1,start = nchar(lid_state$V1)-45,stop = nchar(lid_state$V1)-0)
# lid_state$V1<- NULL
### - before writing function try 1)
# extract bounding box
# bounding box to polygon
# create tile_id
# create state
# create gpkg
# export layer to gpkg

## ----------- end of this function
# F drive is labelled as FPA_lidar1
tile_path <- "./TX_2018/"
tile_list <- list.files("./TX_2018/",pattern = "*.laz$",full.names = FALSE)
gpkg_path <- "./TX_2018_TileEXt/"
#lidar_path<- lid_state$sub_proj
#------------------------------------------------------------------------------#
# Load the required packages
library(doParallel)
library(parallel)


####------------ function 1 for rbind list ------------------------------------#
fastrbindsf = function(x, check = FALSE) {
  # https://gist.github.com/kadyb
  # https://github.com/r-spatial/sf/issues/798
  if (length(x) == 0) stop("Empty list")
  if (isTRUE(check)) {
    ref_crs = sf::st_crs(x[[1]])
    ref_colnames = colnames(x[[1]])
    for (i in seq_len(length(x))) {
      if (isFALSE(sf::st_crs(x[[i]]) == ref_crs)) stop("Diffrent CRS")
      if (isFALSE(all(colnames(x[[i]]) == ref_colnames))) stop("Diffrent columns")
    }
  }
  geom_name = attr(x[[1]], "sf_column")
  x = collapse::unlist2d(x, idcols = FALSE, recursive = FALSE)
  x[[geom_name]] = sf::st_sfc(x[[geom_name]], recompute_bbox = TRUE)
  x = sf::st_as_sf(x)
}

##----------------- Function 2 for getting extent -----------------------------#

# Modify the get_tile_extent function to use parallel processing
get_tile_extent_parallel <- function(tile_list, state, lidar_path, gpkg_path, num_cores = 4) {
  
  # Create a cluster using all all cores except two
  if (num_cores >= parallel::detectCores()) {
    warning("specified core number is equal to or more than available core\n
              defaulting to n-1 cores")
    num_cores <- parallel::detectCores() - 1
  } else {
    num_cores <- num_cores
    print(paste0("total of ", " ", num_cores, " are used, out of ", parallel::detectCores(), " ", "cores available"))
  }
  cl <- makeCluster(num_cores)
  
  # Register the parallel backend
  registerDoParallel(cl)
  

  # Parallel processing of each tile
 future_gpkg_files<-  foreach(tile = tile_list, .packages = c("sf", "lidR", "collapse")) %dopar% {
    ind <- which(tile_list == tile)
    print(paste0("processing file: ", ind, " remaining files are:", " ", length(tile_list) - ind))
    
    # Read the header of the lidar file to get the bounding box coordinates
    lid_head <- lidR::readLASheader(paste0(lidar_path, tile))
    
    # Extract the bounding box coordinates
    x1 <- lid_head@PHB$`Min X`
    x2 <- lid_head@PHB$`Max X`
    y1 <- lid_head@PHB$`Min Y`
    y2 <- lid_head@PHB$`Max Y`
    
    # Define a function to create a square polygon from the bounding box coordinates
    square <- function(x1, y1, x2, y2) {
      sf::st_polygon(list(matrix(c(x1, x2, x2, x1, x1, y1, y1, y2, y2, y1), 5)))
    }
    
    # Create a polygon representing the tile's bounding box
    tile_box <- square(x1, y1, x2, y2)
    
    # Get the coordinate reference system (CRS) EPSG code from the lidar header
    #crs_epsg <- sf::st_crs(lid_head)$input
    
    # Create a sfc of the tile's bounding box with the given CRS
    geometry <- sf::st_sfc(tile_box,crs = st_crs(lid_head))
    
    # Convert the sfc to a sf object
    bbox_poly <- sf::st_sf(geometry)
    
    # Add two columns to the spatial dataframe
    bbox_poly$state <- state
    bbox_poly$`tile_id` <- tile
    bbox_poly$`file_year`<- lid_head$`File Creation Year`
    bbox_poly$`file_month`<- round(lid_head$`File Creation Day of Year`/30.41,0)
   # bbox_poly$`native_prj`<- lid_head@VLR[["GeoAsciiParamsTag"]][["tags"]]
    
    # Transform the spatial dataframe to a different CRS (EPSG:5070)
    # This coordinates is USA_Albers_EqualArea_conformal Conic. 
    # This is going to be helpful, as we can use same coordinate system for our
    # study area.
    bbox_poly_trans <- sf::st_transform(bbox_poly, crs = "EPSG:6342")
    # Create the layer name
    layer_name <- paste0("tile_", ind)
    
    # Write the tile extent to a separate GeoPackage file for each core
    # st_write(bbox_poly_trans,
    #          dsn = paste0(gpkg_path, "BB_", state, "_", ind, ".gpkg"),
    #          layer = layer_name, driver = 'GPKG')
    return(bbox_poly_trans)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  registerDoSEQ()
  
  # file_list <- list.files(path = gpkg_path, pattern = paste0("BB_", state, "_[0-9]+.gpkg"), full.names = TRUE)
  final_gpkg <- paste0(gpkg_path, "BB_", state, ".gpkg")
  #sf::st_layers(future_gpkg_files[[1]])  # To ensure consistency among layers in all files
  sf::st_write(fastrbindsf(future_gpkg_files), dsn = final_gpkg, driver = "GPKG")
}



# Call the modified function with provided arguments
get_tile_extent_parallel(tile_list = tile_list,
                         state = state,
                         lidar_path = lidar_path,
                         gpkg_path = gpkg_path,
                         num_cores = 3)
