############# --------------- TREE CROWN based on lidar
################## ---------- Lubbock County, Texas, USA

####- 1. Global Option -####
# -----------------------------------------------------------------------------#
options(timeout   = 15,     prompt    = "< ",
        continue  = "...",  digits    = 5,
        width     =  80,    papersize = "letter" )

####- 2. Clear workspace -####
# -----------------------------------------------------------------------------#
rm(list = ls())

####- 3. Load Package -####
# -----------------------------------------------------------------------------#
library(sf)         # sf
library(lidR)       # lpc data 
library(RMCC)       # require for mcc algorithm
library(terra)      # shade/raster data

####- 4. Read Data -####
# -----------------------------------------------------------------------------#
drive         <- "F:/"
base_path     <- paste0(drive, "Mukti_TTU_CC/LBB_Project/downloadLidar/")
folder        <- paste0(base_path, "TX_2011/") # TX_2019 
lidar_list    <- list.files(folder,pattern = "*.laz$", full.names = FALSE)
out_folder    <- paste0(base_path, "LB_2011_Out/") # "LB_2019_Out/"
# -----------------------------------------------------------------------------#
if(!dir.exists(out_folder)){
  dir.create(out_folder)
}
# set lidr threads I have 24 threads
set_lidr_threads(threads = 20)

process_lidar <- function(lidar_list = lidar_list, base_path = base_path,
                          res = 0.6,         # 0.6 m to match 2019 NAIP image
                          max_ht = 25.0,     # define maximum height to remove some tall buildings 
                          out_folder = out_folder
                         ) { # start function
  
  for (file in lidar_list) {
    ind <- which(lidar_list == file)
    print(paste0(
      "processing file: ", ind, ": remaining files are :", " ",
      length(lidar_list) - ind
    ))
    #----------------------------------------------------------------------------#
    dat   <- lidR::readLAS(files = paste0(folder, file), filter = "-drop_withheld")
    ldat  <- st_transform(dat,crs = "epsg:26914")
    ldat@data$Z <- ldat@data$Z*0.304800609601219           # convert height from feet to meter
    ldat  <- las_update(las = ldat)
    ldat@data <- ldat@data[(Classification%in%c(1,2,3,4,5,6,8,9))]
    #-------------------------------------------------------------------------#
    
    las2 <- normalize_height(ldat, res = res, algorithm = tin())
    ## read shapefile
    # While we tried using a building shapefile to mask out building  points; due to spatial differences and building layer was not 
    # up to date, we decided to remove them later using slope and ndvi using NAIP image
    #-------------------------------------------------------------------------------#
    # 
    #     building <- st_buffer(building, dist = 2, endCapStyle = "ROUND")
    #     building <- building %>% st_transform(., crs = st_crs(ldat))
    #     # 
    #     las3 <- merge_spatial(las2, building, "inpoly")
    #     # 
    #     las3$Classification[las3$inpoly == TRUE] <- 6L
    #     las4 <- las3[las3$Classification != c(6)]
    
    # chm <- rasterize_canopy(
    #   las = las2,
    #   res = res, algorithm = dsmtin(max_edge = 20)
    # )
    # chm = rasterize_canopy(las = las2, res = res,
    #                        algorithm = pitfree(c(0,5,10,15,20,25,30,40,50,60,70,80,90,100), c(0, 2))
    # )
    
    chm = rasterize_canopy(las = las2, res = res,
                           algorithm = dsmtin())
    
    chm[chm > max_ht]<- 0.01
    w <- matrix(1, 3, 3)
    chm1 <- terra::focal(chm, w, fun = median, na.rm = TRUE)
    
    file_name <- substring(file, 1, last = nchar(file) - 4)
        print(paste0("writing chm ", file_name))
    plot(chm1)
    terra::writeRaster(chm1,filename =paste0(out_folder, file_name, ".tif"), overwrite = TRUE )
  }
}

# process lidar files (separate)
process_lidar(lidar_list = lidar_list,base_path = base_path,res = 0.5,max_ht = 100,out_folder = out_folder)
