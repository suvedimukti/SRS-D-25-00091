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
library(imager)     # cimg based image analysis
library(lidR)       # lpc data 
library(RMCC)       # require for mcc algorithm
library(terra)      # shade/raster data
library(dplyr)

####- 4. Read Data -####
#------------------------------------------------------------------------------#
drive = "K:/"
tif_path  <- paste0(drive, "Mukti_TTU_CC/LBB_Project/downloadLidar/LB_2011_Out/")
out_path  <- paste0(drive, "Mukti_TTU_CC/LBB_Project/downloadLidar/")
rast_vect <- list.files(tif_path,pattern = "*.tif$",full.names = FALSE)


# fun1. validate inputs arguments
#------------------------------------------------------------------------------#
validate_input <- function(input_raster, min_height, drop_zero, min_win_neib) {
  # convert raster (from raster to terra supported raster)
  if ("RasterLayer" %in% class(input_raster)) {
    input_raster <- rast(input_raster)
  }
  
  # check_neigbhorhood type
  # queen or rook contiguity
  valid_neighborhoods <- c("queen", "rook")
  if (!tolower(min_win_neib) %in% valid_neighborhoods) {
    stop("Invalid input for 'min_win_neib'. Set to 'queen' or 'rook'.")
  }
  
  # check values are positive
  if (!is.null(min_height) && min_height <= 0) {
    stop("Minimum height must be set to a positive value.")
  }
  
  # # this is just for insurence, as min-height should get this
  #-----------------------------------------------------------------------------#
  if (drop_zero) {
    input_raster[input_raster < 0] <- 0
  }
  
  return(input_raster)
}

# predefine window for speed
# this vwf takes slighly longer time to run
# in the future, may be parallel options should be opted
# fun2. generate window mask
#--------------------------------------------------------------------------------#
gen_win_mask <- function(rast_lyr, window_radii, min_win_neib, win_diam) {
  window_masks <- lapply(window_radii, function(radius) {
    if (is.na(radius)) return(NULL)
    
    cell_size <- terra::res(rast_lyr)[1]
    # for small radius uses upto three cells, use a 3*3 matrix
    if (radius <= cell_size * 3) {
      if (tolower(min_win_neib) == "queen") {
        window_matrix <- terra::focalMat(rast_lyr, 3, type = "circle")
      } else if (tolower(min_win_neib) == "rook") {
        window_matrix <- terra::focalMat(rast_lyr, 3, type = "rectangle", edges = TRUE)
      }
    } else {
      # for larger-radius 
      radius_cells  <- round(radius / cell_size)
      window_matrix <- terra::focalMat(rast_lyr, radius_cells, type = "circle")
    }
    
    # enable padding if required
    pad <- (win_diam - ncol(window_matrix)) / 2
    pad <- ifelse(pad < 0, 0, pad)
    padded_window <- terra::extend(terra::rast(window_matrix), pad, fill = 0)
    
    # logical vector for padded window
    as.vector(padded_window[] != 0)
  })
  
  names(window_masks) <- window_radii
  return(window_masks)
}

#--------------- Check local maxima
# with window values and windows
# fun3. check local maxima
#------------------------------------------------------------------------------#
local_maxima <- function(window_values, windows) { # fun starts
  # calculate the center index for the current focal window 
  center_index  <- (length(window_values) + 1) / 2
  central_value <- window_values[center_index]
  
  # Return NA if the central value is NA
  if (is.na(central_value)) return(NA)
  
  # select the largest window mask from the pre defined list
  largest_mask <- windows[[length(windows)]]
  mask_indices <- which(largest_mask)
  
  # Compare the central value with all values in the neighborhood defined by the mask
  if (max(window_values[mask_indices], na.rm = TRUE) == central_value)
    return(central_value)
  else
    return(NA)
} # fun end

# function with variabel window filter
# process_local_maxima
# output polygon (sf), or terra(raster), or both as list
# fun 4. Process local maxima using variable window
#------------------------------------------------------------------------------#
process_local_maxima <- function(input_raster, 
                                 search_radius = 3, 
                                 min_height = 2.0,
                                 warnings = TRUE, 
                                 min_win_neib = "queen", 
                                 drop_zero = TRUE,
                                 id_field = "tree_id", 
                                 output_format = c("sf", "raster", "both")) {
  # check output format
  output_format <- match.arg(output_format)
  
  # vaalidate and convert the input raster
  # see funtion above validate_input
  rast_lyr <- validate_input(input_raster, min_height, drop_zero, min_win_neib)
  
  # check if raster cells are square round at 5 decimal points
  # this is hard coded
  pix_res <- round(terra::res(rast_lyr), 5)
  if (pix_res[1] != pix_res[2]) {
    stop("Input raster does not have square cells.")
  }
  if (pix_res[1] == 0) {
    stop("The map units of the input raster are too small (resolution equals 0).")
  }
  
  # Compute the raster range to ensure the absence of non-finite values
  raster_range <- terra::minmax(rast_lyr, compute = TRUE)[, 1, drop = TRUE]
  if (any(!is.finite(raster_range))) {
    stop("Input raster contains non-finite values.")
  }
  
  # filter out min-height cells, set cells with values below min_height to zero
  if (!is.null(min_height)) {
    if (min_height >= raster_range["max"]) {
      if (warnings) warning("'min_height' is set higher than the maximum cell value in the input raster")
      min_height <- 0
    }
    rast_lyr[rast_lyr < min_height] <- 0
  }
  
  # Determine the cell resolution and compute the sequence of window radii
  cell_size    <- pix_res[1]
  seq_floor    <- cell_size
  seq_ceiling  <- cell_size * ceiling(search_radius / cell_size)
  window_radii <- seq(seq_floor, seq_ceiling, by = cell_size)
  
  # window diameter in cell size
  win_diam <- ceiling((max(window_radii) / cell_size) * 2)
  if (win_diam %% 2 == 0) win_diam <- win_diam + 1
  
  # generate window for all radii
    window_masks <- gen_win_mask(rast_lyr, window_radii, min_win_neib, win_diam)
  
  # apply the focal operation on the raster using a fixed window (win_diam x win_diam matrix)
  # The custom function uses our helper to check for a local maximum in the neighborhood
  lmr <- terra::focal(
    rast_lyr,
    matrix(1, nrow = win_diam, ncol = win_diam),
    fun = function(values) local_maxima(values, windows = window_masks)
  )
  
  # Convert the resulting raster of local maxima into spatial points (sf object)
  lmp <- sf::st_as_sf(as.points(lmr, na.rm = TRUE, values = TRUE))
  
  # Rename the primary attribute to "height"
  names(lmp)[1] <- "height"
  
  # Filter out non-positive values and remove duplicate entries
  lmp <- lmp %>%
    dplyr::filter(height > 0) %>%
    dplyr::distinct(dplyr::across(-geometry), .keep_all = TRUE)
  
  # Add an ID field with unique identifiers
  lmp[[tolower(id_field)]] <- seq_len(nrow(lmp))
  lmp <- lmp[, c(tolower(id_field), "height", "geometry")]
  
  # Return output in the format specified by the user
  if (output_format == "sf") {
    return(lmp)
  } else if (output_format == "raster") {
    return(lmr)
  } else if (output_format == "both") {
    return(list(points = lmp, raster = lmr))
  }
}

# fun 5. watershed segmentation
#### -4.2 watershed segmentation -####
marker_watershed <- function(tree_tops, CHM, min_height = 2, format = "raster", 
                             ID_field = "tree_id",slope_neighbors = 8,slope_unit = "degrees",
                             min_slope = 0) {
  # Convert tree_tops to an sf object if needed
  if (inherits(tree_tops, "SpatialPointsDataFrame"))
    tree_tops <- sf::st_as_sf(tree_tops)
  
  # Validate tree_tops: must be an sf object with POINT geometry
  if (!inherits(tree_tops, "sf") || sf::st_geometry_type(tree_tops, 
                                                        by_geometry = FALSE) != "POINT")
    stop("Invalid input: 'tree_tops' should be an 'sf' object with 'POINT' geometry.")
  
  # Convert CHM to a SpatRaster if necessary
  if (inherits(CHM, "RasterLayer"))
    CHM <- terra::rast(CHM)
  if (!inherits(CHM, "SpatRaster"))
    stop("Invalid input: 'CHM' should be a 'SpatRaster' object.")
  
  # Validate 'format'
  if (!toupper(format) %in% c("RASTER", "POLYGONS", "POLYGON", "POLY","sf"))
    stop("Invalid input: 'format' must be set to either 'raster' or 'polygons'.")
  
  # Get CHM maximum and check usability
  CHM_max <- terra::global(CHM, fun = max, na.rm = TRUE)[1, "max"]
  if (!is.finite(CHM_max))
    stop("'CHM' does not contain any usable values.")
  if (min_height > CHM_max)
    stop("'min_height' is set higher than the highest cell value in 'CHM'.")
  
  # Extract treetop heights from CHM and filter tree_tops
  treetop_heights <- terra::extract(CHM, tree_tops)[, 2]
  tree_tops <- tree_tops[!is.na(treetop_heights) & (treetop_heights >= min_height), ]
  if (nrow(tree_tops) == 0)
    stop("No usable tree_tops. They are either outside of CHM's extent or below the 'min_height' threshold.")
  
  # Validate ID field in tree_tops
  if (!ID_field %in% names(tree_tops))
    stop("Could not find ID field '", ID_field, "' in the 'tree_tops' object.")
  tree_tops[[ID_field]] <- as.integer(tree_tops[[ID_field]])
  if (any(tree_tops[[ID_field]] == 0))
    stop("ID field cannot be equal to 0.")
  if (any(is.na(tree_tops[[ID_field]])))
    stop("ID field cannot contain NA values.")
  if (any(duplicated(tree_tops[[ID_field]])))
    warning("ID field contains duplicated values. Please check your tree_tops data.")
  
  # Create a mask: cells in CHM below min_height (or NA) are masked, then set to 0
  CHM_mask <- CHM < min_height
  CHM_mask[is.na(CHM)] <- TRUE
  CHM[CHM_mask] <- 0
  
  # Rasterize tree_tops using the ID field; background set to 0
  tree_tops_ras <- terra::rasterize(tree_tops, CHM, field = ID_field, background = 0)
  
  # Convert CHM and tree_tops raster to matrices (wide = TRUE preserves orientation)
  CHM_img   <- imager::as.cimg(terra::as.matrix(CHM, wide = TRUE))
  ttops_img <- imager::as.cimg(terra::as.matrix(tree_tops_ras, wide = TRUE))
  
  # Perform watershed segmentation using tree_tops as markers
  ws_img <- imager::watershed(ttops_img, CHM_img)
  
  # Convert watershed result back to a SpatRaster using CHM's extent and CRS
  ws_ras <- terra::rast(as.matrix(ws_img), extent = terra::ext(CHM), crs = terra::crs(CHM))
  
  # Reapply CHM mask to set areas below min_height to NA
  ws_ras[CHM_mask] <- NA
  
  if (toupper(format) %in% c("POLYGONS", "POLYGON", "POLY", "sf")) {
    # Polygonize the watershed raster
    polys <- sf::st_as_sf(terra::as.polygons(ws_ras, dissolve = TRUE))
    polys <- sf::st_cast(polys, "MULTIPOLYGON")
    polys <- sf::st_cast(polys, "POLYGON", warn = FALSE)
    if (nrow(polys) == 0)
      stop("No segments created.")
    names(polys)[1] <- ID_field
    row.names(polys) <- 1:nrow(polys)
    
    # Retain only polygons that intersect with tree_tops
    polys_out <- polys[lengths(sf::st_intersects(polys, tree_tops)) > 0, ]
    polys_out <- polys_out[!is.na(polys_out[[ID_field]]), ]
    poly_vfil <- sf::st_as_sf(fillHoles(terra::vect(polys_out)))
    slope_rast <- terra::terrain(CHM,v = "slope",
                                       neighbors = slope_neighbors,
                                       unit = slope_unit)
        # Extract the mean slope within each polygon
        mean_slope  <- terra::extract(slope_rast, poly_vfil, fun = mean, na.rm = TRUE)
        mean_height <- terra::extract(CHM, poly_vfil, fun = mean, na.rm = TRUE)
        min_ht      <- terra::extract(CHM, poly_vfil, fun = min, na.rm = TRUE)
        names(mean_slope)[2]    <- "mean_slope"
        names(mean_height)[2]   <- "mean_ht"
        names(min_ht)[2]        <- "min_ht"
        attbr = mean_slope %>% dplyr::left_join(mean_height, by = c("ID" = "ID")) %>%
          dplyr::left_join(min_ht, by = c("ID" = "ID")) #%>%
          # dplyr::left_join(mean_ndvi, by = c("ID" = "ID")) %>%
          # dplyr::left_join(min_ndvi, by = c("ID" = "ID"))

        # Add a "slope" column to the polygons
        #polys$slope <- mean_slopes
        polys_final <- cbind(poly_vfil, attbr %>% dplyr::select(-1, everything()))
        # Filter polygons based on the mean slope threshold
        if (!is.null(min_slope)) {
          polys_final <- polys_final[polys_final$slope >= min_slope, ]
          if (nrow(polys_final) == 0) {
            stop("No polygons meet the slope threshold")
          }
        }
    if (any(duplicated(polys_final[[ID_field]])))
      stop("Multiple segments with the same ID field. Check for duplicated tree_tops.")
    
    return(polys_final)
  } else {
    return(ws_ras)
  }
}


out_folder = "LB_2011_shp_lid_imp"

if(!dir.exists(paste0(out_path, out_folder))){
  dir.create(paste0(out_path,   out_folder))
}
#-----------------------------------------------------------------------------#
 # fun6: process all data based on 
 # fun 4, and fun5.
 # local maxima and watershed segmantation
 #-----------------------------------------------------------------------------#
##### - 5. Create Function - ####
process_all_data <- function(rast_vect, min_height = 2,
                             drop_zero = TRUE, out_folder)
{
  for(file in rast_vect){
    ind <- which(rast_vect == file)
    cat(sprintf("processing file: %s (%d/%d)\n",
                file, ind, length(rast_vect)))
    
    #-------------- read raster (file)
    img <- terra::rast(paste0(tif_path, file))
    
    # check raster condition, if raster's max < min_height, skip to next file
    # In agriculture field this is possible
    
    raw_max    <- terra::minmax(img)["max", ]
    if(raw_max < min_height){
      cat(sprintf("skipping: raw max = %.2f < min_height = %.2f\n",
                  raw_max, min_height))
      next
    }
    
    #-- filter out everything below 2m 
    mask         <- img > min_height
    filtered_chm <- img * terra::values(mask, ifelse = NA)
    
    # #--- detect tree_tops on the **filtered** raster
    ttops <- process_local_maxima(input_raster    = filtered_chm,
                                  search_radius   = 5,
                                  min_height      = min_height,
                                  min_win_neib    = "queen",
                                  output_format   = "sf")
    if(nrow(ttops) == 0){
      cat("no tree_tops found above threshold processing will happen to next file \n")
      next
    }
    cat(sprintf("found %d tree_tops\n", nrow(ttops)))
    
    # perform marker watershed with tryCatch
    #---------------------------------------
    poly <- tryCatch(
      {
        marker_watershed(tree_tops  = ttops,
                         CHM       = filtered_chm,
                         min_height = min_height,
                         format    = "POLY",
                         ID_field   = "tree_id",
						 slope_neighbors = 8,
						 slope_unit = "degrees",
						 min_slope = NULL)
      },
      error = function(e){
        cat(" mc_watershed failed:", e$message, "\n")
        return(NULL)
      }
    )
    if(is.null(poly) || nrow(poly) == 0){
      cat("no polygons to write moving to next one \n")
      next
    }
    
    # 6) write out
    name <- sub("\\.[^.]+$", "", file)
    out_file <- file.path(out_path, out_folder, paste0(name, "_", ind, ".shp"))
    cat("....writing:", out_file, "\n")
    # geopackge would have been awesome; however, as other processing were to be done in 
    # ArcGIS shapefile was used
    sf::st_write(poly, out_file, driver = "ESRI Shapefile", quiet = TRUE)
  }
}

# now process all data iteratively
#-----------------------------------------------------------#
process_all_data(rast_vect = rast_vect,
                 min_height = 2,
                 out_folder = out_folder,
                 drop_zero = TRUE)
