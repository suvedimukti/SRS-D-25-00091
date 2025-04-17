# SRS-D-25-00091
Evaluation of trends in urban tree cover in the city of Lubbock, Texas, using LiDAR and NAIP imagery fusion

This repository has the following files:

1. 00_01_2011_lidar_data.txt
   This file contains the 2011 lidar data link 
2. 00_02_lidar_data_2019.txt
   This file contains the 2019 lidar data file link
3. 00_get_lidar_bounding_boxes.R
   Read the lidar data and get the bounding boxes such that only relevant lidar tiles are read and processed
4. 01_process_lidar_data.R
   The document contains a variable window filter function for local maxima detection
6. 02_process_chm_lcl_maxima_ws_seg.R
   This file contains a function for marker-controlled watershed segmentation.
   
