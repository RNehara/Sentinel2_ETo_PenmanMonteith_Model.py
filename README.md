# PET Estimation Using Sentinel-2 and ASCE-PM Method    

This repository provides a Python script for estimating **daytime potential evapotranspiration (PET)** using the **Penman-Monteith ASCE** approach. It uses **Sentinel-2 L2A imagery** already with atmospheric and terrain correction and incorporates NDVI and surface albedo.    

**Features**   
- Calculates **potential evapotranspiration (PET)** using the ASCE-EWRI method (Allen et al., 1998; ASCE-EWRI, 2005)  
- Supports **NDVI** and **surface albedo** estimation from Sentinel-2 bands (Bonafoni & Sekertekin, 2020).
- Resamples 20m bands (B11, B12) to 10m resolution for unified processing  
- Exports raster outputs of PET, NDVI, and albedo  
- Allows for interactive input of meteorological and geographic data from a weather station

**Required python packages**  
import numpy as np  
import os  
import rasterio  
from rasterio.enums import Resampling  
import glob  
from rasterio.warp import reproject    

**Required Sentinel-2 Bands**  
- B02, B03, B04, B08 (10m)  
- B11, B12 (20m, automatically resampled)

**Edit the paths at the beginning of the script**  
dir10m = "path/to/Sentinel2/10m_bands/"  
dir20m = "path/to/Sentinel2/20m_bands/"  
output = "path/to/output/folder/"    

**Meteorological Input (manual entry)**  
Air temperature (°C) – reference temperature at the time of image acquisition  
Relative humidity (decimal) – e.g., 0.5 for 50%  
Atmospheric pressure (kPa) – e.g., 101.3 kPa for standard pressure  
Tmax / Tmin (°C) –  maximum and minimum air temperature at the time of image acquisition  
Wind speed (m/s) – measured at 2 meters above the ground  
Solar radiation (MJ/m²) – incoming shortwave radiation  
Day of year (DOY) – e.g., 15 for January 15  
Elevation (m) – altitude of the weather station above sea level  

**Outputs**  
- `PET.tif`: Estimated potential evapotranspiration (mm/hour) raster  
- `NDVI.tif`: NDVI raster  
- `albedo.tif`: Surface albedo raster    

**References**
- ASCE-EWRI. The ASCE standardized reference evapotranspiration equation (2005). In: ALLEN RG et al. (eds.), Report 0-7844-0805-X St. Louis: American Society of Civil Engineers, Environmental Water Resources Institute, p. 1- 69.
- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300(9), D05109.
- Bonafoni, S., & Sekertekin, A. (2020). Albedo retrieval from Sentinel-2 by new narrow-to-broadband conversion coefficients. IEEE Geoscience and remote sensing letters, 17(9), 1618-1622.

**Author**  
Rodrigo Nehara  
r.nehara@usp.br  
GitHub: @RNehara    

---    

*Note: It is recommended to use Sentinel-2 images and meteorological data from the same day and hour (or as close as possible) for consistency. Meteorological data should be from the same location and the hour immediately after the satellite acquisition time (e.g., if the image was acquired at 13:15 PM, use data from 14:00 PM).*
