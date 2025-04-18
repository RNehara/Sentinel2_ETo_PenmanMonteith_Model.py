# ======================================================================
#                                        Rodrigo Nehara
#                                        r.nehara@usp.br
#                                        github
#                                        11-21-2024 (mm/dd/yyyy)
# ASCE-PM Model for PET estimation 
# Also albedo and NDVI estimation
# Sentinel-2 L2A image
# running Python 3.x
# ======================================================================

# ------------------------------------------------------------
# libraries
import numpy as np  # For mathematical operations
import os  # For working with file paths
import rasterio  # For working with raster data
from rasterio.enums import Resampling # For resampling raster data
import glob  # For pattern matching with file names
from rasterio.warp import reproject # For raster transform

# To calculate PET
sigma = 4.903e-9/24  # Stefan-Boltzmann constant (W m^-2 K^-4)
lapse = 0.0065 # Lapse rate (temperature decrease with height) in °C/m
G = 0 # Soil heat flux: Typically small in daily scales and can be assumed zero if measurements are unavailable (ALLEN et al., 1998).
cnh = 0.37 # Constants related to grass during the day (hourly) (K * mm * s³ / Mg / h)
cdh = 0.24 # Constants related to grass during the day (hourly) (s / m)

# ------------------------------------------------------------
# edit here
# ------------------------------------------------------------

    # This function requires Sentinel-2 L2A bands 2, 3, 4, 8, 11, and 12.
    # Keep only the necessary bands in each directory.
    # dir: the directory (folder paths) where the images are stored.
    # output: the directory (folder paths) where the images should be exported
    # It is recommended not to rename the bands or the metadata files in the folder.

# Define paths analysis
dir10m = "C:/Users/nehar/OneDrive/Documentos/Sentinel_Article/Images/Image4/10m/"
dir20m = "C:/Users/nehar/OneDrive/Documentos/Sentinel_Article/Images/Image4/20m/"

output = "Downloads/"
     
# -------------- edition ends here ---------------------------


print("------------------------------------------------------------")
print('Uploading Sentinel-2 10m bands')
# Find all .jp2 files in the 10m directory
jp2_files_10m = glob.glob(os.path.join(dir10m, "*.jp2"))

# List of expected band names with leading zeros (matching real file names)
expected_10m_bands = ["B02", "B03", "B04", "B08"]

# Function to find the file path for a given band code
def find_band_path(file_list, band_code):
    for path in file_list:
        filename = os.path.basename(path)
        if f"_{band_code}." in filename or f"_{band_code}_" in filename:
            return path
    raise FileNotFoundError(f"Band {band_code} not found in provided files.")

# Build a dictionary mapping band names to their corresponding file paths
band_10m_paths = {
    band: find_band_path(jp2_files_10m, band) for band in expected_10m_bands
}

# Optional: print loaded band mapping
for band, path in band_10m_paths.items():
    print(f"{band} → {path}")
    
print('It is done')
print("------------------------------------------------------------")

print('Uploading Sentinel-2 20m bands')
# Find all .jp2 files in the 20m directory
jp2_files_20m = glob.glob(os.path.join(dir20m, "*.jp2"))

# List of expected band names at 20m resolution (with leading zeros)
expected_20m_bands = ["B11", "B12"]

# Reuse the same helper function (already handles B01 format)
def find_band_path(file_list, band_code):
    for path in file_list:
        filename = os.path.basename(path)
        if f"_{band_code}." in filename or f"_{band_code}_" in filename:
            return path
    raise FileNotFoundError(f"Band {band_code} not found in provided files.")

# Build dictionary mapping each band name to its actual file path
band_20m_paths = {
    band: find_band_path(jp2_files_20m, band) for band in expected_20m_bands
}

# Print result for verification
for band, path in band_20m_paths.items():
    print(f"{band} → {path}")

print('It is done')
print("------------------------------------------------------------")

print("Resampling Sentinel-2 20m bands to 10m resolution")

# Use band B4 (10m) as reference for resolution and dimensions
ref_band_path = band_10m_paths["B04"]
with rasterio.open(ref_band_path) as ref:
    ref_transform = ref.transform
    ref_crs = ref.crs
    ref_height = ref.height
    ref_width = ref.width

# Dictionary to store resampled 20m bands in memory
resampled_bands_20m = {}

# Loop through each 20m band and resample
for band_name, input_path in band_20m_paths.items():
    with rasterio.open(input_path) as src:
        data = src.read(1)  # Read original 20m data
        resampled = np.empty((ref_height, ref_width), dtype=src.dtypes[0])  # Allocate 10m-sized array

        # Perform resampling
        reproject(
            source=data,
            destination=resampled,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_transform,
            dst_crs=ref_crs,
            resampling=Resampling.bilinear
        )

        # Store result in dictionary
        resampled_bands_20m[band_name] = {
            "data": resampled,
            "meta": {
                "crs": ref_crs,
                "transform": ref_transform,
                "width": ref_width,
                "height": ref_height,
                "dtype": src.dtypes[0]
            }
        }
print('It is done')
print("------------------------------------------------------------")


print("Stacking Sentinel-2 bands")
# Define the band order for the stack
band_order = ["B02", "B03", "B04", "B08", "B11", "B12"]

# Create a list to hold each band as a 2D array
stack_list = []

for band_name in band_order:
    if band_name in band_10m_paths:
        # Read from 10m band directly
        with rasterio.open(band_10m_paths[band_name]) as src:
            band_array = src.read(1)
    elif band_name in resampled_bands_20m:
        # Use resampled 20m band already stored in memory
        band_array = resampled_bands_20m[band_name]["data"]
    else:
        raise ValueError(f"Band {band_name} not found in any source.")

    stack_list.append(band_array)

# Stack bands into a single 3D array: shape = (bands, height, width)
stack = np.stack(stack_list)
print("It is done")
print("------------------------------------------------------------")

print('Starting processing')
# Unpack each band from the resampled stack array
B02  = stack[0]  # Band 2
B03  = stack[1]  # Band 3
B04  = stack[2]  # Band 4
B08  = stack[3]  # Band 8
B11 = stack[4]  # Band 11
B12 = stack[5]  # Band 12

# Calculate NDVI with protection against divide-by-zero and invalid operations
with np.errstate(divide='ignore', invalid='ignore'):
    ndvi = (B08.astype("float32") - B04.astype("float32")) / (B08.astype("float32") + B04.astype("float32"))
    ndvi = np.where((ndvi < -1) | (ndvi > 1), np.nan, ndvi)
print('It is done')
print("------------------------------------------------------------")

print('Calculation Land Surface Albedo estimation')
  # This function calculates land surface albedo (Bonafoni & Sekertekin (2020).
  # image: pre-processed image with bands B2, B4, B8, B11

albedo = ((0.2266 * B02) + (0.1236 * B03) + (0.1573 * B04) + (0.3417 * B08) + (0.1170 * B11) + (0.0338 * B12)) / 10000
albedo = np.where((albedo < -1) | (albedo > 1), np.nan, albedo)

print('It is done')
print("------------------------------------------------------------")

print('Starting Potential Evapotranspiration estimation')

DOY = input(
    "Please type the day of the year from Sentinel-2 image aquisition (01/01/2024 (mm/dd/yyyy) e.g., 1: ")
print("------------------------------------------------------------")

welev = input(
    "Please type the elevation of the station above sea level (m): ")
print("------------------------------------------------------------")

t_air = float(input(
    "Please type the reference air temperature from weather station (°C): "))
print("------------------------------------------------------------")

rel_hum = float(input(
    "Please type the reference relative humidity from weather station (percentage, e.g., 0.5 for 50%): "))
print("------------------------------------------------------------")

p_atm = float(input(
    "Please type the reference atmospheric pressure from weather station (KPa, e.g., 10KPa for 100HPa): "))
print("------------------------------------------------------------")

Tmax = float(input(
    "Please type the maximum temperature value in the weather station (ºC): "))
print("------------------------------------------------------------")

Tmin = float(input(
    "Please type the minimum temperature value in the weather station (ºC) : "))
print("------------------------------------------------------------")

wind = float(input(
    "Please type the wind speed value in the weather station (height of the 2 m (m/s) : "))
print("------------------------------------------------------------")

rs = float(input(
    "Please type the radiation value in the weather station ((MJ/m², e.g., 1 MJ/m² for 1000 KJ/m²) : "))
print("------------------------------------------------------------")

Hour = float(input(
    "Please type the latitude in the weather station (HOUR, use a negative value for locations in the Southern Hemisphere): "))
print("------------------------------------------------------------")

Min = float(input(
    "Please type the latitude in the weather station (MIN) : "))
print("------------------------------------------------------------")

print('Calculating Psychrometric Coefficient (kPa/°C)')
# Calculate the psychrometric coefficient (kPa/°C)
y = 0.665 * 0.01 * p_atm
print('It is done')
print("------------------------------------------------------------")

print('Calculating slope of the vapor pressure curve with respect to temperature (kPa/°C)')
# Calculate the slope of the vapor pressure curve with respect to temperature (kPa/°C)
delta = (4098 * (0.6108 * np.exp((17.27 * t_air) / (t_air + 237.3))) / ((t_air + 237.3)**2))
print('It is done')
print("------------------------------------------------------------")

print('Calculating Saturated Vapor Pressure (kPa)')
# Calculate the Saturated Vapor Pressure (kPa)
es = 0.6107 * np.exp(17.269 * t_air / (t_air + 273.3))
print('It is done')
print("------------------------------------------------------------")

print('Calculating Vapor Pressure (kPa)')
# Calculate vapor pressure (ea) in kPa
ea = (rel_hum * es) / 100
print('It is done')
print("------------------------------------------------------------")

print('Calculating relative distance factor of Earth to the Sun')
# Calculate Relative distance factor of Earth to the Sun
Dr = 1 + 0.033 * np.cos(((2 * np.pi) / 365) * int(DOY))
print('It is done')
print("------------------------------------------------------------")

print('Calculating latitude in degrees')
# Calculate latitude in degrees
latitude_deg = Hour + Min / 60.0  # Use 60.0 to ensure floating-point division
print('It is done')
print("------------------------------------------------------------")

print('Calculating latitude from degrees to radians')
# Convert latitude from degrees to radians
phi = np.radians(latitude_deg)
print('It is done')
print("------------------------------------------------------------")

print('Calculating solar declination (δ) in radians')
# Calculate solar declination (δ) in radians
δ = 0.409 * np.sin((((2 * np.pi) / 365) * int(DOY)) - 1.39)

print('Calculating X')
# Calculate X
X = 1 - ((np.tan(phi)**2) * (np.tan(δ)**2))
print('It is done')
print("------------------------------------------------------------")

print('Calculating sunset hour angle in radians')
# Calculate sunset hour angle in radians)
omega_s = (np.pi / 2) - np.arctan((-np.tan(phi) * np.tan(δ)) / np.sqrt(X))
print('It is done')
print("------------------------------------------------------------")

print('Calculating extraterrestrial radiation')
# Calculate extraterrestrial radiation
Ra = (118.08 / np.pi) * Dr * (omega_s * np.sin(phi) * np.sin(δ) + np.cos(phi) * np.cos(δ) * np.sin(omega_s))
print('It is done')
print("------------------------------------------------------------")

print('Calculating net shortwave radiation')
# Calculate net shortwave radiation
Rns = (1 - albedo) * rs
print('It is done')
print("------------------------------------------------------------")

print('Calculating clear-sky solar radiation')
# Calculate clear-sky solar radiation
Rso = (0.75 + 2e-5 * int(welev)) * Ra
print('It is done')
print("------------------------------------------------------------")

print('Calculating net longwave radiation')
# Calculate net longwave radiation
Rnl = sigma * ((((Tmax + 273.15)**4 + (Tmin + 273.15)**4) / 2) * ((0.94 - 0.14 * np.sqrt(ea)) * (1.35 * (rs / Rso) - 0.35)))
print('It is done')
print("------------------------------------------------------------")

print('Calculating net radiation')
# Calculate net radiation
Rn = Rns - Rnl
print('It is done')
print("------------------------------------------------------------")

print('Calculating potential evapotranspiration using Penman-Monteith FAO-52 method with ASCE-EWRI (2005) reference values')
# Calculate Potential evapotranspiration (PET)
PET = ((0.408 * delta * (Rn - G)) + ((y * cnh * wind * (es - ea)) / (t_air + 273))) / (delta + y * (1 + cdh * wind))
print('It is done')
print("------------------------------------------------------------")

print("Exporting rasters")
# Dictionary of rasters in memory (filename → array)
rasters = {
    "albedott.tif": albedo,
    "PETtt.tif": PET,
    "NDVItt.tif": ndvi
}

# Save each raster using the correct CRS and transform from the ref area
for name, array in rasters.items():
    path = os.path.join(output, name)
    height, width = array.shape
    dtype = array.dtype

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=dtype,
        crs=ref_crs,               # Use original CRS (from reference band)
        transform=ref_transform # Use transform from original raster
    ) as dst:
        dst.write(array, 1)
