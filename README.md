# Hybrid_LST
Hybrid_LST is a generic program to visualize, manipulate, and hybridize satellite-based Land Surface Temperature Data High (1 km) Resolution MODIS (Aqua - Terra),  VIIIRS, Sentinel3A- SLSTR and hyper (<100m) LANDSAT, ECOSTRESS. As a result, spatially and temporally continuous LST maps can be automatically obtained for basin or regional scale studies using this generic program.  
The program supports downloading Sentinel 3A L2 LST products and categorizes data with time and Orbit Direction (Descending, Ascending). Code can handle offline files. Code generates a txt file and saves offline files to it. Rerun your code a couple of times (1-2 hours later!) until the txt file is empty. 

### **Download Sentinel SLSTR LST Data**

```
from Sentinel3A_Download.py import SentinelLST

us_name, pasword = "username", "password"  # authentication of sentinel open access hub.
download_loc = "D:\\LST_Satellites\\Sentinel_3A"  # full path of download location of the LST products.
period = ("2023-01-01", "2023-01-11")  # [begin_date,end_date) ex.("2022-12-02", "2022-12-03") listed 1 day
# Rectangle of the study area (west,south,east,north)
"""  Shape file implementation  
shp_area = "shape_file_location"
gdf = geopandas.read_file(shp_area)
w, s = gdf.bounds[["minx", "miny"]].min()
e, n = gdf.bounds[["maxx", "maxy"]].max()
"""

# Change area list depends on your area of interest (w,s,e,n) if there is no shapefile.
area = [26, 36, 30, 40]  # Rectangle of the study area (west,south,east,north)
# time line can be "Non Time Critical" or "Near Real Time"
lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area, timeline="Non Time Critical")
lst.download_all()
lst.extract_zip()
```
<!-- #### Code Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/xxxx.PNG?raw=true)


#### Folder Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/loc.PNG?raw=true) -->

### Mask, regrid, spatial mask, and visualize, Sentinel LST Data
L2 Sentinel SLSTR data is not conventional rasterized LST products and is not directly open using GIS programs. Sentinel_Data_Process.py can be utilized for masking, re-gridding, and visualization of Sentinel LST data. Inputs are the LST data folder location and shape file or rectangular coordinates of the study area. First you need to be create SentinelLst class. Quality definitions are given in "flags.nc" and "LST.nc" inside lst_path for Sentinel LST data.
```
from Sentinel_Data_Process.py import SentinelLST
lst_path = "lst_folder_loc"
shape_path = "shafile_loc"
sentinel_class = SentinelLST(lst_path,shpfile=shape_path)
# It creates one day and one ascending or descending order LST class.  flags_in clean cloud contaminated data. Uncertainty masks control data quality using the threshold value "uncertainty."
lst_data = sentinel_class.sentinel_lst_data(flags_in=True, uncertainty_mask=True, exception_mask=True,uncertainty=2)
lst_df = sentinel_class.boundary_mask(lst_data) # basically clip lst_data with boundaries of the study are output is pandas dataframe.
sentinel_class.visualize_data_point(lst_df)
```
<!-- ![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/visualize_point.png?raw=true) -->

Data points can be used to interpolate and re-grid data for NetCDF or raster format. Code uses scipy. interpolate NearestNDInterpolator when cloud contamination high direct usage of the interpolated map can be problematic; hence if there are points inside some raster, those values interpolated. Geographic coordinate system is WGS84 and resolution of the grid data is optinal.
```
y,x = sentinel_class.extract_wgs_coord
lon, lat, lst_raster = sentinel_class.re_griding(y,x_df, res=0.01)
plt.pcolormesh(lon,lat,lst_raster)
plt.show()
```
<!-- ![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/regrid2.png?raw=true) -->
Shapefile can be used to mask regrid data. "shp_file_masking" use matplotlib.path import Path for pixel masking. It returns 2d the boolean numpy array.
```
mask = sentinel_class.shp_file_masking(lat, lon)
lst_raster[mask == False] = np.nan
plt.pcolormesh(lon, lat, lst_raster)
plt.show()
```
<!-- ![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/mask.png?raw=true) -->

Re-grid data can be saved to .tif file.
```
sentinel_class.save_to_tiff(lon,lat,lst_raster,save_loc = "loc", product = "day") 
# save_loc where tif file will be saved!
# Naming covention is "SLSTR_LST" + time + "day or night" + .tif file.
```
### Terra LST class.
It is easy to get MODIS LST products using AρρEEARS - NASA tool. TerraLST class access lat, lon coordinates and LST data. Variables [""LST_Day_1km",""LST_Night_1km","QC_Day","QC_Night"] have to be inside Terra nc file.(It supoorts only netcdf file!). Quality variables of both day, night products provide LST errors which are [0 - Average LST error <= 1K, (65,73,81) - Average LST error <= 2K and (129,145) - Average LST error <= 3K). Use uncertainty input params to choose your quality flag. Other values will be masked. Straigthfoward usage of TerraClass 
```
from Terra_Data_Process.py import TerraLST
terra_path = "Full path of Terra L3 product is taken from AρρEEARS"
terra = TerraLST(terra_path)
lat, lon = terra.get_wgs_coord() # arrays of latitudes and longitudes.
terra_lst = terra.get_lst()     # 3d numpy array (time, lat, lon)
```


