# Hybrid_LST
Hybrid_LST is a new initiative to create a generic program to visualize, manipulate, and hybridize satellite-based Land Surface Temperature Data. As a result, spatially and temporally continuous LST maps can be automatically obtained for basin or regional scale using the generic program.  
First version of the program supports downloading Sentinel 3A L2 LST products and categorizes data with time and Orbit Direction (Descending, Ascending). Afterward, the code automatically extracts NetCDF files from the zip for parent (time) folders. Code can handle offline files for specific time intervals. Each run code generates a txt file and saves offline files to it. Rerun your code a couple of times (1-2 hours later!) until the txt file is empty. 

### **Getting Started**

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
#### Code Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/xxxx.PNG?raw=true)


#### Folder Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/loc.PNG?raw=true)

### Mask, regrid, spatial mask, and visualize, Sentinel LST Data
L2 Sentinel SLSTR data is not conventional rasterized LST products and is not directly open using GIS programs. Sentinel_Data_Process.py can be utilized for masking, re-gridding, and visualization of Sentinel LST data. Inputs are the LST data folder location and shape file or rectangular coordinates of the study area. First you need to be create SentinelLst class. Quality definitions are given in "flags.nc" and "LST.nc" inside lst_path for Sentinel LST data.
```
from Sentinel_Data_Process.py import SentinelLST
lst_path = "lst_folder_loc"
shape_path = "shafile_loc"
sentinel_class = SentinelLST(lst_path,shpfile=shape_path)
# It creates one day and one ascending or descending order LST class.
lst_data = sentinel_class.sentinel_lst_data(flags_in=True, uncertainty_mask=True, exception_mask=True,uncertainty=2)
#  flags_in clean cloud contaminated data. Uncertainty masks control data quality using the threshold value "uncertainty."
lst_df = sentinel_class.boundary_mask(lst_data) # basically clip lst_data with boundaries of the study are output is pandas dataframe.
# Study are and LST data can be easliy visualize 
sentinel_class.visualize_data_point(lst_df)
```
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/visualize_plot.png?raw=true)
