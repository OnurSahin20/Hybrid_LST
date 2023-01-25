# Hybrid_LST
Hybrid_LST is a new initiative to create a generic program to visualize, manipulate, and hybridize satellite-based Land Surface Temperature Data. As a result, spatially and temporally continuous LST maps can be automatically obtained for basin or regional scale using the generic program.  
First version of the program supports downloading Sentinel 3A L2 LST products and categorizes data with time and Orbit Direction (Descending, Ascending). Afterward, the code automatically extracts NetCDF files from the zip for parent (time) folders. Code can handle offline files for specific time intervals. Each run code generates a txt file and saves offline files to it. Rerun your code a couple of times (1-2 hours later!) until the txt file is empty. 

### **Getting Started**

```
from Sentinel3A_Download.py import SentinelLST
us_name, pasword = "username", "password"  # authentication of sentinel open access hub.
download_loc = "location"  # full path of download location of the LST products.
period = ("2022-01-01", "2022-02-01")  # [begin_date,end_date) ex.("2022-12-02", "2022-12-03") listed 1 day

shp_area = "shapefile location"
gdf = geopandas.read_file(shp_area)
w, s = gdf.bounds[["minx", "miny"]].min()
e, n = gdf.bounds[["maxx", "maxy"]].max()
area = [w, s, e, n]  # Rectangle of the study area (west,south,east,north)
# It can be directly pass boundries of the study area without shapefile
# Spatial bounds also can be implemented using shapefile of study area (Later release!)
# time line can be "Non Time Critical" or "Near Real Time"
lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area, timeline="Non Time Critical")
# lst.download_all()
lst.extract_zip()
```
#### Code Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/messages.PNG?raw=true)

#### Folder Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/loc.PNG?raw=true)
