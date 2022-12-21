# Hybrid_LST
Hybrid_LST is a new initiative to create a generic program to visualize, manipulate, and hybridize satellite-based Land Surface Temperature Data. As a result, spatially and temporally continuous LST maps can be automatically obtained for basin or regional scale using the generic program.  
First version of the program supports downloading Sentinel 3A L2 LST products and categorizes data with time and Orbit Direction (Descending, Ascending). Afterward, the code automatically extracts NetCDF files from the zip for parent (time) folders.

### **Getting Started**

```
from Download_SentinelL2_LST.py import SentinelLST

us_name, pasword = "guest", "guest"  # authentication of sentinel open access hub.
download_loc = "D:\\Sentinel_3A"  # full path of download location of the LST products.
period = ("2022-12-02", "2022-12-03")  # [begin_date,end_date) ex.("2022-12-02", "2022-12-03") listed 1 day
area = [26, 36, 28, 38]  # Rectangle of the study area (west,south,east,north)  

# Spatial bounds also can be implemented using shapefile of study area (Later release!)
# Create Class SentinelLST

lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area) # It interacts with Copernicus Open Access Hub

lst.download_all()

lst.extract_zip()  # Extracting data from the zip files.
```
#### Code Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/messages.PNG?raw=true)

#### Folder Outputs
![solarized palettes](https://github.com/OnurSahin20/Hybrid_LST/blob/main/loc.PNG?raw=true)
