import time

from Sentinel3A_Download import SentinelLST
import geopandas
import os

us_name, pasword = "onursahin", "4M8CVNGK"
download_loc = "E:\LST_Satellites\Sentinel_3A\WAB\\2021"
period = ("2021-04-19", "2021-04-20")

shp_area = os.path.join(os.getcwd(), "WAB", "WAB_area.shp")
print(shp_area)
gdf = geopandas.read_file(shp_area)
w, s = gdf.bounds[["minx", "miny"]].min()
e, n = gdf.bounds[["maxx", "maxy"]].max()

# Change area list depends on your area of interest (w,s,e,n) if there is no shapefile.
area = [w, s, e, n]  # Rectangle of the study area (west,south,east,north)
# timeline can be "Non-Time Critical" or "Near Real Time"
for i in range(100):
    print("--------------")
    print(f"attemp {i}")

    lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area, timeline="Non Time Critical")
    lst.download_all()
    time.sleep(1200)
