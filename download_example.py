from Sentinel3A_Download import SentinelLST
import geopandas

us_name, pasword = "usname", "password"
download_loc = "download_loc"
period = ("2022-04-12", "2023-01-13")

shp_area = "shape area"
gdf = geopandas.read_file(shp_area)
w, s = gdf.bounds[["minx", "miny"]].min()
e, n = gdf.bounds[["maxx", "maxy"]].max()

# Change area list depends on your area of interest (w,s,e,n) if there is no shapefile.
area = [w, s, e, n]  # Rectangle of the study area (west,south,east,north)
# time line can be "Non Time Critical" or "Near Real Time"
lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area, timeline="Non Time Critical")
lst.download_all()
