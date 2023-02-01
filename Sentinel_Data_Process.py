import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import pandas as pd
from shapely.geometry import Point
import geopandas
from matplotlib.path import Path


class SentinelLST:

    def __init__(self, folder_loc, shpfile="", rec=()):
        self.folder_loc, self.shpfile, self.rec = folder_loc, shpfile, rec
        if len(self.shpfile) == 0 and len(self.rec) == 0:
            raise ValueError('One of the variable "shpfile" or "rec" must be filled')
        self.sen_lst_file, self.sen_coord = "LST_in.nc", "geodetic_in.nc"
        self.sen_cartesian, self.flag, self.met = "cartesian_in.nc", "flags_in.nc", "met_tx.nc"
        self.lst_file = self.folder_loc + "\\" + self.sen_lst_file
        self.coord_file = self.folder_loc + "\\" + self.sen_coord
        self.cartesian_file = self.folder_loc + "\\" + self.sen_cartesian
        self.flag_file = self.folder_loc + "\\" + self.flag
        self.met_file = self.folder_loc + "\\" + self.met

    def shp_file_masking(self, lat, lon):

        if len(self.shpfile) == 0:
            raise ValueError('"shapefile" must be declared')
        else:
            gridx, gridy = np.meshgrid(lon, lat)
            gdf = geopandas.read_file(self.shpfile)
            r, c = gridx.shape
            poly = list(gdf.unary_union.exterior.coords)
            poly_path = Path(poly)
            coors = np.hstack((gridx.reshape((r * c, 1)), gridy.reshape((r * c, 1))))
            mask = poly_path.contains_points(coors).reshape((r, c))
        return mask.astype(int)

    @staticmethod
    def access_nc_file(file):
        # method opens the .nc extension file.
        nc_file = netCDF4.Dataset(file)
        return nc_file

    def sentinel_lst_data(self, flags_in=True, uncertainty_mask=True,
                          exception_mask=True, uncertainty=2):
        """ Method extracts LST from Sentinel 3A data from the folder location.Two type flag method is available.
        First one is the "flags_in" which converts pixel values nan if pixel labeled with any flag inside "flags_in.nc"
        Second is "met_tx" which mask values higher than "the cloud fraction" values. Predefined is "flags_in strictly
        masks if data labeled as cloud contaminated"""
        nc_file, lst_var, exc_var = self.access_nc_file(self.lst_file), "LST", "exception"
        nc_data, exc_data = nc_file[lst_var], nc_file[exc_var]
        flag_mask, exception = exc_data.flag_masks, np.array(exc_data)
        lst_data = np.array(nc_data)
        lst_data[lst_data == nc_data._FillValue] = np.nan
        if exception_mask:
            for flag in flag_mask:
                lst_data[exception == int(flag)] = np.nan
        if flags_in:
            masks = self.flags_in()
            for mask in masks:
                lst_data[mask != 0] = np.nan
        if uncertainty_mask:
            mask = self.uncertainty_mask(uncertainty=uncertainty)
            lst_data[mask] = np.nan
        return lst_data - 272.15  # it converts Kelvin to Celsius

    def flags_in(self):
        nc_file, f1, f2, f3 = self.access_nc_file(self.flag_file), "cloud_in", "pointing_in", "bayes_in"
        return [np.array(nc_file[f1]), np.array(nc_file[f2]), np.array(nc_file[f3])]

    def extract_wgs_coord(self):
        """Method extracts WGS84 coordinates from Sentinel 3A data from the folder location"""
        lat_var, lon_var, nc_file = "latitude_in", "longitude_in", self.access_nc_file(self.coord_file)
        lat, lon = np.array(nc_file[lat_var]).astype(np.float32), np.array(nc_file[lon_var]).astype(np.float32)
        return lat, lon

    def extract_cartesian_coord(self):
        """Method extracts cartesian coordinates"""
        y_var, x_var, nc_file = "y_in", "x_in", self.access_nc_file(self.cartesian_file)
        y, x = np.array(nc_file[y_var]), np.array(nc_file[x_var])
        yfil, xfil = nc_file[y_var]._FillValue, nc_file[x_var]._FillValue
        y[y == yfil] = np.nan
        x[x == xfil] = np.nan
        return y, x

    def visualize_data_point(self, df):
        """ Method transforms arrays to data frame and visualize the land surface temperature as points"""
        geometry = [Point(xy) for xy in zip(df["lon"], df["lat"])]
        geo_df = geopandas.GeoDataFrame(df, geometry=geometry)
        if len(self.shpfile) > 0:
            fig, ax = plt.subplots(1, 1, figsize=(10,6))
            gdf = geopandas.read_file(self.shpfile)
            base = geo_df.plot(column="lst", legend=True, ax=ax)
            gdf.plot(ax=ax, color="white", edgecolor="black", alpha=0.1)
        else:
            base = geo_df.plot(column="lst", legend=True)
        plt.ylabel("Latitude"), plt.xlabel("Longitude")
        plt.show()

    def uncertainty_mask(self, uncertainty):
        """The method determines the quality data of given data set based on the threshold which masks."""
        lst_unc = "LST_uncertainty"  # LST uncertainty variable name in the sen_lst_file
        nc_data = self.access_nc_file(self.lst_file)[lst_unc]
        fill = nc_data._FillValue
        lst_unc = np.array(nc_data)
        lst_unc[lst_unc == fill] = np.nan
        return np.abs(lst_unc) > uncertainty

    def boundary_mask(self, data):
        """lst_data could be the raw data or the quality checked LST depends on the users requirements.
        If shape_mask is False, method uses rectangle coordinates (west,south,east,north).
        If shape_mask is True specify full path of the shapefile file with extension .shp to mask the data
        Method is the first basic rectangular masking
        """
        if len(self.rec) == 0:
            gdf = geopandas.read_file(self.shpfile)
            w, s = gdf.bounds[["minx", "miny"]].min()
            e, n = gdf.bounds[["maxx", "maxy"]].max()
        else:
            if len(self.rec) < 4:
                raise ValueError("input tuple has to be (west,south,east,north)")
            w, s, e, n = self.rec
        lat, lon = self.extract_wgs_coord()
        # y, x = self.extract_cartesian_coord() # later!!
        # Actual equivalent of lat,lon can be found using cartesian coordinates x,y (Later implementation)
        r, c = lat.shape
        df = pd.DataFrame({"lon": lon.reshape((r * c)), "lat": lat.reshape((r * c)),
                           "lst": data.reshape((r * c))}).dropna()
        df_mask = df.loc[(df.lat >= s) & (df.lat <= n) & (df.lon >= w) & (df.lon <= e)]
        return df_mask

    def re_griding(self, mask_df, res=0.01):
        from scipy.interpolate.interpnd import _ndim_coords_from_arrays
        from scipy.spatial import KDTree
        from scipy.interpolate import NearestNDInterpolator
        """ Only have one interpolation option (nearest) and it can be improved (RBF implementation)"""
        if len(self.rec) == 0:
            gdf = geopandas.read_file(self.shpfile)
            w, s = gdf.bounds[["minx", "miny"]].min()
            e, n = gdf.bounds[["maxx", "maxy"]].max()
        else:
            if len(self.rec) < 4:
                raise ValueError("input tuple has to be (west,south,east,north)")
            w, s, e, n = self.rec

        x, y = np.arange(w, e + res, res), np.arange(n, s - res, -res)
        grid, data_points = np.meshgrid(x, y), mask_df.loc[:, ["lon", "lat"]].values
        tree, xi = KDTree(data_points), _ndim_coords_from_arrays(tuple(grid), ndim=data_points.shape[1])
        dists, indexes = tree.query(xi)
        interp = NearestNDInterpolator(list(zip(data_points[:, 0], data_points[:, 1])), mask_df["lst"].values)
        model = interp(*grid)
        dists[dists > np.sqrt(2) * res] = 0
        dists[dists != 0] = 1
        model[dists.astype(int) == 0] = np.nan
        return x, y, model
    
    def save_to_tiff(self, x, y, grid_data, save_loc, product="Day"):
        import rasterio
        from rasterio.transform import Affine
        xres, yres = x[1] - x[0], y[0] - y[1]
        crs = "epsg:4326"
        """save_loc is location where tiff file will save. Products is day (Day) or night (night)"""
        fname = save_loc + "\\" +  "SLSTR_LST_" + self.folder_loc.split("____")[-1].split("T")[0] \
                + "_" + product + ".tif"
        transform = Affine.translation(x[0] - xres / 2, y[-1] - yres / 2) * Affine.scale(xres, yres)
        try:
            with rasterio.open(
                    fname,
                    mode="w",
                    driver="GTiff",
                    height=grid_data.shape[0],
                    width=grid_data.shape[1],
                    count=1,
                    dtype=grid_data.dtype,
                    crs=crs,
                    transform=transform,
                    nodata = -9999
            ) as new_dataset:
                new_dataset.write(grid_data[::-1], 1)
            return "LST data saved to tif file successfully"

        except:
            return "tiff file can't save !!!"

if __name__ == '__main__':
    lst_path = "LST path"
    shape_path = "shapefile_path"
    sentinel_class = SentinelLST(lst_path, shpfile=shape_path)
    lst_data = sentinel_class.sentinel_lst_data(flags_in=True, uncertainty_mask=True, exception_mask=True,
                                                uncertainty=2)
    lst_df = sentinel_class.boundary_mask(lst_data)
    lon, lat, lst_raster = sentinel_class.re_griding(lst_df, res=0.01)
    mask = sentinel_class.shp_file_masking(lat, lon)
    lst_raster[mask == False] = np.nan
    plt.pcolormesh(lon, lat, lst_raster)
    plt.show()
    sentinel_class.save_to_tiff(lon,lat,lst_raster,"loc","day")
