import os
import numpy as np
from Sentinel_Data_Process import SentinelLST
import geopandas
import netCDF4


class SentinelDaily:
    def __init__(self, fold_path, shape_path, output_dir, res):
        self.fold_path = fold_path  # path to day folder!
        self.shape_path = shape_path  # shape file!
        self.gdf = geopandas.read_file(self.shape_path)
        self.res = res
        self.e, self.n = self.gdf.bounds[["maxx", "maxy"]].max()
        self.w, self.s = self.gdf.bounds[["minx", "miny"]].min()
        self.x = np.arange(self.w, self.e + self.e / 2, self.res, dtype=float)
        self.y = np.arange(self.n, self.s - self.s / 2, -res, dtype=float)
        self.r, self.c = self.y.shape[0], self.x.shape[0]
        self.mask = self.shp_file_masking()
        self.nans = np.sum(self.mask == 0)  # constant nans (sea)
        self.count_pixels = self.r * self.c - self.nans  # total pixel numbers for area of interest
        self.output_dir = output_dir

    def shp_file_masking(self):
        from matplotlib.path import Path
        if len(self.shape_path) == 0:
            raise ValueError('"shapefile" must be declared')
        else:
            gridx, gridy = np.meshgrid(self.x, self.y)
            gdf = geopandas.read_file(self.shape_path)
            r, c = gridx.shape
            poly = list(gdf.unary_union.exterior.coords)
            poly_path = Path(poly)
            coors = np.hstack((gridx.reshape((r * c, 1)), gridy.reshape((r * c, 1))))
            mask = poly_path.contains_points(coors).reshape((r, c))
        return mask.astype(int)

    def lst_instance(self, file):
        path = os.path.join(self.fold_path, file)
        sentinel = SentinelLST(path, self.shape_path)
        lst = sentinel.sentinel_lst_data(exception_mask=True, uncertainty=2)
        lst_df = sentinel.boundary_mask(lst)
        lst_data = sentinel.re_griding(self.y, self.x, lst_df, res=self.x[1] - self.x[0])
        lst_data[self.mask == 0] = np.nan
        return lst_data

    def get_lst_am_pm(self, product="D"):
        files = os.listdir(self.fold_path)
        list_files = [file for file in files if file.split("_")[-1] == product and ".zip" not in file]
        if len(list_files) == 0:
            lst = np.zeros((self.r, self.c)) * np.nan
        elif len(list_files) == 1:
            lst = self.lst_instance(list_files[0])
        else:
            lst_data = np.zeros((len(list_files), self.r, self.c))
            for d in range(len(list_files)):
                lst_data[d] = self.lst_instance(list_files[d])

            if np.all(np.isnan(lst_data)):
                lst = np.zeros((self.r, self.c)) * np.nan
            else:
                lst = np.nanmean(lst_data, axis=0)
        return lst

    def merge_write_instances(self):
        lst_d = self.get_lst_am_pm(product="D")
        self.write_nc(lst_d, "D")
        lst_a = self.get_lst_am_pm(product="A")
        self.write_nc(lst_a, "A")
        day = self.fold_path.split("_")[-1]
        print(f"Descending and ascending lst maps are written for the day {day}")
        return 1

    def write_nc(self, data, product="D"):
        percent = int(np.sum(np.isnan(data) == 0) / self.count_pixels * 100 + 0.7)
        nc_name = self.fold_path.split("_")[-1] + "_" + product + "_P" + str(percent)
        lat, lon = self.y.copy(), self.x.copy()
        ds = netCDF4.Dataset(self.output_dir + "\\" + "WAB_Sentinel_LST_" + nc_name + ".nc4", "w", format='NETCDF4')
        ds.createDimension('lat', self.r)
        ds.createDimension('lon', self.c)
        yy = ds.createVariable('lat', 'f4', ('lat',))
        xx = ds.createVariable('lon', 'f4', ('lon',))
        yy[:] = lat
        xx[:] = lon
        values = ds.createVariable("lst", 'f4', ('lat', 'lon',), fill_value=-9999)
        values.unit = "celsius"
        values[:, :] = data
        ds.close()
        return 1
