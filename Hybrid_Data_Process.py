import netCDF4
import numpy as np
from pandas import date_range
import os
from Terra_Data_Process import TerraLST
from Sentinel_Data_Process import SentinelLST


class HybridLst:
    def __init__(self, terra_nc_path, sentinel_direc, time_period, shp_file=""):
        self.shp_file = shp_file
        self.time_period = list((date_range(time_period[0], time_period[1])).strftime('%Y-%m-%d'))[0:-1]
        self.terra_nc_path = terra_nc_path
        self.sentinel_direc = sentinel_direc
        self.terra_class = TerraLST(self.terra_nc_path)
        self.terra_lat, self.terra_lon = self.terra_class.get_wgs_coord()
        self.all_terra_lst_day = self.terra_class.get_lst(product="day")
        self.all_terra_lst_night = self.terra_class.get_lst(product="night")
        self.mask_file = "spatial_mask_array.txt"

        if self.mask_file not in os.listdir(os.getcwd()):
            self.masking_array = SentinelLST.shp_file_masking(self.shp_file, self.terra_lat, self.terra_lon)
            # how many empty cell due to spatial masking.
            np.savetxt(self.mask_file, self.masking_array, fmt="%10.0f")
        else:
            self.masking_array = np.loadtxt(self.mask_file)

        self.count_empty = np.sum(self.masking_array == False)

    def get_sentinel_instant_lst(self, cur_date, date="D"):
        full_path = self.sentinel_direc + "\\Sentinel_LST_" + cur_date
        empty_data = np.zeros((len(self.terra_lat), len(self.terra_lon))) * np.nan
        if "Sentinel_LST_" + cur_date not in os.listdir(self.sentinel_direc):
            return empty_data

        all_files = os.listdir(full_path)
        zip_filter = [file for file in all_files if ".zip" not in file]  # it filter .zip files and correct products.
        product_filter = [file for file in zip_filter if file.split("_")[-1] == date]
        if len(product_filter) == 0:
            return empty_data
        lst_data = np.zeros((len(product_filter), len(self.terra_lat), len(self.terra_lon))) * np.nan
        # scan_time = "No Data"
        for f in range(len(product_filter)):
            # scan_time = product_filter[f].split("_")[7].split("T")[-1]
            sen_cla = SentinelLST(full_path + "\\" + product_filter[f], self.shp_file)
            lst_df = sen_cla.boundary_mask(sen_cla.sentinel_lst_data())
            if len(lst_df) == 0:
                lst_data[f, :, :] = empty_data
            else:
                x, y, model = sen_cla.re_griding(self.terra_lon, self.terra_lat)
                lst_data[f, :, :] = model
        if len(product_filter) > 1:
            if np.all(np.isnan(lst_data[1, :, :])):
                lst = lst_data[0, :, :]
                lst[self.masking_array == False] = np.nan
                return lst
            else:
                lst_mean = np.nanmean(lst_data, axis=0)
                lst_mean[self.masking_array == False] = np.nan
                return lst_mean
        else:
            lst = lst_data[0, :, :]
            lst[self.masking_array == False] = np.nan
            return lst_data[0, :, :]

    def get_terra_instant_lst(self, cur_date, date="day"):
        index = self.time_period.index(cur_date)
        if date == "day":
            terra_lst = self.all_terra_lst_day[index, :, :]
        elif date == "night":
            terra_lst = self.all_terra_lst_night[index, :, :]
        else:
            raise ValueError("Only terra and slstr instance ara avaliable right now !!!!")
        terra_lst[self.masking_array == False] = np.nan
        return terra_lst

    def build_time_series(self, product, date):
        lst_data = np.zeros((len(self.time_period), len(self.terra_lat), len(self.terra_lon))) * np.nan
        for t in range(len(self.time_period)):
            if product == "terra":
                lst_data[t, :, :] = self.get_terra_instant_lst(self.time_period[t], date)
            elif product == "slstr":
                lst_data[t, :, :] = self.get_sentinel_instant_lst(self.time_period[t], date)
        return lst_data

    def save_lsts_netcdf(self, loc, nc_name):
        ds = netCDF4.Dataset(loc + "\\" + nc_name + ".nc4", "w", format='NETCDF4')
        ds.createDimension('lat', len(self.terra_lat))
        ds.createDimension('lon', len(self.terra_lon))
        ds.createDimension('time', len(self.time_period))
        x, y = ds.createVariable('lon', 'f4', ('lon',)), ds.createVariable('lat', 'f4', ('lat',))
        t = ds.createVariable('time', 'f4', ('time',))
        x[:], y[:], t[:] = self.terra_lon, self.terra_lat, [_ for _ in range(len(self.time_period))]
        for product in ["terra", "slstr"]:
            for d1, d2 in [("D", "day"), ("A", "night")]:
                if product == "terra":
                    name = product + "_" + d2
                    values = ds.createVariable(name, 'f4', ("time", 'lat', 'lon',), fill_value=-9999)
                    values[:, :, :] = self.build_time_series(product, d2)
                else:
                    name = product + "_" + d2
                    values = ds.createVariable(name, 'f4', ("time", 'lat', 'lon',), fill_value=-9999)
                    values[:, :, :] = self.build_time_series(product, d1)
        ds.close()
        return "day and night LSTs save to .nc4 file successfully."

if __name__ == '__main__':
    sentinel_path = "D:\\LST_Satellites\\Sentinel_3A\\KMB\\August_2021"
    shpp_file = "C:\\Users\\Gungor\\PycharmProjects\\sentinel_api\\KMB\\KMB.shp"
    terra_path = "KMB_LST_Values_August.nc"
    hybrid = HybridLst(terra_path, sentinel_path, ("2022-08-01", "2022-09-01"), shpp_file)
