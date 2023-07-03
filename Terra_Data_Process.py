import numpy as np
import netCDF4


class TerraLST:
    def __init__(self, full_nc_path, uncertainty=2):
        self.full_nc_path = full_nc_path
        self.lst_varibs = {"day": "LST_Day_1km", "night": "LST_Night_1km"}
        self.qc_varibs = {"day": "QC_Day", "night": "QC_Night"}
        self.view_time = {"day": "Day_view_time", "night": "Night_view_time"}
        self.uncertanity = uncertainty
        self.data_set = netCDF4.Dataset(self.full_nc_path)

    def get_wgs_coord(self):
        lat, lon = np.array(self.data_set["lat"]), np.array(self.data_set["lon"])
        return lat, lon

    def get_time(self):
        return np.array(self.data_set["time"])

    def get_lst(self, product="day"):
        if product not in ["day", "night"]:
            raise ValueError('Only products "day" and "night" parameters for MODIS LST')

        lst, qc = np.array(self.data_set[self.lst_varibs[product]]), np.array(self.data_set[self.qc_varibs[product]])
        lst[lst == self.data_set[self.lst_varibs[product]]._FillValue] = np.nan
        flags = []
        if self.uncertanity == 2:
            flags = [129, 145]
        elif self.uncertanity == 1:
            flags = [65, 73, 81, 129, 145]
        for flag in flags:
            lst[qc == flag] = np.nan
        return lst - 273.15

