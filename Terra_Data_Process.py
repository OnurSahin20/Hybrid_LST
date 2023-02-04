from numpy import array
from numpy import nan
from netCDF4 import Dataset


class TerraLST:
    def __init__(self, full_nc_path, uncertainty=2):
        self.full_nc_path = full_nc_path
        self.lst_varibs = {"day": "LST_Day_1km", "night": "LST_Night_1km"}
        self.qc_varibs = {"day": "QC_Day", "night": "QC_Night"}
        self.view_time = {"day": "Day_view_time", "night": "Night_view_time"}
        self.uncertanity = uncertainty

    def get_wgs_coord(self):
        nc_file = Dataset(self.full_nc_path)
        lat, lon = array(nc_file["lat"]), array(nc_file["lon"])
        return lat, lon

    def get_lst(self, product="day"):
        if product not in ["day", "night"]:
            raise ValueError('Only products "day" and "night" parameters for MODIS LST')

        nc_file = Dataset(self.full_nc_path)
        lst, qc = array(nc_file[self.lst_varibs[product]]), array(nc_file[self.qc_varibs[product]])
        lst[lst == nc_file[self.lst_varibs[product]]._FillValue] = nan
        flags = []
        if self.uncertanity == 2:
            flags = [129, 145]
        elif self.uncertanity == 1:
            flags = [65, 73, 81, 129, 145]
        for flag in flags:
            lst[qc == flag] = nan
        return lst - 273.15


if __name__ == '__main__':
    terra_path = "KMB_LST_Values_August.nc"
    terra = TerraLST(terra_path)
    lat, lon = terra.get_wgs_coord()
    terra_lst = terra.get_lst()
    print(terra_lst)
