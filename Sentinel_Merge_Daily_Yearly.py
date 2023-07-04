import numpy as np
import os
import netCDF4

direc = "E:\\WAB_LST\\Sentinel_Daily\\2022\\"


def daily_to_yearly(in_direc, out_direc):
    all_files = os.listdir(in_direc)
    dfiles = sorted(list(filter(lambda x: "D" == x.split("_")[-2], all_files)))
    afiles = sorted(list(filter(lambda x: "A" == x.split("_")[-2], all_files)))
    for files in [dfiles, afiles]:
        c = 0
        fill_rate = []
        for file in files:
            if c % 10 == 0:
                print(file)

            fill_rate.append(file.split("_")[-1].split(".")[0])
            nc_file = netCDF4.Dataset(in_direc + file)
            if c == 0:
                lat, lon = np.array(nc_file["lat"]), np.array(nc_file["lon"])
                t, row, col = len(files), lat.shape[0], lon.shape[0]
                yearly_data = np.ones((t, row, col)) * -9999

            nc_data = np.array(nc_file["lst"])
            yearly_data[c] = nc_data
            c += 1

        yearly_data[yearly_data == -9999.0] = np.nan
        direction = file.split(".")[0].split("_")[-2]
        ds = netCDF4.Dataset(out_direc + "\\" + "WAB_Sentinel_LST_Yearly" + direction + ".nc", "w", format='NETCDF4')
        ds.createDimension('lat', row)
        ds.createDimension('lon', col)
        ds.createDimension('time', t)
        yy = ds.createVariable('lat', 'f4', ('lat',))
        xx = ds.createVariable('lon', 'f4', ('lon',))
        tt = ds.createVariable('time', 'i4', ('time',))
        yy[:] = lat
        xx[:] = lon
        tt[:] = [_ for _ in range(0, t)]
        values = ds.createVariable("lst", 'f4', ("time", 'lat', 'lon',), fill_value=-9999)
        values.unit = "celsius"
        values[:, :, :] = yearly_data
        ds.close()

    return "netcdf file is written"


if __name__ == '__main__':
    daily_to_yearly(direc, os.getcwd())
