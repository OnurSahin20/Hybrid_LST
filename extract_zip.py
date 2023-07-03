from Sentinel3A_Download import SentinelLST
import os

direc = "E:\LST_Satellites\Sentinel_3A\WAB\\2020_zips\\"
unzip_folder = "E:\LST_Satellites\Sentinel_3A\WAB\\2020_extracted\\"
days = sorted(os.listdir(direc))
print(days)
for day in days:
    path = os.path.join(direc + day)
    SentinelLST.extract_zip(path,unzip_folder)
