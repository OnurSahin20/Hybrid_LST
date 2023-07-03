from Sentinel_Extract_Daily import SentinelDaily
import os

shp_area = "C:\\Users\\gungor-sahin\\PycharmProjects\\Hybrid_LST_modified\\WAB\\WAB_area.shp"

fold_path = "E:\\LST_Satellites\\Sentinel_3A\WAB\\2022_extracted\\"
output_dir = os.getcwd()
folders = sorted(os.listdir(fold_path))

for fold in folders:
    sentinel = SentinelDaily(fold_path + fold, shp_area, output_dir, 0.08)
    sentinel.merge_write_instances()

