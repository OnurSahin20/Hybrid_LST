from Sentinel3A_Download import SentinelLST

direc = "D:\\LST_Satellites\\Sentinel_3A\\WAB\\"

months = ["january", "february", "march", "april", "may", "june",
          "july", "august", "september", "october", "november", "december"]


folders = [direc + month + "_" + "2022" for month in months]
print(folders)
for fold in folders:
    SentinelLST.extract_zip(fold)
