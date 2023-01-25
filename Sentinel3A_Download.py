import requests
import pandas as pd
from bs4 import BeautifulSoup
import os
import geopandas
import webbrowser

""" Script is created by Onur Gungor Sahin
Institution : Izmir Institute of Technology
Department : International Water Resources
emails : onursahin@iyte.edu.tr, ogungorsahin@gmail.com
Script is first release version 0.1 for downloading Sentinel 3A Land Surface Temperature (LST) data.
This project extraction of LST for specific study area and post processing of it.
1) It creates parent folders (dates) for Sentinel 3A 
   and downloads each days LST to the right folders and labels them with pass direction (Ascending,Descending).
"""


class SentinelLST:

    def __init__(self, download_path, user_name, password, begindate, enddate, coordinates,
                 timeline="Non Time Critical", del_zips=True):
        self.down_path = download_path
        self.user_name = user_name
        self.password = password
        self.begin_date = begindate
        self.enddate = enddate
        self.timeline = timeline
        self.rec = coordinates
        self.url_start = "https://scihub.copernicus.eu/dhus/search?start=0&rows=100&q="
        self.platform = 'platformname:Sentinel-3'
        self.file_name = 'filename:S3A_SL_2_LST*'
        self.time_range = list((pd.date_range(self.begin_date, self.enddate)).strftime('%Y-%m-%d'))
        self.all_queries = self.open_search_query()
        self.offline_products = []

    def open_search_query(self):  # it returns to request get object.
        w, s, e, n = self.rec
        loc = f'footprint:"Intersects(POLYGON(({w} {s}, {e} {s}, {e} {n},{w} {n},{w} {s})))"'
        all_queries = []
        for i in range(len(self.time_range) - 1):
            date = 'beginposition:' + "[" + self.time_range[i] + "T00:00:00.000Z" + " TO " + self.time_range[
                i + 1] + "T00:00:00.000Z" + "]"
            all_queries.append(self.url_start + self.platform + " AND " + f"timeliness:{self.timeline}" + " AND " +
                               date + " AND " + loc + " AND " + self.file_name + '&orderby=ingestiondate desc')
        return all_queries

    def get_request(self, index):
        return requests.get(self.all_queries[index], auth=(self.user_name, self.password))

    def download_day(self, index) -> None:

        # above comment is creates child directories for descending and ascending products.
        req = self.get_request(index)  # getting request for the day.
        if req.status_code != 200:
            raise ValueError("Request to open access hub can't created. Check your queries")
        soup = BeautifulSoup(req.text, features="xml")
        direction = [str(tag).split(">")[1].split("<")[0] for tag in
                     soup.find_all("str", attrs={"name": "passdirection"})]  # ascending, descending directions.
        down_links = [tag.get("href") for tag in
                      soup.find_all("link", attrs={"rel": "alternative"})]  # download links from the query for day
        title_tags = soup.findAll("title")[1:]  # first title is definition and others are the name of the files
        print(f"{len(title_tags)} file is found for {self.time_range[index]}")
        for i in range(len(down_links)):
            link = down_links[i]
            name = title_tags[i].string
            online = requests.get(link + "Online/$value", auth=(self.user_name, self.password)).text
            if online == "true":
                par_dir = "Sentinel_LST_" + self.time_range[index]
                full_path = self.down_path + "\\" + par_dir  # direc for request day
                if par_dir not in os.listdir(self.down_path):
                    os.makedirs(full_path, exist_ok=True)  # it is creating parent directory for current day

                zip_file = name + "_" + direction[i][0].upper() + ".zip"
                if zip_file not in os.listdir(full_path):
                    with open(full_path + "\\" + zip_file, "wb") as file:
                        print(f"Download process of {name} begun. If it took to long check your internet speed")
                        response = requests.get(link + "$value", auth=(self.user_name, self.password))
                        file.write(response.content)
                    print("Done!")
                else:
                    print("File exist in the folder.")
            else:
                d1 = requests.get(link + "$value", auth=(self.user_name, self.password)).status_code
                if int(d1) == 202:
                    print(name + " is offline and triggered")
                    self.offline_products.append([name, link, direction[i][0].upper()])
                else:
                    self.offline_products.append([name, link, direction[i][0].upper()])

    def download_all(self):
        import csv
        for i in range(len(self.time_range) - 1):
            self.download_day(i)

        with open("Request_" + "OfflineFiles_" + self.time_range[0] + "-" + self.time_range[-2] + ".txt", "w",
                  newline="") as txt_file:
            writer = csv.writer(txt_file, delimiter=",")
            for off in self.offline_products:
                writer.writerow(off)
        print("All online files are downloaded and offline files are saved to txt")

    def extract_zip(self):
        import zipfile
        par_direcs = [self.down_path + "\\" + path for path in os.listdir(self.down_path)]  # folders to extract
        print("unzipping is started !")
        if len(par_direcs) == 0:
            print("There is no zip file to extract.")
            return -1
        for par in par_direcs:
            files = os.listdir(par)
            for z in files:
                if ".zip" in z:
                    with zipfile.ZipFile(par + "\\" + z, "r") as zip_obj:
                        zip_obj.extractall(path=par)

                pdirec = z.split(".")[0][-1]
                fname = z.split(".")[0][0:-2] + ".SEN3"
                os.rename(par + "\\" + fname, par + "\\" + fname + "_" + pdirec)
            break
        print("Files are extracted successfully !")


if __name__ == '__main__':
    us_name, pasword = "onursahin", "4M8CVNGK"  # authentication of sentinel open access hub.
    download_loc = "D:\\LST_Satellites\\Sentinel_3A"  # full path of download location of the LST products.
    period = ("2022-01-01", "2022-02-01")  # [begin_date,end_date) ex.("2022-12-02", "2022-12-03") listed 1 day

    shp_area = "shapefile location"
    gdf = geopandas.read_file(shp_area)
    w, s = gdf.bounds[["minx", "miny"]].min()
    e, n = gdf.bounds[["maxx", "maxy"]].max()
    area = [w, s, e, n]  # Rectangle of the study area (west,south,east,north)
    # It can be directly pass boundries of the study area without shapefile
    # Spatial bounds also can be implemented using shapefile of study area (Later release!)
    # time line can be "Non Time Critical" or "Near Real Time"
    lst = SentinelLST(download_loc, us_name, pasword, period[0], period[1], area, timeline="Non Time Critical")
    # lst.download_all()
    lst.extract_zip()
