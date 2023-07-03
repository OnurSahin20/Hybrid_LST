import requests
import pandas as pd
from bs4 import BeautifulSoup
import os
import geopandas


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
        self.logic = False

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
            raise ValueError("Request to open access hub can't created. Check your authentication")
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
            # if name.split("_")[-3] == "O":
            #     continue
            online = requests.get(link + "Online/$value", auth=(self.user_name, self.password)).text
            par_dir = "Sentinel_LST_" + self.time_range[index]
            full_path = self.down_path + "\\" + par_dir  # direc for request day
            zip_file = name + "_" + direction[i][0].upper() + ".zip"

            if par_dir not in os.listdir(self.down_path):
                os.makedirs(full_path, exist_ok=True)  # it is creating parent directory for current day

            if zip_file in os.listdir(full_path):
                if os.path.getsize(full_path + "\\" + zip_file) == 0:
                    os.remove(full_path + "\\" + zip_file)
                    print(f"Empty zip file  = {zip_file} removed !.")
            if online == "true":
                if zip_file not in os.listdir(full_path):
                    print(f"Download process of {name} begun. If it took to long check your internet speed")
                    with open(full_path + "\\" + zip_file, "wb") as file:
                        final = link + "$value"
                        response = requests.get(final, auth=(self.user_name, self.password))
                        if response.status_code == 200:
                            success = True
                            file.write(response.content)
                            print("Done!")
                        else:
                            success = False
                            print("File is online but it can't be downloaded - HTTP ERROR 500!")
                            self.offline_products.append([name, link, direction[i][0].upper()])
                    if not success:
                        os.remove(full_path + "\\" + zip_file)
                else:
                    print("File exist in the folder!")
            else:
                if zip_file not in os.listdir(full_path):
                    d1 = requests.get(link + "$value", auth=(self.user_name, self.password)).status_code
                    if int(d1) == 202:
                        print(name + " is offline and triggered")
                        self.offline_products.append([name, link, direction[i][0].upper()])
                    elif int(d1) == 403:
                        print("The request is not accepted because the number of submitted "
                              "requested exceeded the allowed user quota")
                        self.offline_products.append([name, link, direction[i][0].upper()])
                        self.logic = True

                    else:
                        print("HTTP ERROR 500!")
                        self.offline_products.append([name, link, direction[i][0].upper()])
                else:
                    print("File is already exist in folder that is why is not triggered")

    def download_all(self):
        for i in range(len(self.time_range) - 1):
            if self.time_range[i] == "2020-01-27":
                dummy = 5
            self.download_day(i)
            if self.logic:
                print("-------------------")
                print("Requested exceeded the allowed user quota. That's why process is stopped. Run this script"
                      " later and repeat several times!")
                break
        # with open("Request_" + "OfflineFiles_" + self.time_range[0] + "-" + self.time_range[-2] + ".txt", "w",
        #           newline="") as txt_file:
        #     writer = csv.writer(txt_file, delimiter=",")
        #     for off in self.offline_products:
        #         writer.writerow(off)
        # print("All online files are downloaded and offline files are saved to txt")

    @staticmethod
    def extract_zip(input_path, unzip_path):
        import zipfile
        files = os.listdir(input_path)
        par_direcs = [input_path + "\\" + file for file in files]
        unzip_name = unzip_path + "Sentinel_LST_" + input_path.split("_")[-1]
        day = input_path.split("_")[-1]
        if "Sentinel_LST_" + day not in os.listdir(unzip_path):
            os.mkdir(unzip_name)
        if len(files) > 0:
            print(f"day = {day} and unzipping is started !")
            if len(par_direcs) == 0:
                print("There is no zip file to extract.")
                return -1
            c = 0
            for par in par_direcs:
                if ".zip" in par:
                    pdirec = par.split(".")[0][-1]
                    fname = files[c].split(".")[0][0:-2] + ".SEN3"
                    fold_name = fname + "_" + pdirec
                    if fold_name not in os.listdir(unzip_name):
                        with zipfile.ZipFile(par, "r") as zip_obj:
                            print(f"Zip files in {par} are extracting right now!")
                            zip_obj.extractall(path=unzip_name)
                            os.rename(unzip_name + "\\" + fname, unzip_name + "\\" + fold_name)
                    else:
                        print("Zip file is extracted before !")
                c += 1
