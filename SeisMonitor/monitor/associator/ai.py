import numpy as np
import pandas as pd
import datetime as dt
import pyproj
from . import utils as ut
from gamma.utils import association
from obspy.core.inventory.inventory import (Inventory,read_inventory)
from tqdm import tqdm


class GaMMAObj():
    def __init__(self,response,region,epsg_proj,
                use_dbscan=True,use_amplitude=True,
                dbscan_eps=10.0,dbscan_min_samples=3,
                vel = {"p": 7.0, "s": 7.0 / 1.75},
                method="BGMM",oversample_factor=20,
                min_picks_per_eq=5,max_sigma11=2.0,
                max_sigma22=1.0,max_sigma12=1.0,
                calculate_amp=True,p_window=10,
                s_window=5,waterlevel=10
                ):
        if isinstance(response,Inventory):
            self.response = response
        else:
            self.response = read_inventory(response)

        self.lon_lims = region[0:2]
        self.lat_lims = region[2:4]
        self.z_lims = region[4:]
        self.epsg_proj = epsg_proj
        self.use_dbscan = use_dbscan
        self.use_amplitude = use_amplitude
        self.vel = vel
        self.dbscan_eps = dbscan_eps
        self.dbscan_min_samples = dbscan_min_samples
        self.method = method
        self.oversample_factor = oversample_factor
        self.min_picks_per_eq = min_picks_per_eq
        self.max_sigma11 = max_sigma11
        self.max_sigma22 = max_sigma22
        self.max_sigma12 = max_sigma12
        self.dims = ["x(km)", "y(km)", "z(km)"]
        self.calculate_amp=calculate_amp
        self.p_window=p_window
        self.s_window=s_window
        self.waterlevel=waterlevel
        self.config = self._get_config()

    # def _get_config(self):

    #     config = self.__dict__
    #     config["degree2km"] = 111.195
    #     config["center"] = [np.mean(self.lon_lims), np.mean(self.lat_lims)]
    #     config["x(km)"] = ((np.array(self.lon_lims) - config["center"][0]) * config["degree2km"]).tolist()
    #     config["y(km)"] = ((np.array(self.lat_lims) - config["center"][1]) * config["degree2km"]).tolist()
    #     config["z(km)"] = np.array(self.z_lims).tolist()
    #     config["bfgs_bounds"] = [list(config[x]) for x in config["dims"]] + [[None, None]]
    #     return config

    # @property
    # def stations(self):
    #     stations = ut.get_stations_GaMMA_df(self.response)
    #     config = self._get_config()
    #     stations["x(km)"] = stations["longitude"].apply(lambda x: (x - config["center"][0]) * config["degree2km"])
    #     stations["y(km)"] = stations["latitude"].apply(lambda x: (x - config["center"][1]) * config["degree2km"])
    #     stations["z(km)"] = stations["elevation(m)"].apply(lambda x: -x / 1e3)
    #     return stations
        

    def _get_config(self):

        config = self.__dict__

        in_proj = pyproj.Proj("EPSG:4326")
        out_proj = pyproj.Proj(self.epsg_proj)
        y_min,x_min = pyproj.transform(in_proj,out_proj,self.lat_lims[0],self.lon_lims[0])
        y_max,x_max = pyproj.transform(in_proj,out_proj,self.lat_lims[1],self.lon_lims[1])


        config["x(km)"] = np.array([x_min,x_max])/1e3
        config["y(km)"] = np.array([y_min,y_max])/1e3
        config["z(km)"] = np.array(self.z_lims)
        config["bfgs_bounds"] = (
                                (config["x(km)"][0] - 1, config["x(km)"][1] + 1),  # x
                                (config["y(km)"][0] - 1, config["y(km)"][1] + 1),  # y
                                (0, config["z(km)"][1] + 1),  # x
                                (None, None),  # t
                            )
        return config

    @property
    def stations(self):
        stations = ut.get_stations_GaMMA_df(self.response)

        in_proj = pyproj.Proj("EPSG:4326")
        out_proj = pyproj.Proj(self.epsg_proj)

        stations["x(km)"] = stations.apply(lambda x: pyproj.transform(in_proj,out_proj,
                                            x["latitude"], x["longitude"])[0] / 1e3, axis=1)
        stations["y(km)"] = stations.apply(lambda x: pyproj.transform(in_proj,out_proj,
                                            x["latitude"], x["longitude"])[1] / 1e3, axis=1)
        stations["z(km)"] = stations["elevation(m)"] / -1e3
        return stations

class GaMMA():
    def __init__(self,picks_csv,xml_path,out_dir):
        self.picks_csv = picks_csv
        self.xml_path = xml_path
        self.response = read_inventory(xml_path)
        self.out_dir = out_dir


    def associator(self,gamma_obj):
        picks_df = ut.get_picks_GaMMa_df(self.picks_csv,
                                self.response,
                                compute_amplitudes=gamma_obj.calculate_amp,
                                p_window=gamma_obj.p_window,
                                s_window=gamma_obj.s_window,
                                waterlevel = gamma_obj.waterlevel)
        stations = list(set(picks_df["station"].to_list()))

        station_df = gamma_obj.stations
        station_df =  station_df[station_df["station_name"].isin(stations)]
        station_df = station_df.reset_index(drop=True)

        config = gamma_obj.__dict__

        # pbar = tqdm(1)
        catalogs, assignments = association(picks_df, station_df, config, 
                                method=gamma_obj.method)
        
        catalog = pd.DataFrame(catalogs)
        assignments = pd.DataFrame(assignments,
                                columns=["pick_idx", "event_idx", "prob_gamma"])
        print(picks_df)
        print(catalog)
        print(assignments)
        return catalog,assignments,station_df

