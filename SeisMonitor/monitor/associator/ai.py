import numpy as np
import pandas as pd
import json
import math
import datetime as dt
import pyproj
import os
from . import utils as ut
from SeisMonitor.utils import isfile
from gamma.utils import association
from obspy.core.inventory.inventory import (Inventory,read_inventory)
from tqdm import tqdm
from obspy import UTCDateTime
from obspy.core.event.event import Event
from obspy.core.event import ResourceIdentifier,Catalog
from obspy.core.event.base import CreationInfo
from obspy.core.event.origin import Pick
from obspy.core.event.base import (QuantityError,
                                WaveformStreamID,
                                CreationInfo,
                                Comment)
from obspy.core.event import ResourceIdentifier
from obspy.core.event.origin import Origin, OriginQuality
from obspy.core.event.base import QuantityError
from obspy.core.event.origin import Arrival

class GaMMAObj():
    def __init__(self,region,epsg_proj,
                use_dbscan=True,use_amplitude=True,
                dbscan_eps=10.0,dbscan_min_samples=3,
                vel = {"p": 7.0, "s": 7.0 / 1.75},
                method="BGMM",oversample_factor=20,
                min_picks_per_eq=5,max_sigma11=2.0,
                max_sigma22=1.0,max_sigma12=1.0,
                calculate_amp=True,p_window=10,
                s_window=5,waterlevel=10
                ):
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
        self.response = None
        self.name = "GaMMA"

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
        
    def add_response(self,response):
        if isinstance(response,Inventory):
            self.response = response
        else:
            self.response = read_inventory(response)

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

        if self.response == None:
            raise Exception("You must add response file to the GaMMA Object")

        stations = ut.get_stations_GaMMA_df(self.response)

        in_proj = pyproj.Proj("EPSG:4326")
        out_proj = pyproj.Proj(self.epsg_proj)

        stations["y(km)"] = stations.apply(lambda x: pyproj.transform(in_proj,out_proj,
                                            x["latitude"], x["longitude"])[0] / 1e3, axis=1)
        stations["x(km)"] = stations.apply(lambda x: pyproj.transform(in_proj,out_proj,
                                            x["latitude"], x["longitude"])[1] / 1e3, axis=1)
        stations["z(km)"] = stations["elevation(m)"] / -1e3

        stations = stations.drop_duplicates(ignore_index=True)
        return stations

def get_gamma_picks(event_picks):
    pick_list = []
    for i, row in event_picks.iterrows():
        # print(row.location)
        # if math.isnan(row.location):
        #     loc = ""
        # else:
        #     loc = "{:02d}".format(int(row.location))
        loc = row.location
        str_id = ".".join((row.network,row.station,
                           loc,row.instrument_type + "Z"))
        
        comment = {'probability': row.prob,
                'GaMMA_probability':row.prob_gamma}

        if row.author == "EQTransformer":
            comment['snr'] = row.snr
            comment['detection_probability'] = row.detection_probability
            comment['event_start_time'] = row.event_start_time.strftime("%Y-%m-%d %H:%M:%S.%f")
            comment['event_end_time'] = row.event_end_time.strftime("%Y-%m-%d %H:%M:%S.%f")
       
        pick_obj = Pick(resource_id=ResourceIdentifier( id= row.pick_id,prefix="pick"),
                                time=UTCDateTime(row.timestamp),
                                time_errors=QuantityError(uncertainty=20/100, # 20 samples before and 20 samples after
                                                        confidence_level=row.prob*100),
                                waveform_id=WaveformStreamID(network_code=row.network,
                                                             station_code=row.station,
                                                             location_code=loc,
                                                             channel_code=row.instrument_type + "Z",
                                                             resource_uri= ResourceIdentifier(id= str_id ),
                                                             seed_string = str_id),
                                phase_hint = row["type"],
                                evaluation_mode = 'automatic',
                                creation_info= CreationInfo(author=row.author,
                                                            creation_time=UTCDateTime.now()),
                                method_id=row.author,
                                comments= [Comment(text=json.dumps(comment))]
                                )
        pick_list.append(pick_obj)
    
    return pick_list

def picks2arrivals(picks):
    arrivals = []
    for pick in picks:
        arrival = Arrival( resource_id =ResourceIdentifier(id=pick.resource_id.id,
                                                            prefix='arrival'),
                           pick_id = pick.resource_id,
                           phase = pick.phase_hint,
                           time_weight = pick.time_errors.confidence_level/100,
                           creation_info=pick.creation_info )
        arrivals.append(arrival)
    return arrivals
        
def get_gamma_origin(catalog_info,event_picks,in_proj="EPSG:3116",
                    out_proj ="EPSG:4326"):
    y,x = catalog_info["y(km)"]*1e3,catalog_info["x(km)"]*1e3
    in_proj = pyproj.Proj(in_proj) #colombia/MAGNA SIRGAS BOGOTÁ/ZONE
    out_proj = pyproj.Proj(out_proj)
    lat,lon = pyproj.transform(in_proj,out_proj,y,x)

    origin = Origin(resource_id = ResourceIdentifier(id=UTCDateTime(catalog_info.time).strftime("%Y%m%d.%H%M%S.%f"),
                                                    prefix='origin'),
                    time = UTCDateTime(catalog_info.time),
                    time_errors = QuantityError(uncertainty=catalog_info.sigma_time),
                    longitude = lon,
                    longitude_errors = QuantityError(),
                    latitude = lat,
                    latitude_errors = QuantityError(),
                    depth = catalog_info["z(km)"],
                    # depth = catalog_info["z(km)"]*1e3,
                    depth_errors = QuantityError(),
                    method_id = ResourceIdentifier(id="GaMMA"),
                    arrivals = picks2arrivals(event_picks),
                    quality = OriginQuality(associated_phase_count=len(event_picks)),
                    evaluation_status = "preliminary",
                    evaluation_mode = "automatic",
                    creation_info = CreationInfo(author="SeisMonitor",
                                                creation_time=UTCDateTime.now())
                     )
    return origin

def get_gamma_catalog(picks_df,catalog_df):

    events = []
    for i, row in catalog_df.iterrows():
        picks = picks_df[picks_df["event_idx"]==i]
        picks = get_gamma_picks(picks)
        origin = get_gamma_origin(row,picks)
        ev = Event(resource_id = ResourceIdentifier(id=origin.resource_id.id,
                                                            prefix='event'),
                    event_type = "earthquake",
                    event_type_certainty = "known",
                    picks = picks,
                    amplitudes = [],
                    focal_mechanisms = [],
                    origins = [origin],
                    magnitudes = [],
                    station_magnitudes = [],
                    creation_info = CreationInfo(author="SeisMonitor",
                                                creation_time=UTCDateTime.now())
                    )
        ev.preferred_origin_id = origin.resource_id.id
        events.append(ev)
    
    catalog = Catalog(events = events,
                    resource_id = ResourceIdentifier(prefix='catalog'),
                    creation_info = CreationInfo(author="SeisMonitor",
                                    creation_time=UTCDateTime.now()) )
    return catalog


class GaMMA():
    def __init__(self,gamma_obj):
        self.gamma_obj = gamma_obj

    def associate(self,picks_csv,xml_path,out_dir):
        self.picks_csv = picks_csv
        self.xml_path = xml_path
        self.response = read_inventory(xml_path)
        self.out_dir = out_dir
        self.xml_out_file = os.path.join(out_dir,"associations.xml")
        self.catalog_out_file = os.path.join(out_dir,"catalog.csv")
        self.picks_out_file = os.path.join(out_dir,"picks.csv")

        picks_df = ut.get_picks_GaMMa_df(self.picks_csv,
                                self.response,
                                compute_amplitudes=self.gamma_obj.calculate_amp,
                                p_window=self.gamma_obj.p_window,
                                s_window=self.gamma_obj.s_window,
                                waterlevel = self.gamma_obj.waterlevel)
        stations = list(set(picks_df["station"].to_list()))

        self.gamma_obj.add_response(self.response)
        station_df = self.gamma_obj.stations
        station_df =  station_df[station_df["station_name"].isin(stations)]
        station_df = station_df.reset_index(drop=True)

        config = self.gamma_obj.__dict__
        pbar = tqdm(1)
        meta = station_df.merge(picks_df["id"], how="right", on="id")

        # picks_df = picks_df.iloc[0:200]
        # print(picks_df) 
        # station_df.info())
        # picks_df.to_csv("./test.csv",index=False)
        catalogs, assignments = association(picks_df, station_df, config, 
                                method=self.gamma_obj.method,pbar=pbar)
        # print(catalogs, assignments,config)
        if not catalogs:
            return Catalog(),pd.DataFrame(),pd.DataFrame()
        
        catalog = pd.DataFrame(catalogs)
        
        catalog["time"] = pd.to_datetime(catalog["time"],format="%Y-%m-%dT%H:%M:%S.%f")
        catalog["event_index"] = catalog["time"].apply(lambda x: x)

        assignments = pd.DataFrame(assignments,
                                columns=["pick_index", "event_idx", "prob_gamma"])
        assignments = assignments.set_index("pick_index")
        

        picks = pd.merge(picks_df, assignments, left_index=True, 
                        right_index=True, how='right')
        isfile(self.catalog_out_file)
        catalog.to_csv(self.catalog_out_file,index=False)
        isfile(self.picks_out_file)
        picks.to_csv(self.picks_out_file,index=False)
       
        obspy_catalog = get_gamma_catalog(picks,catalog)
        isfile(self.xml_out_file)
        obspy_catalog.write(self.xml_out_file,format="SC3ML")


        return obspy_catalog,catalog,picks

if __name__ == "__main__":
    csv = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/associator/test.csv"
    cat_csv = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/associator/cat_test.csv"
    df = pd.read_csv(csv)
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    cat_df = pd.read_csv(cat_csv)
    cat_df["time"] = pd.to_datetime(cat_df["time"])

    # print(df)
    # print(cat_df)
    get_gamma_catalog(df,cat_df)