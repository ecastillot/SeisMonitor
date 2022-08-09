from typing import Union
import pandas as pd
import SeisMonitor.utils as sut
from obspy import read_inventory
from obspy.core.event.catalog import Catalog, read_events

def resp2df(resp):
    """
    Parameters:
    -----------
    resp: str
        RESP filepath

    Returns: DataFrame
        Dataframe with the next columns
        network,station,latitude,longitude,elevation
    """
    networks = []
    stations = []
    longitudes = []
    latitudes = []
    elevations = []
    inv = read_inventory(resp)
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
            elevations.append(sta.elevation)
            stations.append(sta.code)
            networks.append(net.code)

    df = {"network":networks,"station":stations,
        "latitude":latitudes,"longitude":longitudes,
        "elevation":elevations}
    df = pd.DataFrame(df)
    return df 

class VelModel():
    def __init__(self,vel_path,
                vp_vs_ratio=1.78,
                compute_vs=True):

        self.vel_path = vel_path
        self.vp_vs_ratio = vp_vs_ratio
        self.compute_vs = compute_vs
        self.vel = pd.read_csv(vel_path)

    def to_nlloc(self,out):
        vm = open(out, 'w')
        msg ="# model layers (LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad)\n"
        vm.write(msg)
        for i,row in self.vel.iterrows():
            if self.compute_vs:
                vs = row.vp/self.vp_vs_ratio
            else:
                vs = row.vs

            if i == len(self.vel)-1:
                enter = ""
            else:
                enter = "\n"

            msg = f"LAYER    {row.depth:<6.2f}    {row.vp:<.2f}    0.00    {vs:<.2f}    0.00    {row.rho:<.2f}    0.00{enter}"
            vm.write(msg)
        vm.close()

class Stations():
    def __init__(self,stations_path) -> None:
        self.stations_path = stations_path
        self.stations = resp2df(stations_path)

    def to_nlloc(self,out):
        vs = open(out, 'w')
        msg = "# GTSRCE label LATLON latSrce longSrce zSrce elev\n"
        vs.write(msg)
        for i,row in self.stations.iterrows():
            elv = row.elevation/1e3
            if i == len(self.stations)-1:
                enter = ""
            else:
                enter = "\n"

            msg = f"GTSRCE  {row.station:<5}  LATLON  {row.latitude:<7.3f}  {row.longitude:<8.3f}  0.000  {elv:<7.3f}{enter}"
            vs.write(msg)
        vs.close()
        
class LocatorBasicInputs():
    def __init__(self,
                vel_model:VelModel,
                stations:Stations):
        self.vel_model = vel_model
        self.stations = stations

        
