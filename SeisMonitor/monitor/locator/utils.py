from typing import Union
import pandas as pd
import datetime as dt
from itertools import groupby
import json
import SeisMonitor.utils as sut
from obspy import read_inventory
from obspy.core.event.catalog import Catalog, read_events
from obspy.core.inventory.inventory import (Inventory,read_inventory)

def changing_picks_info(ev,ref_picks):
    true_picks = []
    picks_dict = {}
    for pick in ev.picks:
        station = pick.waveform_id.station_code
        phasehint = pick.phase_hint
        try:
            time = pick.time.strftime("%Y%m%dT%H%M%S")
            true_pick = ref_picks[station+"_"+phasehint+"_"+time]
        except:
            try:
                time = (pick.time-dt.timedelta(seconds=1)).strftime("%Y%m%dT%H%M%S")
                true_pick = ref_picks[station+"_"+phasehint+"_"+time]
            except:
                try:
                    time = (pick.time+dt.timedelta(seconds=1)).strftime("%Y%m%dT%H%M%S")
                    true_pick = ref_picks[station+"_"+phasehint+"_"+time]
                except:
                    sut.printlog("warning","NLLOC",f'no information could be extracted from the mext pick in the input catalog: {station+"_"+phasehint+"_"+time}')

        # print(list(ref_picks.keys()))

        # picks_dict[pick.resource_id.id] = true_pick.resource_id.id
        picks_dict[pick.resource_id.id] = true_pick
        true_picks.append(true_pick)

    ori_pref  = ev.preferred_origin()
    true_arrivals = []
    for arrival in ori_pref.arrivals:
        try:
            arrival_prob = json.loads(picks_dict[arrival.pick_id.id].comments[0].text)
            arrival.time_weight = arrival_prob["probability"]
        except:
            print("picks don't have probability")
        arrival.comments = picks_dict[arrival.pick_id.id].comments
        # print(picks_dict[arrival.pick_id.id])
        arrival.pick_id.id = picks_dict[arrival.pick_id.id].resource_id.id
        true_arrivals.append(arrival)

    

    ori_pref.arrivals = true_arrivals
    ev.picks = true_picks
    # for arrival in ori_pref.arrivals:
    return ev

def get_picks(catalog):
    picks = {}
    for ev in catalog:
        for pick in ev.picks:
            station = pick.waveform_id.station_code
            phasehint = pick.phase_hint
            time = pick.time.strftime("%Y%m%dT%H%M%S")
            picks[station+"_"+phasehint+"_"+time] = pick
    return picks

def get_bad_and_good_events(reloc_catalog):

    def all_equal(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)
    bad_events = []
    good_events = []
    for ev in reloc_catalog:
        pref_origin = ev.preferred_origin()
        arrivals = pref_origin.arrivals
        azimuths = [x.azimuth for x in arrivals] 
        equal = all_equal(azimuths)
        if equal:
            bad_events.append(ev)
                
        else:
            good_events.append(ev)
    return good_events,bad_events

def filter_arrivals_by_distance(events,distance,
                                        min_P_phases=3,
                                        min_S_phases=2):
    new_events = []
    for ev in events:
        pref_origin = ev.preferred_origin()
        arrivals = pref_origin.arrivals

        picks = {}
        for pick in ev.picks:
            picks[pick.resource_id.id] = pick

        new_arrivals = []
        counts = {"P":[],"S":[]}
        for arrival in arrivals:
            if arrival.distance <= distance:
                new_arrivals.append(arrival)
            else:
                picks.pop(arrival.pick_id.id, None)

            counts[arrival.phase].append(arrival.phase)
        
        counts["P"] = len(counts["P"])
        counts["S"] = len(counts["S"])

        if (counts["P"] < min_P_phases) or\
            (counts["S"] < min_S_phases) :
            continue

        for i,origin in enumerate(ev.origins):
            if origin.resource_id.id == ev.preferred_origin_id:
                ev.origins[i].arrivals = new_arrivals 
            else:
                continue

        ev.picks = list(picks.values())
        new_events.append(ev)
    return new_events

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

    if isinstance(resp,str):
        inv = read_inventory(resp)
    else:
        inv = resp

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
                model_name=None,
                vp_vs_ratio=1.78,
                compute_vs=True):
        self.model_name = model_name
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
    def __init__(self,stations) -> None:
        self.stations = resp2df(stations)
        # else:
        #     raise Exception("Error: NLLOC Stations class")

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

        
