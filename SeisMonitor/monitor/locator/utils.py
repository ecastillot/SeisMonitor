from typing import Union
import pandas as pd
import datetime as dt
from itertools import groupby
import json
import SeisMonitor.utils as sut
from obspy import read_inventory
from obspy.core.event.catalog import Catalog, read_events
from obspy.core.inventory.inventory import Inventory

def changing_picks_info(ev, ref_picks):
    """Update event picks with information from reference picks.

    Args:
        ev (obspy.core.event.Event): Event object to modify
        ref_picks (dict): Dictionary of reference picks with station_phase_time keys

    Returns:
        obspy.core.event.Event: Modified event with updated picks and arrivals
    """
    true_picks = []
    picks_dict = {}
    
    for pick in ev.picks:
        station = pick.waveform_id.station_code
        phasehint = pick.phase_hint
        time = pick.time.strftime("%Y%m%dT%H%M%S")
        key = f"{station}_{phasehint}_{time}"
        
        # Try exact time match first, then ±1 second
        for offset in [0, -1, 1]:
            try:
                adj_time = (pick.time + dt.timedelta(seconds=offset)).strftime("%Y%m%dT%H%M%S")
                true_pick = ref_picks[f"{station}_{phasehint}_{adj_time}"]
                picks_dict[pick.resource_id.id] = true_pick
                true_picks.append(true_pick)
                break
            except KeyError:
                if offset == 1:  # Log warning only after all attempts fail
                    sut.printlog(
                        "warning", "NLLOC",
                        f"No information could be extracted for pick: {key}"
                    )

    ori_pref = ev.preferred_origin()
    true_arrivals = []
    
    for arrival in ori_pref.arrivals:
        try:
            arrival_prob = json.loads(picks_dict[arrival.pick_id.id].comments[0].text)
            arrival.time_weight = arrival_prob["probability"]
            arrival.comments = picks_dict[arrival.pick_id.id].comments
            arrival.pick_id.id = picks_dict[arrival.pick_id.id].resource_id.id
            true_arrivals.append(arrival)
        except (KeyError, json.JSONDecodeError):
            print("Picks don't have probability information")
        
    ori_pref.arrivals = true_arrivals
    ev.picks = true_picks
    return ev


def get_picks(catalog):
    """Extract picks from a catalog into a dictionary.

    Args:
        catalog (obspy.core.event.Catalog): Catalog containing events with picks

    Returns:
        dict: Dictionary of picks keyed by station_phase_time strings
    """
    picks = {}
    for ev in catalog:
        for pick in ev.picks:
            station = pick.waveform_id.station_code
            phasehint = pick.phase_hint
            time = pick.time.strftime("%Y%m%dT%H%M%S")
            picks[f"{station}_{phasehint}_{time}"] = pick
    return picks


def get_bad_and_good_events(reloc_catalog):
    """Separate events based on arrival azimuth uniformity.

    Args:
        reloc_catalog (obspy.core.event.Catalog): Relocated event catalog

    Returns:
        tuple: (good_events, bad_events) where good_events have varying azimuths
               and bad_events have uniform azimuths
    """
    def all_equal(iterable):
        """Check if all elements in an iterable are equal."""
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

    bad_events = []
    good_events = []
    
    for ev in reloc_catalog:
        pref_origin = ev.preferred_origin()
        azimuths = [x.azimuth for x in pref_origin.arrivals]
        if all_equal(azimuths):
            bad_events.append(ev)
        else:
            good_events.append(ev)
    
    return good_events, bad_events


def filter_arrivals_by_distance(events, distance, min_P_phases=3, min_S_phases=2):
    """Filter events based on arrival distance and minimum phase counts.

    Args:
        events (list): List of event objects
        distance (float): Maximum distance for arrivals in degrees
        min_P_phases (int): Minimum number of P phases required
        min_S_phases (int): Minimum number of S phases required

    Returns:
        list: Filtered list of events meeting criteria
    """
    new_events = []
    
    for ev in events:
        pref_origin = ev.preferred_origin()
        arrivals = pref_origin.arrivals
        picks = {pick.resource_id.id: pick for pick in ev.picks}
        
        new_arrivals = []
        counts = {"P": [], "S": []}
        
        for arrival in arrivals:
            if arrival.distance <= distance:
                new_arrivals.append(arrival)
            else:
                picks.pop(arrival.pick_id.id, None)
            counts[arrival.phase].append(arrival.phase)
        
        counts["P"] = len(counts["P"])
        counts["S"] = len(counts["S"])
        
        if counts["P"] < min_P_phases or counts["S"] < min_S_phases:
            continue
        
        for i, origin in enumerate(ev.origins):
            if origin.resource_id.id == ev.preferred_origin_id:
                ev.origins[i].arrivals = new_arrivals
        
        ev.picks = list(picks.values())
        new_events.append(ev)
    
    return new_events


def resp2df(resp):
    """Convert station response inventory to DataFrame.

    Args:
        resp (Union[str, Inventory]): Path to RESP file or Inventory object

    Returns:
        pandas.DataFrame: DataFrame with columns:
                          network, station, latitude, longitude, elevation
    """
    networks = []
    stations = []
    longitudes = []
    latitudes = []
    elevations = []
    
    inv = read_inventory(resp) if isinstance(resp, str) else resp
    
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
            elevations.append(sta.elevation)
            stations.append(sta.code)
            networks.append(net.code)
    
    df = pd.DataFrame({
        "network": networks,
        "station": stations,
        "latitude": latitudes,
        "longitude": longitudes,
        "elevation": elevations
    })
    return df


class VelModel:
    """Class representing a velocity model."""

    def __init__(self, vel_path, model_name=None, vp_vs_ratio=1.78, compute_vs=True):
        """Initialize velocity model.

        Args:
            vel_path (Union[str, pandas.DataFrame]): Path to velocity CSV or DataFrame
            model_name (str, optional): Name of the velocity model
            vp_vs_ratio (float): Vp/Vs ratio for computing Vs if not provided
            compute_vs (bool): Whether to compute Vs from Vp
        """
        self.model_name = model_name
        self.vel_path = vel_path
        self.vp_vs_ratio = vp_vs_ratio
        self.compute_vs = compute_vs
        
        self.vel = pd.read_csv(vel_path) if isinstance(vel_path, str) else vel_path

    def to_nlloc(self, out):
        """Write velocity model in NonLinLoc format.

        Args:
            out (str): Output file path
        """
        with open(out, 'w') as vm:
            vm.write("# model layers (LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad)\n")
            for i, row in self.vel.iterrows():
                vs = row.vp / self.vp_vs_ratio if self.compute_vs else row.vs
                enter = "" if i == len(self.vel) - 1 else "\n"
                vm.write(
                    f"LAYER    {row.depth:<6.2f}    {row.vp:<.2f}    0.00    "
                    f"{vs:<.2f}    0.00    {row.rho:<.2f}    0.00{enter}"
                )


class Stations:
    """Class representing seismic stations."""

    def __init__(self, stations):
        """Initialize stations from DataFrame or inventory.

        Args:
            stations (Union[pandas.DataFrame, str, Inventory]): Station data source

        Raises:
            Exception: If stations parameter is invalid (commented out in original)
        """
        if isinstance(stations, pd.DataFrame):
            self.stations = stations
        else:
            self.stations = resp2df(stations)

    def to_nlloc(self, out):
        """Write stations in NonLinLoc format.

        Args:
            out (str): Output file path
        """
        with open(out, 'w') as vs:
            vs.write("# GTSRCE label LATLON latSrce longSrce zSrce elev\n")
            for i, row in self.stations.iterrows():
                elv = row.elevation / 1e3  # Convert to kilometers
                enter = "" if i == len(self.stations) - 1 else "\n"
                vs.write(
                    f"GTSRCE  {row.station:<5}  LATLON  {row.latitude:<7.3f}  "
                    f"{row.longitude:<8.3f}  0.000  {elv:<7.3f}{enter}"
                )


class LocatorBasicInputs:
    """Class holding basic inputs for seismic location."""

    def __init__(self, vel_model: VelModel, stations: Stations):
        """Initialize locator inputs.

        Args:
            vel_model (VelModel): Velocity model object
            stations (Stations): Stations object
        """
        self.vel_model = vel_model
        self.stations = stations