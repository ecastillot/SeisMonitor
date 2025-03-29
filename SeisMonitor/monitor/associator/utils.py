import os
import datetime as dt
import pandas as pd
import obspy
from obspy.core.inventory.inventory import Inventory
from obspy.io.xseed.parser import Parser
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.utils import printlog, isfile
import warnings
import math

warnings.filterwarnings('ignore')

# Constants
SEISMONITOR_COLUMNS = [
    "pick_id", "arrival_time", "probability", "phasehint",
    "network", "station", "location", "instrument_type", "author",
    "creation_time", "event_start_time", "event_end_time",
    "detection_probability", "snr", "station_lat", "station_lon",
    "station_elv", "file_name"
]

EQT_COLUMNS = [
    "file_name", "network", "station", "instrument_type",
    "station_lat", "station_lon", "station_elv",
    "event_start_time", "event_end_time", "detection_probability",
    "detection_uncertainty", "p_arrival_time", "p_probability",
    "p_uncertainty", "p_snr", "s_arrival_time",
    "s_probability", "s_uncertainty", "s_snr"
]

GAMMA_PICKS_COLUMNS = ["id", "timestamp", "type", "prob", "amp"]

PAZ_WA = {
    'sensitivity': 2800,
    'zeros': [0j],
    'gain': 1,
    'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]
}


def get_picks_gamma_df(picks, response, compute_amplitudes=True, p_window=10, s_window=5, waterlevel=10):
    """Convert picks to GaMMA-compatible DataFrame format.
    
    Args:
        picks (str): Path to picks CSV file
        response (Inventory or Parser): Station response information
        compute_amplitudes (bool): Whether to compute amplitudes
        p_window (int): P-wave window length in seconds
        s_window (int): S-wave window length in seconds
        waterlevel (int): Waterlevel for amplitude calculation
        
    Returns:
        pandas.DataFrame: Formatted picks DataFrame
    """
    df = pd.read_csv(picks, dtype={'location': str})
    date_cols = ["arrival_time", "creation_time", "event_start_time", "event_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    def convert_loc(loc):
        """Convert location code to two-digit string format."""
        if pd.isna(loc):
            return ''
        try:
            return "{:02d}".format(int(loc))
        except ValueError:
            return ''

    df['location'] = df['location'].apply(convert_loc)
    
    if compute_amplitudes:
        df = get_seismonitor_amplitudes(
            df, response,
            p_window=p_window,
            s_window=s_window,
            waterlevel=waterlevel
        )

    id_func = lambda x: ".".join(x.split("-")[-1].split(".")[0:2])
    df["id"] = df["pick_id"].apply(id_func)
    df = df.rename(columns={
        "arrival_time": "timestamp",
        "probability": "prob",
        "phasehint": "type",
        "amplitude": "amp"
    })
    df = df.drop_duplicates(ignore_index=True)
    return df


def get_stations_gamma_df(response):
    """Convert station response to GaMMA-compatible DataFrame format.
    
    Args:
        response (Inventory): Station response inventory
        
    Returns:
        pandas.DataFrame: Formatted stations DataFrame
    """
    id_list = []
    station_list = []
    longitude_list = []
    latitude_list = []
    elevation_list = []

    for network in response:
        for station in network:
            net = network.code
            sta = station.code
            elv = station.elevation
            lat = station.latitude
            lon = station.longitude
            inst = station[0].code[:-1]
            loc = station[0].location_code
            scr_id = ".".join((net, sta))

            id_list.append(scr_id)
            station_list.append(sta)
            longitude_list.append(lon)
            latitude_list.append(lat)
            elevation_list.append(elv)

    df = pd.DataFrame({
        "id": id_list,
        "longitude": longitude_list,
        "latitude": latitude_list,
        "elevation(m)": elevation_list,
        "station_name": station_list
    })
    df = df.drop_duplicates(subset="id", ignore_index=True)
    return df


def get_paz_from_response(seed_id, response, datetime=None):
    """Extract Poles and Zeros (PAZ) information from response.
    
    Args:
        seed_id (str): SEED identifier
        response (Inventory or Parser): Response information
        datetime (UTCDateTime, optional): Specific time for response
        
    Returns:
        dict: PAZ dictionary or None if not found
    """
    if isinstance(response, Parser):
        try:
            paz = response.get_paz(seed_id, datetime)
        except:
            return None
    elif isinstance(response, Inventory):
        try:
            response = response.get_response(seed_id, datetime)
            paz_stage = response.get_paz()
            paz = {
                'poles': paz_stage.poles,
                'zeros': paz_stage.zeros,
                'gain': paz_stage.normalization_factor,
                'sensitivity': response.instrument_sensitivity.value
            }
        except:
            return None
    return paz


def get_amplitudes_from_pick(st, picktime, phasehint, p_window=10, s_window=5):
    """Calculate amplitude from a seismic trace for a given pick.
    
    Args:
        st (Stream): Obspy Stream object
        picktime (datetime): Pick time
        phasehint (str): Phase type ('P' or 'S')
        p_window (int): P-wave window length in seconds
        s_window (int): S-wave window length in seconds
        
    Returns:
        float: Maximum amplitude value
    """
    _st = st.copy()
    if phasehint.upper() == "P":
        trimmedtime = dt.timedelta(seconds=p_window)
        _st.trim(UTCDateTime(picktime), UTCDateTime(picktime) + trimmedtime)
        tr = _st.select(component="Z")[0] if len(_st) > 0 else st[0]
        ampl = max(abs(tr.data))
    elif phasehint.upper() == "S":
        trimmedtime = dt.timedelta(seconds=s_window)
        _st.trim(UTCDateTime(picktime), UTCDateTime(picktime) + trimmedtime)
        if len(_st) > 0:
            tr_n = _st.select(component="N")[0]
            tr_e = _st.select(component="E")[0]
            ampl = max(abs(tr_n.data), abs(tr_e.data))
        else:
            ampl = max(abs(st[0].data))
    return ampl


def get_amplitudes_from_local_st(file_name, df, response, p_window=10, s_window=5, waterlevel=10):
    """Calculate amplitudes from local seismic data file.
    
    Args:
        file_name (str): Path to seismic data file
        df (pandas.DataFrame): Picks DataFrame
        response (Inventory or Parser): Response information
        p_window (int): P-wave window length in seconds
        s_window (int): S-wave window length in seconds
        waterlevel (int): Waterlevel for simulation
        
    Returns:
        pandas.DataFrame: DataFrame with added amplitudes
    """
    st = obspy.read(file_name)
    for tr in st:
        paz = get_paz_from_response(tr.id, response, tr.stats.starttime)
        if paz is None:
            print(f"\t->No response found: {tr.id}-{tr.stats.starttime}")
        tr.simulate(paz_remove=paz, paz_simulate=PAZ_WA, water_level=waterlevel)

    ampl_func = lambda x: get_amplitudes_from_pick(st, x.arrival_time, x.phasehint, p_window, s_window)
    df["amplitude"] = df.apply(ampl_func, axis=1)
    return df


def get_seismonitor_amplitudes(df, response, out=None, p_window=10, s_window=5, waterlevel=10):
    """Add amplitude calculations to SeisMonitor picks.
    
    Args:
        df (pandas.DataFrame): Input picks DataFrame
        response (Inventory or Parser): Response information
        out (str, optional): Output CSV file path
        p_window (int): P-wave window length in seconds
        s_window (int): S-wave window length in seconds
        waterlevel (int): Waterlevel for simulation
        
    Returns:
        pandas.DataFrame: DataFrame with amplitudes
    """
    date_cols = ["arrival_time", "creation_time", "event_start_time", "event_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    gdf = df.groupby(by="file_name")
    dfs = [get_amplitudes_from_local_st(file_name, df, response, p_window, s_window, waterlevel)
           for file_name, df in gdf.__iter__()]
    
    df = pd.concat(dfs)
    fourth_column = df.pop('amplitude')
    df.insert(3, 'amplitude', fourth_column)
    df = df.sort_values(by="arrival_time", ignore_index=True)
    
    if out:
        isfile(out)
        df.to_csv(out, index=False)
    
    return df


def link(p_row, s_df, tt, st=2, et=4):
    """Link P and S phases for a given pick.
    
    Args:
        p_row (pandas.Series): P pick row
        s_df (pandas.DataFrame): S picks DataFrame
        tt (float): Time threshold in seconds
        st (float): Start time offset in seconds
        et (float): End time offset in seconds
        
    Returns:
        pandas.Series: Updated pick row with linked phases
    """
    p_time = p_row.arrival_time
    s_minus_p = (s_df["arrival_time"] - p_time).dt.total_seconds()
    s_df["tt"] = s_df["arrival_time"][abs(s_minus_p) <= tt]
    s_df = s_df.dropna(subset=['tt']).sort_values(by=["tt"], ignore_index=True)

    if len(s_df) >= 1:
        p_row.update({
            "s_arrival_time": s_df["tt"][0],
            "s_probability": s_df["probability"][0],
            "s_pick_id": s_df["pick_id"][0],
            'event_start_time': p_time - dt.timedelta(seconds=st),
            'event_end_time': p_time + dt.timedelta(seconds=et),
            'detection_probability': (p_row["probability"] + s_df["probability"][0]) / 2
        })
    else:
        p_row.update({
            "s_arrival_time": None,
            "s_probability": None,
            "s_pick_id": None,
            'event_start_time': p_time - dt.timedelta(seconds=st),
            'event_end_time': p_time + dt.timedelta(seconds=et),
            'detection_probability': p_row["probability"] / 2
        })
    
    p_row.update({
        'file_name': p_row.mseed_name,
        'station_lat': p_row.station_lat,
        'station_lon': p_row.station_lon,
        'station_elv': p_row.station_elv
    })
    return p_row


def make_eqt_dirs(eqt_df, folder):
    """Create EQTransformer-style directory structure and files.
    
    Args:
        eqt_df (pandas.DataFrame): EQTransformer formatted DataFrame
        folder (str): Output directory path
    """
    dfbystation = eqt_df.groupby(by=["station"])
    for name, df in dfbystation.__iter__():
        dirname = name + "_outputs"
        dirpath = os.path.join(folder, dirname)
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        
        filepath = os.path.join(dirpath, "X_prediction_results.csv")
        df.to_csv(filepath, index=False)


def link_seismonitor_phases(df, tt=30, st=2, et=4):
    """Link P and S phases in SeisMonitor format.
    
    Args:
        df (pandas.DataFrame): SeisMonitor picks DataFrame
        tt (float): Time threshold in seconds for S pick search
        st (float): Start time offset in seconds
        et (float): End time offset in seconds
        
    Returns:
        pandas.DataFrame: EQTransformer formatted DataFrame
    """
    columns = [
        "pick_id", "arrival_time", "probability", 'phasehint', 'network',
        'station', 'location', 'instrument_type', 'author',
        "station_lat", "station_lon", "station_elv",
        'creation_time', "mseed_name"
    ]
    df = df[columns]
    df["arrival_time"] = pd.to_datetime(df["arrival_time"])
    dfbystation = df.groupby(by=["station"])

    eqt_df = pd.DataFrame()
    for _, df in dfbystation.__iter__():
        p_df = df[df["phasehint"] == "P"]
        s_df = df[df["phasehint"] == "S"]
        for _, row in p_df.iterrows():
            p_row = link(row, s_df, tt, st, et)
            eqt_df = eqt_df.append(p_row, ignore_index=True)

    eqt_df = eqt_df.assign(
        p_snr=None,
        p_uncertainty=None,
        s_snr=None,
        s_uncertainty=None,
        detection_uncertainty=None
    ).rename(columns={
        "arrival_time": "p_arrival_time",
        "probability": "p_probability"
    })[EQT_COLUMNS]
    
    return eqt_df


def link_eqt_phases(df):
    """Link P and S phases in EQTransformer format.
    
    Args:
        df (pandas.DataFrame): Input picks DataFrame
        
    Returns:
        pandas.DataFrame: EQTransformer formatted DataFrame
    """
    eqt_and_seismonitor_cols = [
        "file_name", "network", "station", "instrument_type",
        "station_lat", "station_lon", "station_elv", "event_start_time",
        "event_end_time", "detection_probability"
    ]
    p_df = df[df["phasehint"] == "P"]
    s_df = df[df["phasehint"] == "S"]
    df = p_df.merge(s_df, how="outer", on=eqt_and_seismonitor_cols, suffixes=("_p", "_s"))

    df = df.assign(
        p_uncertainty=None,
        s_uncertainty=None,
        detection_uncertainty=None
    ).rename(columns={
        'pick_id_p': "p_pick_id", 'arrival_time_p': "p_arrival_time",
        'probability_p': "p_probability", 'phasehint_p': "p_phasehint",
        'location_p': "p_location", 'author_p': "p_author",
        'creation_time_p': "p_creation_time", 'snr_p': "p_snr",
        'pick_id_s': "s_pick_id", 'arrival_time_s': "s_arrival_time",
        'probability_s': "s_probability", 'phasehint_s': "s_phasehint",
        'location_s': "s_location", 'author_s': "s_location",
        'creation_time_s': "s_creation_time", 'snr_s': "s_snr"
    })[EQT_COLUMNS]
    
    return df


def seismonitor_picks_to_eqt_fmt(seismonitor_picks, eqt_folder, tt=30, st=2, et=4):
    """Convert SeisMonitor picks to EQTransformer format.
    
    Args:
        seismonitor_picks (str): Path to SeisMonitor picks CSV
        eqt_folder (str): Output directory for EQTransformer files
        tt (float): Time threshold in seconds
        st (float): Start time offset in seconds
        et (float): End time offset in seconds
    """
    df = pd.read_csv(seismonitor_picks)
    date_cols = ["arrival_time", "creation_time", "event_start_time", "event_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    gdf = df.groupby(by="author")
    eqts = []
    for author, df in gdf.__iter__():
        eqt = link_eqt_phases(df) if author == "EQTransformer" else link_seismonitor_phases(df, tt, st, et)
        eqts.append(eqt)
    
    eqt = pd.concat(eqts)
    make_eqt_dirs(eqt, eqt_folder)
    printlog('info', 'pnet_to_eqt_fmt', 'ok')