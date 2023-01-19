from locale import D_FMT
import os
import datetime as dt
import pandas as pd
import obspy
from obspy.core.inventory.inventory import Inventory
from obspy.io.xseed.parser import Parser
from SeisMonitor.utils import printlog,isfile
from obspy.core.utcdatetime import UTCDateTime
import warnings
warnings.filterwarnings('ignore')
### EQTransformer util associator

SeisMonitor_columns = ["pick_id","arrival_time","probability","phasehint",
                    "network","station","location","instrument_type","author",
                    "creation_time","event_start_time","event_end_time",
                    "detection_probability","snr","station_lat","station_lon",
                    "station_elv","file_name"]

EQT_COLUMNS = ["file_name","network",'station',
    'instrument_type','station_lat','station_lon','station_elv',
    'event_start_time','event_end_time','detection_probability',
    'detection_uncertainty','p_arrival_time','p_probability',
    'p_uncertainty','p_snr','s_arrival_time',
    's_probability','s_uncertainty','s_snr']

GaMMA_picks_columns = ["id","timestamp","type","prob","amp"]

paz_wa = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}

def get_picks_GaMMa_df(picks,response,compute_amplitudes=True
                        ,p_window=10,
                         s_window=5,waterlevel=10):
    df = pd.read_csv(picks, dtype={'location': str})
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time"]
    # df["location"] = df["location"].apply(lambda x: "{:02d}".format(x))
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    if compute_amplitudes:
        df = get_seismonitor_amplitudes(df,response,
                                        p_window=p_window,
                                        s_window=s_window,
                                        waterlevel = waterlevel)
    else:
        df = df

    id_func = lambda x: x.split("-")[-1]
    df["id"] = df["pick_id"].apply(id_func)
    df = df.rename(columns={"arrival_time":"timestamp",
                        "probability":"prob",
                        "phasehint":"type",
                        "amplitude":"amp"})
    df = df.drop_duplicates(ignore_index=True)
    # print(df)
    return df

def get_stations_GaMMA_df(response):
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
            # components = []
            # sensitivities = []
            # for channel in station:
            #     component = channel.code[-1]
            #     sensitivity = channel.response.get_paz().normalization_factor
                
            #     components.append(component)
            #     sensitivities.append(sensitivity)


            # unit = channel.response.response_stages[0].input_units
            inst = station[0].code[:-1]
            loc = station[0].location_code

            scr_id = ".".join((net,sta,loc,inst))

            id_list.append(scr_id)
            station_list.append(sta)
            longitude_list.append(lon)
            latitude_list.append(lat)
            elevation_list.append(elv)

    df = {"id":id_list,
        "longitude":longitude_list,
        "latitude":latitude_list,
        "elevation(m)":elevation_list,
        "station_name":station_list,
        }
    df = pd.DataFrame(df)
    df = df.drop_duplicates(subset="id",ignore_index=True)
    # print(df)
    return df

def get_paz_from_response(seed_id,response,
                            datetime=None):
    if isinstance(response,Parser):
        try:
            paz = response.get_paz(seed_id,datetime)
        except:
            return None

    elif isinstance(response,Inventory):
        try:
            response = response.get_response(seed_id,
                                    datetime)
        
            paz_stage = response.get_paz()
        except:
            return None

        sensitivity = response.instrument_sensitivity.value
        gain = paz_stage.normalization_factor
        poles = paz_stage.poles
        zeros = paz_stage.zeros

        paz = {'poles': poles,
                'zeros': zeros,
                'gain': gain,
                'sensitivity': sensitivity}
    return paz

def get_amplitudes_from_pick(st,picktime,phasehint,
                            p_window=10,s_window=5):

    _st = st.copy()

    if phasehint.upper() == "P":
        trimmedtime = dt.timedelta(seconds=p_window)
        _st.trim(UTCDateTime(picktime),UTCDateTime(picktime)+trimmedtime)

        if len(_st)== 0:
            tr_z = _st.select(component="Z")[0]
            ampl = max(abs(tr_z.data))
        else:
            tr = st[0]
            ampl = max(abs(tr.data))

    elif phasehint.upper() == "S":
        trimmedtime = dt.timedelta(seconds=s_window)
        _st.trim(UTCDateTime(picktime),UTCDateTime(picktime)+trimmedtime)

        if len(_st)== 0:
            tr_n = _st.select(component="N")[0]
            ampl_n = max(abs(tr_n.data))
            tr_e = _st.select(component="E")[0]
            ampl_e = max(abs(tr_e.data))
            ampl = max(ampl_n, ampl_e)
        else:
            tr = st[0]
            ampl = max(abs(tr.data))

    return ampl

def get_amplitudes_from_local_st(file_name,df,response,
                                p_window=10,s_window=5,
                                waterlevel=10):
    st = obspy.read(file_name)

    for tr in st:
        paz = get_paz_from_response(tr.id,response,
                                tr.stats.starttime)
        if paz == None:
            print(f"\t->No response found: {tr.id}-{tr.stats.starttime}")

        tr.simulate(paz_remove=paz, 
                    paz_simulate=paz_wa, 
                    water_level=waterlevel)


    ampl_func = lambda x: get_amplitudes_from_pick(st,x.arrival_time,
                                            x.phasehint,p_window,s_window)
    df["amplitude"] = df.apply(ampl_func,axis=1)

    return df
    
def get_seismonitor_amplitudes(df,response,out=None,
                         p_window=10,s_window=5,
                        waterlevel=10):
    # df = pd.read_csv(picks)
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    gdf = df.groupby(by="file_name")
    dfs = []
    for file_name,df in gdf.__iter__():
        df = get_amplitudes_from_local_st(file_name,df,response,
                        p_window,s_window,waterlevel)
        dfs.append(df)

    df = pd.concat(dfs)

    fourth_column = df.pop('amplitude')
    df.insert(3, 'amplitude', fourth_column)

    df = df.sort_values(by="arrival_time",ignore_index=True)
    
    if out != None:
        isfile(out)
        df.to_csv(out,index=False)

    return df

def link(p_row,s_df,tt,st=2,et=4):
    p_time = p_row.arrival_time
    s_minus_p = (s_df["arrival_time"]-p_time).dt.total_seconds()
    s_df["tt"] = s_df["arrival_time"][abs(s_minus_p) <= tt]
    s_df = s_df.dropna(subset=['tt'])

    s_df = s_df.sort_values(by=["tt"],ignore_index=True)
    # print(s_df["tt"])
    if len(s_df) >= 1:
        p_row["s_arrival_time"] = s_df["tt"][0]
        p_row["s_probability"] = s_df["probability"][0]
        p_row["s_pick_id"] = s_df["pick_id"][0]

        p_row['event_start_time'] = p_time-dt.timedelta(seconds=st)
        p_row['event_end_time'] = p_time+dt.timedelta(seconds=et)
        p_row['detection_probability'] = (p_row["probability"] + p_row["s_probability"])/2
        
    else:
        p_row["s_arrival_time"] = None
        p_row["s_probability"] = None
        p_row["s_pick_id"] = None
    
        p_row['event_start_time'] = p_time-dt.timedelta(seconds=st)
        p_row['event_end_time'] = p_time+dt.timedelta(seconds=et*1)
        p_row['detection_probability'] = p_row["probability"]/2
    p_row['file_name'] = p_row.mseed_name
    p_row['station_lat'] = p_row.station_lat
    p_row['station_lon'] = p_row.station_lon
    p_row['station_elv'] = p_row.station_elv
    return p_row

def make_eqt_dirs(eqt_df,folder):
    dfbystation = eqt_df.groupby(by=["station"])

    for name,df in dfbystation.__iter__():
        dirname = name + "_outputs"
        dirpath = os.path.join(folder,dirname)

        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        
        filepath = os.path.join(dirpath,"X_prediction_results.csv")
        df.to_csv(filepath,index=False)

def link_seismonitor_phases(df,tt=30,st=2,et=4):
    """
    For each P pick, it searchs the corresponding S pick
    according tt parameter.
    Parameters:
    -----------
    df: pd.DataFrame
        seismonitor dataframe
    tt: float
        Each P pick search S picks in a tt-second radius.
        If there are more than one S pick, it chooses the closest.
    st: float
        Seconds before the P pick. It is used to set up
        the event start time
    et: float
        Seconds after the S pick. It is used to set up
        the event end time. If there isn't P pick, the event
        end time is configured et*2-seconds adfter the P phase time.
    Returns:
    --------
    df: pd.DataFrame
        P and S phases are linked by
    """
    columns = ["pick_id","arrival_time","probability",'phasehint','network',
             'station', 'location','instrument_type','author', 
             "station_lat","station_lon","station_elv",
            'creation_time',"mseed_name"]

    df = df [columns]
    df["arrival_time"] = pd.to_datetime(df["arrival_time"])
    dfbystation = df.groupby(by=["station"])

    eqt_df = pd.DataFrame()
    for i,df in dfbystation.__iter__():
    
        p_df = df[df["phasehint"] == "P"]
        s_df = df[df["phasehint"] == "S"]
        
        for j,row in p_df.iterrows():
            p_row = link(row,s_df,tt,st,et)
            eqt_df = eqt_df.append(p_row,ignore_index=True)

    eqt_df["p_snr"] = None
    eqt_df["p_uncertainty"] = None
    eqt_df["s_snr"] = None
    eqt_df["s_uncertainty"] = None
    eqt_df["detection_uncertainty"] = None

    eqt_df = eqt_df.rename(columns={"arrival_time":"p_arrival_time",
                            "probability":"p_probability"
                            })
    eqt_df = eqt_df[EQT_COLUMNS]
    return eqt_df

def link_eqt_phases(df):

    eqt_and_seismonitor_cols = ["file_name","network","station","instrument_type",
                "station_lat","station_lon","station_elv","event_start_time",
                "event_end_time","detection_probability"]

    p_df = df[df["phasehint"] == "P"]
    s_df = df[df["phasehint"] == "S"]

    df = p_df.merge(s_df,how="outer",on=eqt_and_seismonitor_cols,suffixes=("_p","_s"))

    df["p_uncertainty"] = None
    df["s_uncertainty"] = None
    df["detection_uncertainty"] = None

    renaming = {'pick_id_p':"p_pick_id", 'arrival_time_p':"p_arrival_time",
                 'probability_p':"p_probability", 'phasehint_p':"p_phasehint",
                 'location_p':"p_location", 'author_p':"p_author",
                'creation_time_p':"p_creation_time", 
                 'snr_p':"p_snr", 'pick_id_s':"s_pick_id", 
                 'arrival_time_s':"s_arrival_time",'probability_s':"s_probability",
                  'phasehint_s':"s_phasehint", 'location_s':"s_location",
                   'author_s':"s_location",'creation_time_s':"s_creation_time",
                'snr_s':"s_snr"}
    df = df.rename(columns=renaming)
    df = df[EQT_COLUMNS]
    return df

def seismonitor_picks_to_eqt_fmt(seismonitor_picks,eqt_folder,
                            tt=30,st=2,et=4):

    df = pd.read_csv(seismonitor_picks)
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    gdf = df.groupby(by="author")
    eqts = []
    for author,df in gdf.__iter__():
        if author == "EQTransformer":
            eqt = link_eqt_phases(df)
        else:
            eqt = link_seismonitor_phases(df,tt,st,et)
        eqts.append(eqt)
    eqt = pd.concat(eqts)

    make_eqt_dirs(eqt,eqt_folder)
    printlog('info','pnet_to_eqt_fmt',
                f'ok')