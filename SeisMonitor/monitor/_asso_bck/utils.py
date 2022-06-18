import os
import datetime as dt
import pandas as pd
from SeisMonitor.utils import printlog,isfile
import warnings
warnings.filterwarnings('ignore')
### EQTransformer util associator

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
    eqt_columns = ["file_name","network",'station',
    'instrument_type','station_lat','station_lon','station_elv',
    'event_start_time','event_end_time','detection_probability',
    'detection_uncertainty','p_arrival_time','p_probability',
    'p_uncertainty','p_snr','s_arrival_time',
    's_probability','s_uncertainty','s_snr']

    eqt_df = eqt_df[eqt_columns]
    return eqt_df

def seismonitor_picks_to_eqt_fmt(seismonitor_picks,eqt_folder,
                            tt=30,st=2,et=4):

    df = pd.read_csv(seismonitor_picks)
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time",
                "mseed_start_time","mseed_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)

    eqt = link_seismonitor_phases(df,tt,st,et)
    printlog('info','link_phasenet_phases',
                f'ok')

    make_eqt_dirs(eqt,eqt_folder)
    printlog('info','pnet_to_eqt_fmt',
                f'ok')