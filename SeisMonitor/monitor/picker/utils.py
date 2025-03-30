# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:49:23
#  * @modify date 2021-12-22 09:49:23
#  * @desc [description]
#  */
import os
import csv
import sys
import json
import math
import time
import obspy
import shutil
import logging
import datetime
import numpy as np
import pandas as pd
import datetime as dt
from obspy import UTCDateTime
from datetime import timedelta
from SeisMonitor.utils import printlog, isfile
from SeisMonitor.monitor.downloader.utils import get_chunktimes
import concurrent.futures as cf
from git import Repo
from obspy.core.event.origin import Pick
from obspy.core.event.base import QuantityError, WaveformStreamID, CreationInfo, Comment
from obspy.core.event import ResourceIdentifier

# Define standard columns for SeisMonitor
SeisMonitor_columns = [
    "pick_id", "arrival_time", "probability", "phasehint",
    "network", "station", "location", "instrument_type", "author",
    "creation_time", "event_start_time", "event_end_time",
    "detection_probability", "snr", "station_lat", "station_lon",
    "station_elv", "file_name"
]


def clone_aipicker(name, output_folder):
    """Clone AI picker repository from GitHub.
    
    Args:
        name (str): Name of the AI picker ('EQTransformer' or 'PhaseNet')
        output_folder (str): Directory path to store the cloned repository
        
    Returns:
        bool: True if successful
        
    Raises:
        Exception: If invalid picker name is provided
    """
    git_urls = {
        "PhaseNet": "https://github.com/ecastillot/PhaseNet.git",
        "EQTransformer": "https://github.com/ecastillot/EQTransformer.git"
    }
    
    if name not in git_urls:
        raise Exception(
            f"Invalid picker name: {name}. Choose 'EQTransformer' or 'PhaseNet'"
        )
    
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    
    Repo.clone_from(git_urls[name], output_folder)
    return True


def eqt_picks_2_seismonitor_fmt(eqt_folder, mseed_folder, out_path):
    """Convert EQTransformer picks to SeisMonitor format.
    
    Args:
        eqt_folder (str): Directory containing EQTransformer pick files
        mseed_folder (str): Directory containing MSEED files
        out_path (str): Output path for the converted CSV file
        
    Returns:
        pandas.DataFrame: Converted picks DataFrame
    """
    date_cols = ["p_arrival_time", "s_arrival_time", "event_start_time", "event_end_time"]
    
    def parse_datetime(datetime_str):
        """Parse datetime string with multiple possible formats."""
        if isinstance(datetime_str, float):
            return None
        for fmt in ["%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"]:
            try:
                return dt.datetime.strptime(datetime_str, fmt)
            except ValueError:
                pass
        return None

    dfs = []
    for dp, dn, filenames in os.walk(eqt_folder):
        for f in filenames:
            if f == "X_prediction_results.csv" :
                search_path = os.path.join(dp, f)
                df = pd.read_csv(search_path,dtype={'station': str,'location':str})
                if not df.empty:
                    # to_date = lambda x: pd.to_datetime(x,format="mixed")
                    # df[date_cols] = df[date_cols].apply(to_date)

                    to_date = lambda x: parse_datetime(x)
                    # df[date_cols] = df[date_cols].apply(to_date)

                    for date_col in date_cols:
                        df[date_col] = df[date_col].apply(parse_datetime)
                    dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs,ignore_index=True)
    df["station"] = df["station"].apply(lambda x: x.strip())
    df = df.sort_values(by="p_arrival_time",ignore_index=True)

    x_results_path = os.path.join(os.path.dirname(out_path),
                                "X_prediction_results_merge.csv")
    df.to_csv(x_results_path,index=False,
                date_format="%Y-%m-%d %H:%M:%S.%f")

    #seismonitor
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time"]

    get_loc = lambda x: x.file_name.split(".")[2]
    df["location"] = df.apply(get_loc,axis=1)
    df["author"] = "EQTransformer"

    filename_func = lambda x: os.path.join(mseed_folder,x)
    df["file_name"] = df["file_name"].apply(filename_func)

    p_cols = df.columns.to_list()
    s_cols = df.columns.to_list()
    dfs = []
    for phase in ["p","s"]:
        if phase.lower() == "p":
            rm = "s"
        else:
            rm = "p"

        phase_cols = df.columns.to_list()
        not_rm_cols = [f'{phase}_arrival_time', f'{phase}_probability',
                         f'{phase}_uncertainty', f'{phase}_snr']
        rm_cols = [f'{rm}_arrival_time', f'{rm}_probability',
                 f'{rm}_uncertainty', f'{rm}_snr']
        for x in rm_cols:
            phase_cols.remove(x)

        phase_df = df[phase_cols]

        new_cols_gen = lambda x: x[2:]
        new_cols = list(map(new_cols_gen,not_rm_cols))

        columns = dict(zip(not_rm_cols, new_cols))
        phase_df = phase_df.rename(columns=columns)
        phase_df["phasehint"] = phase.upper()
        phase_df = phase_df.dropna(subset=["arrival_time"])

        phase_df["creation_time"] = dt.datetime.now()
        pick_id = lambda x: id_maker(x.arrival_time, x.network, 
                            x.station, x.location, 
                            x.instrument_type, x.phasehint)
        phase_df["pick_id"] = phase_df.apply(pick_id,axis=1) 


        phase_df[date_cols] = phase_df[date_cols].apply(pd.to_datetime)

        dfs.append(phase_df)
        

    df = pd.concat(dfs,ignore_index=True)
    df = df[SeisMonitor_columns]
    df = df.sort_values(by="arrival_time",ignore_index=True)
    df.to_csv(out_path,index=False,
            date_format="%Y-%m-%d %H:%M:%S.%f")
    return df


def get_filenames(mseed, filter_net=None, filter_sta=None, filter_cha=None):
    """Retrieve MSEED filenames with optional filtering.
    
    Args:
        mseed (str): Path to MSEED files directory
        filter_net (list, optional): Networks to exclude
        filter_sta (list, optional): Stations to exclude
        filter_cha (list, optional): Channels to include (if specified)
        
    Returns:
        list: Filtered filenames
        
    Raises:
        Exception: If filter parameters are not lists
    """
    filter_net = filter_net or []
    filter_sta = filter_sta or []
    filter_cha = filter_cha or []
    
    for param in (filter_net, filter_sta, filter_cha):
        if not isinstance(param, list):
            raise Exception(f"{param} must be a list")

    filenames = []
    file_list = [f for f in os.listdir(mseed) if f.lower().endswith('.mseed')]
    for filename in file_list:
        net, sta, loc, cha = filename.split('__')[0].split('.')
        if net in filter_net or sta in filter_sta or (filter_cha and cha not in filter_cha):
            continue
        filenames.append(filename)
    return filenames


def mv_mseed2onefolder(mseed_folder, one_folder):
    """Move MSEED files from subfolders to a single folder.
    
    Args:
        mseed_folder (str): Source directory with MSEED files
        one_folder (str): Destination directory
    """
    for dp, _, filenames in os.walk(mseed_folder):
        for f in filenames:
            if f.endswith(".mseed"):
                src = os.path.join(dp, f)
                dst = os.path.join(one_folder, f)
                printlog("debug", "mv_mseed2onefolder", f"{src} -> {dst}")
                shutil.move(src, dst)


def mv_mseed2stationfolder(one_folder, mseed_folder):
    """Move MSEED files to station-specific subfolders.
    
    Args:
        one_folder (str): Source directory with MSEED files
        mseed_folder (str): Destination base directory
    """
    for dp, _, filenames in os.walk(one_folder):
        for f in filenames:
            if f.endswith(".mseed"):
                src = os.path.join(dp, f)
                net, sta, _, _ = f.split("__")[0].split(".")
                dst = os.path.join(mseed_folder, sta, f)
                if not os.path.isdir(os.path.dirname(dst)):
                    os.makedirs(os.path.dirname(dst))
                printlog("debug", "mv_mseed2stationfolder", f"{src} -> {dst}")
                shutil.move(src, dst)


def make_dataframe(mseed, json_path, filter_net=None, filter_sta=None, filter_cha=None):
    """Create DataFrame from MSEED files and JSON metadata.
    
    Args:
        mseed (str): Directory containing MSEED files
        json_path (str): Path to JSON metadata file
        filter_net (list, optional): Networks to exclude
        filter_sta (list, optional): Stations to exclude
        filter_cha (list, optional): Channels to include
        
    Returns:
        pandas.DataFrame: DataFrame with waveform metadata
    """
    with open(json_path, "r") as jfile:
        datajson = json.load(jfile)

    filenames = get_filenames(mseed, filter_net, filter_sta, filter_cha)
    mseed_paths = []
    
    for onefile in filenames:
        df = {}
        path = os.path.join(mseed, onefile)
        st = obspy.read(path)
        stats = st[0].stats
        
        for key, val in st._get_common_channels_info().items():
            net, sta, loc, instrument_type = key
            channels = list(val["channels"].keys())

        df["fname"] = path
        for channel in channels:
            direction = channel[-1]
            if len(channels) == 1:
                df["E"] = df["N"] = df["Z"] = channel
            else:
                df[direction] = channel

        instrument = instrument_type[:-1]
        instruments = [x[:-1] for x in datajson[sta]['channels']]
        idx = next((i for i, x in enumerate(instruments) if x == instrument), None)
        if idx is None:
            continue
        
        sampling_rate = datajson[sta]['sampling_rate'][idx]
        lat, lon, elv = datajson[sta]['coords']
        
        df.update({
            'mseed_start_time': stats.starttime.strftime("%Y-%m-%d %H:%M:%S.%f"),
            'mseed_end_time': stats.endtime.strftime("%Y-%m-%d %H:%M:%S.%f"),
            "network": net,
            "station": sta,
            "location": loc,
            "instrument_type": instrument,
            "sampling_rate": sampling_rate,
            "sta_lat": lat,
            "sta_lon": lon,
            "sta_elv": elv
        })
        mseed_paths.append(df)
    
    return pd.DataFrame(mseed_paths)


def make_phasenet_datalist(
        datadir, json_path, datalist_path, channel_list=None,
        filter_network=None, filter_station=None, **kwargs
    ):
    """Create PhaseNet datalist CSV file.
    
    Args:
        datadir (str): Directory with MSEED files
        json_path (str): Path to JSON metadata
        datalist_path (str): Output path for datalist CSV
        channel_list (list, optional): Channels to include
        filter_network (list, optional): Networks to exclude
        filter_station (list, optional): Stations to exclude
        kwargs: Additional keyword arguments
        
    Returns:
        str: Directory containing the datalist
    """
    df = make_dataframe(datadir, json_path, filter_network, filter_station, channel_list)
    datalist_dir = os.path.dirname(datalist_path)
    isfile(datalist_path)
    df.to_csv(datalist_path, index=False)
    return datalist_dir


def create_datalist(self, all_in_folder=False):
    """Create PhaseNet datalist(s) for processing.
    
    Args:
        self: Object with datalist_dir and mseed_storage attributes
        all_in_folder (bool): If True, create single datalist; otherwise, per station
    """
    tic = time.time()
    logger = logging.getLogger('PhaseNet: datalist')
    logger.info("Running to create datalist")
    
    if not os.path.exists(self.datalist_dir):
        os.makedirs(self.datalist_dir)

    if all_in_folder:
        datalist = os.path.join(self.datalist_dir, 'fname.csv')
        make_phasenet_datalist(
            datadir=self.all_mseed,
            datalist_path=datalist,
            json_path=self.json_path,  # Assuming json_path is an attribute
            groupby='{network}.{station}'
        )
    else:
        for sta in os.listdir(self.mseed_storage):
            datadir = os.path.join(self.mseed_storage, sta)
            datalist_path = os.path.join(self.datalist_dir, sta, 'fname.csv')
            make_phasenet_datalist(
                datadir=datadir,
                datalist_path=datalist_path,
                json_path=self.json_path,  # Assuming json_path is an attribute
                groupby='{network}.{station}.{location}'
            )
    
    toc = time.time()
    exetime = timedelta(seconds=toc - tic)
    logger.info(f'Total time of execution: {exetime.total_seconds()} seconds')


def phasenet_from_console(pnet_obj, msg_author):
    """Run PhaseNet from console with specified parameters.
    
    Args:
        pnet_obj: Object with PhaseNet configuration attributes
        msg_author (str): Author identifier for logging
    """
    run = os.path.join(pnet_obj.phasenet_path, "run.py")
    base_command = (
        f"python \"{run}\" --mode=pred "
        f"--model_dir=\"{pnet_obj.model_dir}\" --data_dir=\"{pnet_obj.data_dir}\" "
        f"--data_list=\"{pnet_obj.data_list}\" --output_dir=\"{pnet_obj.output_dir}\" "
        f"--batch_size={pnet_obj.batch_size} --tp_prob={pnet_obj.tp_prob} "
        f"--ts_prob={pnet_obj.ts_prob}"
    )
    
    command = (
        f"{base_command} --input_mseed" if os.name == 'nt' else
        f"{base_command} --resampling={pnet_obj.one_single_sampling_rate} --input_mseed"
    )
    
    if pnet_obj.plot_figure:
        command += " --plot_figure"
    if pnet_obj.save_result:
        command += " --save_result"
    
    printlog("info", msg_author, "Console-> " + command)
    printlog("info", msg_author, "Running...")
    os.system(command)


def get_pickframe(path, csv_filename):
    """Collect picks from CSV files in a directory structure.
    
    Args:
        path (str): Root directory containing station folders with picks
        csv_filename (str): Name of CSV file containing picks
        
    Returns:
        pandas.DataFrame: Combined picks DataFrame
    """
    pick_paths = [
        os.path.join(dp, f) for dp, _, filenames in os.walk(path)
        for f in filenames if f == csv_filename
    ]
    
    dfs = [pd.read_csv(pick_path, index_col=None, header=0) for pick_path in pick_paths]
    return dfs[0] if len(dfs) == 1 else pd.concat(dfs, axis=0, ignore_index=True)


def merge_picks(path, csv_filename, output_path, sort=None):
    """Merge picks from multiple CSV files into one.
    
    Args:
        path (str): Root directory containing station folders with picks
        csv_filename (str): Name of CSV file containing picks
        output_path (str): Path for output merged CSV
        sort (str, optional): Column name to sort by
    """
    pickframe = get_pickframe(path, csv_filename)
    if sort:
        pickframe = pickframe.sort_values(sort)
    
    output_dir = os.path.dirname(output_path)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    pickframe.to_csv(output_path, index=False)


def rm_phasenet_duplicate_picks(path, output_path):
    """Remove duplicate picks from PhaseNet output.
    
    Args:
        path (str): Path to input picks CSV
        output_path (str): Path for output cleaned CSV
        
    Returns:
        pandas.DataFrame: Cleaned picks DataFrame
    """
    def select_pick(df):
        """Select picks based on segment type and probability."""
        logger = logging.getLogger('PhaseNet: Select_picks')
        
        df_overlap = df[df['segment_type'] == 'overlap']
        df_single = df[df['segment_type'] == 'single']
        
        if df_overlap.empty:
            logger.debug('Single selected')
            return df_single
        if df_single.empty:
            logger.debug('Overlap selected')
            return df_overlap
        
        single_prob_avg = np.average(df_single['probability'])
        overlap_prob_avg = np.average(df_overlap['probability'])
        
        choice = df_single if single_prob_avg >= overlap_prob_avg else df_overlap
        logger.debug(
            f"{'Single' if choice is df_single else 'Overlap'} selected: "
            f"Single:{single_prob_avg} vs Overlap:{overlap_prob_avg}"
        )
        return choice

    df = pd.read_csv(path, dtype={"location": str})
    if df.empty:
        return df

    timecols = ['arrival_time', 'creation_time', 'mseed_start_time', 'mseed_end_time']
    df[timecols] = df[timecols].apply(pd.to_datetime)
    
    group_cols = ['mseed_start_time', 'mseed_end_time', 'network', 'station', 'instrument_type']
    mseed_groups = df.groupby(group_cols, as_index=False)
    
    picks_ok = []
    for ind in mseed_groups.indices.keys():
        mseed_df = mseed_groups.get_group(ind)
        mseed_starttime = mseed_df['mseed_start_time'].iloc[0]
        mseed_endtime = mseed_df['mseed_end_time'].iloc[0]
        mseed_sr = mseed_df['sampling_rate'].iloc[0]
        
        overlap_in_sec = 1500 / mseed_sr
        overlap_times = get_chunktimes(mseed_starttime, mseed_endtime, overlap_in_sec, 0)
        
        for ovl_starttime, ovl_endtime in overlap_times:
            ovl_df = mseed_df[
                (mseed_df['arrival_time'] <= ovl_endtime) &
                (mseed_df['arrival_time'] >= ovl_starttime)
            ]
            if not ovl_df.empty:
                picks_ok.append(select_pick(ovl_df))
    
    df_picks_ok = pd.concat(picks_ok, ignore_index=True, sort=False)
    df_picks_ok = df_picks_ok.sort_values('arrival_time').reset_index(drop=True)
    df_picks_ok.to_csv(output_path, index=False)
    return df_picks_ok


def picks2df(picks):
    """Convert list of Pick objects to DataFrame.
    
    Args:
        picks (list): List of Pick objects
        
    Returns:
        pandas.DataFrame: Converted DataFrame
    """
    df_data = []
    for pick in picks:
        comment = json.loads(pick.comments[0].text)
        seismonitor = {
            "pick_id": pick.resource_id.id,
            "arrival_time": pick.time.strftime("%Y-%m-%d %H:%M:%S.%f"),
            "probability": pick.time_errors.uncertainty,
            "phasehint": pick.phase_hint,
            "network": pick.waveform_id.network_code,
            "station": pick.waveform_id.station_code,
            "instrument_type": pick.waveform_id.channel_code[0:2],
            "author": pick.method_id.id,
            "creation_time": pick.creation_info.creation_time.strftime("%Y-%m-%d %H:%M:%S.%f"),
            # "event_start_time": comment.get("event_start_time"),
            # "event_end_time": comment.get("event_end_time"),
            # "mseed_start_time": comment.get("mseed_start_time"),
            # "mseed_end_time": comment.get("mseed_end_time"),
            # "detection_probability": comment.get("detection_probability"),
            # "snr": comment.get("snr"),
            # "sampling_rate": comment.get("sampling_rate"),
            # "sample": comment.get("sample"),
            # "segment": comment.get("segment"),
            # "segment_type": comment.get("segment_type"),
            # "station_lat": comment.get("station_lat"),
            # "station_lon": comment.get("station_lon"),
            # "station_elv": comment.get("station_elv"),
            # "file_name": comment.get("file_name")
        }
        
        for key, value in comment.items():
            seismonitor[key] = value
        
        loc = pick.waveform_id.location_code
        seismonitor["location"] = (
            "" if loc in [None, ""] or (isinstance(loc, (int, float)) and math.isnan(loc))
            else f"{int(loc):02d}" if isinstance(loc, (int, float)) else str(loc)
        )
        seismonitor.update(pick.comments[0])
        df_data.append(seismonitor)
    
    df = pd.DataFrame(df_data)
    df["location"] = df["location"].astype(str)
    date_cols = ["arrival_time", "creation_time", "event_start_time", "event_end_time",
                "mseed_start_time", "mseed_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)
    return df


def get_picks(
    datapicks, datalist, min_p_prob=0.3, min_s_prob=0.3,
    one_single_sampling_rate=-1, mode='df_obj', export=None
):
    """Read PhaseNet picks and convert to desired format.
    
    Args:
        datapicks (str): Path to PhaseNet picks CSV
        datalist (str): Path to datalist CSV
        min_p_prob (float): Minimum P-phase probability
        min_s_prob (float): Minimum S-phase probability
        one_single_sampling_rate (float): Override sampling rate (-1 for auto)
        mode (str): Output mode ('pick_obj' or 'df_obj')
        export (str, optional): Path to export DataFrame as CSV
        
    Returns:
        list or pandas.DataFrame: Picks in specified format
    """
    df = pd.read_csv(datapicks).astype({'itp': 'string', 'tp_prob': 'string',
                                       'its': 'string', 'ts_prob': 'string'})
    df = df.replace("[]", np.nan).dropna(subset=['itp', 'tp_prob', 'its', 'ts_prob'], how="all")
    
    picks = []
    def _get_picks(irow):
        i, row = irow
        wf_name = row["fname"]
        
        p_picks = (
            pick_constructor(datalist, row["itp"].strip('[]').strip().split(),
                           row["tp_prob"].strip('[]').strip().split(), wf_name, 'P',
                           min_p_prob, one_single_sampling_rate)
            if not pd.isna(row["itp"]) else []
        )
        
        s_picks = (
            pick_constructor(datalist, row["its"].strip('[]').strip().split(),
                           row["ts_prob"].strip('[]').strip().split(), wf_name, 'S',
                           min_s_prob, one_single_sampling_rate)
            if not pd.isna(row["its"]) else []
        )
        
        picks.append(p_picks + s_picks)

    with cf.ThreadPoolExecutor() as executor:
        executor.map(_get_picks, df.iterrows())
    
    picks = [x for sublist in picks for x in sublist]
    
    if mode == 'pick_obj':
        return picks
    if mode == 'df_obj' or export:
        if not picks:
            return pd.DataFrame()
        picks = picks2df(picks)
        
        logger = logging.getLogger('PhaseNet: picks2df')
        if not isinstance(picks, pd.DataFrame):
            logger.warning('There are no picks to convert.')
            return pd.DataFrame()
        # print(picks)
        pnet_cols = ["mseed_start_time", "mseed_end_time", "sampling_rate", "sample",
                    "segment", "segment_type"]
        picks = picks[SeisMonitor_columns + pnet_cols].sort_values('arrival_time').reset_index(drop=True)
        logger.info('Conversion pick_sample in pick_time ok. See your results in pick_df.csv.')
        
        if export:
            picks.to_csv(export, index=False, date_format="%Y-%m-%d %H:%M:%S.%f")
    return picks


def pick_constructor(datalist, picks, prob, wf_name, ph_type, min_prob, one_single_sampling_rate=-1):
    """Construct Pick objects from PhaseNet data.
    
    Args:
        datalist (str): Path to datalist CSV
        picks (list): Pick sample points
        prob (list): Pick probabilities
        wf_name (str): Waveform name
        ph_type (str): Phase type ('P' or 'S')
        min_prob (float): Minimum probability threshold
        one_single_sampling_rate (float): Override sampling rate
        
    Returns:
        list: List of Pick objects
        
    Raises:
        Exception: If datalist file is missing
    """
    if not os.path.isfile(datalist):
        logging.getLogger("Datalist").error("Error in phasenet_pick2date()")
        raise Exception("No datalist")

    datalist_df = pd.read_csv(datalist, dtype={"location": str})
    datalist_df[['mseed_start_time', 'mseed_end_time']] = datalist_df[
        ['mseed_start_time', 'mseed_end_time']
    ].apply(pd.to_datetime)
    datalist_df = datalist_df.set_index("fname")

    picks_list = []
    if picks != ['']:
        wfpath, segment = wf_name.split(".mseed_")
        wfpath += ".mseed"
        segment = int(segment)
        
        row = datalist_df.loc[wfpath]
        net, sta, loc, ch = row["network"], row["station"], row["location"], row["instrument_type"]
        
        if isinstance(loc, (int, float)) and math.isnan(loc):
            loc = ""
        sampling_rate = (
            one_single_sampling_rate if one_single_sampling_rate != -1 else row["sampling_rate"]
        )
        
        mseed_starttime, mseed_endtime = row['mseed_start_time'], row['mseed_end_time']
        sta_lat, sta_lon, sta_elv = row["sta_lat"], row["sta_lon"], row["sta_elv"]
        dtt = (mseed_endtime - mseed_starttime).seconds
        
        for pick, p_prob in zip(picks, prob):
            p_prob = float(p_prob)
            if p_prob < min_prob:
                continue
            
            pick_time, creation_time, obs = sample2time(pick, mseed_starttime, sampling_rate, segment, dtt)
            pick_id = id_maker(pick_time, net, sta, loc, ch, ph_type)
            evaluation_mode = 'manual' if p_prob >= 0.95 else 'automatic'
            
            pick_obj = Pick(
                resource_id=ResourceIdentifier(id=pick_id),
                time=pick_time,
                time_errors=QuantityError(uncertainty=p_prob, confidence_level=p_prob * 100),
                waveform_id=WaveformStreamID(
                    network_code=net, station_code=sta, location_code=loc,
                    channel_code=ch, resource_uri=ResourceIdentifier(id=pick_id)
                ),
                phase_hint=ph_type,
                evaluation_mode=evaluation_mode,
                creation_info=CreationInfo(author="PhaseNet_picker", creation_time=UTCDateTime.now()),
                method_id="PhaseNet",
                comments=[Comment(text=json.dumps({
                    "event_start_time": None,
                    "event_end_time": None,
                    "detection_probability": None,
                    "snr": None,
                    "station_lat": sta_lat,
                    "station_lon": sta_lon,
                    "station_elv": sta_elv,
                    "sampling_rate": sampling_rate,
                    "file_name": wf_name,
                    "mseed_start_time": mseed_starttime.strftime("%Y-%m-%d %H:%M:%S.%f"),
                    "mseed_end_time": mseed_endtime.strftime("%Y-%m-%d %H:%M:%S.%f"),
                    "sample": pick,
                    "segment": segment,
                    "segment_type": obs
                }))]
            )
            picks_list.append(pick_obj)
    return picks_list


def id_maker(pick_time, net, station, loc, ch, phasehint):
    """Create SeisComP-style pick PublicID.
    
    Args:
        pick_time (UTCDateTime): Pick time
        net (str): Network code
        station (str): Station code
        loc (str): Location code
        ch (str): Channel code
        phasehint (str): Phase hint ('P' or 'S')
        
    Returns:
        str: SeisComP pick PublicID
    """
    date_id = pick_time.strftime('%Y%m%d.%H%M%S.%f')[:-4]
    return f"{date_id}-{'P' if phasehint == 'P' else 'S'}-{net}.{station}.{loc}.{ch}"


def sample2time(sample, to, df, segment, dtt):
    """Convert PhaseNet sample count to time.
    
    Args:
        sample (str): Sample point of pick
        to (datetime): Waveform start time
        df (float): Sampling rate
        segment (int): Segment number
        dtt (float): Waveform duration in seconds
        
    Returns:
        tuple: (pick_time, creation_time, observation_type)
    """
    init_time = to
    obs = 'single'
    if segment != 0:
        if segment >= dtt * df:
            segment -= dtt * df
            init_time -= datetime.timedelta(seconds=1500 / df)
            obs = 'overlap'
    
    add_seconds = (segment + float(sample)) / df
    pick_time = init_time + dt.timedelta(seconds=add_seconds)
    creation_time = UTCDateTime()
    return pick_time, creation_time, obs