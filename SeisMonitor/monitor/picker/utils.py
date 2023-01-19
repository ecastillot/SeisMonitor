# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:49:23
#  * @modify date 2021-12-22 09:49:23
#  * @desc [description]
#  */

SeisMonitor_columns = ["pick_id","arrival_time","probability","phasehint",
                    "network","station","location","instrument_type","author",
                    "creation_time","event_start_time","event_end_time",
                    "detection_probability","snr","station_lat","station_lon",
                    "station_elv","file_name"]

                    #phasenet other cols
                    # "mseed_start_time","mseed_end_time","sampling_rate","sample",
                    # "segment","segment_type"]

"""
Utils file to use EQTransformer Monitor
"""


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
from SeisMonitor.utils import printlog,isfile
from SeisMonitor.monitor.downloader.utils import get_chunktimes
import concurrent.futures as cf
from obspy.core.event.origin import Pick
from obspy.core.event.base import (QuantityError,
                                WaveformStreamID,
                                CreationInfo,
                                Comment)
from git import Repo
from obspy.core.event import ResourceIdentifier


def git_clone_aipicker(name,repo_dir):
    """
    Params:
    -------
    name: str
        EQTransformer or PhaseNet
    repo_dir: str
        Directory path to place the repository
    """
    if name == "PhaseNet":
        git_url = "https://github.com/ecastillot/PhaseNet.git"
    elif name == "EQTransformer":
        git_url = "https://github.com/ecastillot/EQTransformer.git"
    else:
        return False

    repo_dir = os.path.join(repo_dir,name)
    if os.path.isdir(repo_dir):
        print("There is alaready")
        return True
    else:
        Repo.clone_from(git_url, repo_dir)
        return True

## EQTransformer functions
def eqt_picks_2_seismonitor_fmt(eqt_folder,mseed_folder,out_path):
    #eqt
    date_cols = ["p_arrival_time","s_arrival_time",
                "event_start_time","event_end_time"]

    dfs = []
    for dp, dn, filenames in os.walk(eqt_folder):
        for f in filenames:
            if f == "X_prediction_results.csv" :
                search_path = os.path.join(dp, f)
                df = pd.read_csv(search_path)
                if not df.empty:
                    df[date_cols] = df[date_cols].apply(pd.to_datetime)
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
## PhaseNet functions
def get_filenames(mseed, filter_net=[],
                filter_sta=[],filter_cha=[]):
    """
    Parameters:
    -----------
    mseed: str
        path of the folder where you need to get all mseed filenames.
    filter_net: list
        avoid mseed if its filename contain one network specified in
        filter net
    filter_sta: list
        avoid mseed if its filename contain one station specified in
        filter net
    filter_cha: list
        avoid mseed if its filename contain one channel specified in
        filter net
    Returns:
    --------
    list
    Get all filenames contained in mseed folder.
    """

    if not isinstance(filter_net, list):
        raise Exception("filter_net is a list")
    elif not isinstance(filter_sta, list):
        raise Exception("filter_sta is a list")
    elif not isinstance(filter_cha, list):
        raise Exception("filter_cha is a list")

    filenames = []
    file_list = [ev for ev in os.listdir(mseed) \
                if ev.split('/')[-1].split('.')[-1].lower() == 'mseed']
    for onefile in file_list:
        filename = onefile.split('/')[-1]
        info = filename.split('__')[0]
        net,sta,loc,cha = info.split('.')
        if net in filter_net:
            pass
        if sta in filter_sta:
            pass
        if len(filter_cha) != 0:
            if cha not in filter_cha:
                pass
        else:
            filenames.append(filename)
    return filenames

def mv_mseed2onefolder(mseed_folder,one_folder):
    # print("aca",mseed_folder)
    # print("aca",one_folder)
    for dp, dn, filenames in os.walk(mseed_folder):
        for f in filenames:
            if f.endswith(".mseed") :
                mseed_file = os.path.join(dp, f)
                moved_mseed_file = os.path.join(one_folder, f)
                # os.copy()
                # print(mseed_file, moved_mseed_file)
                printlog("debug","mv_mseed2onefolder",f"{mseed_file} -> {moved_mseed_file}")
                shutil.move(mseed_file, moved_mseed_file)

def mv_mseed2stationfolder(one_folder,mseed_folder):
    for dp, dn, filenames in os.walk(one_folder):
        for f in filenames:
            if f.endswith(".mseed") :
                mseed_file = os.path.join(dp, f)

                src_id = f.split("__")[0]
                net,sta,loc,cha = src_id.split(".")

                moved_mseed_file = os.path.join(mseed_folder, sta,f)
                if not os.path.isdir(os.path.dirname(moved_mseed_file)):
                    os.makedirs(os.path.dirname(moved_mseed_file))
                # os.copy()
                printlog("debug","mv_mseed2stationfolder",f"{mseed_file} -> {moved_mseed_file}")
                # print(mseed_file, moved_mseed_file)
                shutil.move(mseed_file, moved_mseed_file)

def make_dataframe(mseed,json_path, filter_net=[],
                filter_sta=[],filter_cha=[]):
    with open(json_path, "r") as jfile:
        datajson = json.load(jfile)

    
    filenames = get_filenames(mseed, filter_net,
                filter_sta,filter_cha)
    logger = logging.getLogger('PhaseNet: datalist')
    
    mseed_paths = []
    for onefile in filenames:
        df = {}
        path = os.path.join(mseed,onefile)
        st = obspy.read(path)
        stats = st[0].stats

        for key,val in st._get_common_channels_info().items():
            net,sta,loc,instrument_type = key
            channels = list(val["channels"].keys())

        df["fname"] = path

        for channel in channels:
            direction = channel[-1]
            if len(channels) == 1:
                df["E"] = channel
                df["N"] = channel
                df["Z"] = channel
            else:
                df[direction] = channel
        

        instrument = instrument_type[:-1]
        instruments = list(map(lambda x: x[:-1], datajson[sta]['channels']))
        my_samplingrate = datajson[sta]['sampling_rate']
        idxs = [i for i,x in enumerate(instruments) if x == instrument]
        if not idxs:
            continue
        else:
            sampling_rate = my_samplingrate[idxs[0]]
        lat,lon,elv = datajson[sta]['coords']

        df['mseed_start_time'] = stats.starttime.strftime("%Y-%m-%d %H:%M:%S.%f")
        df['mseed_end_time'] = stats.endtime.strftime("%Y-%m-%d %H:%M:%S.%f")
        df["network"] = net
        df["station"] = sta
        df["location"] = loc
        df["instrument_type"] = instrument
        df["sampling_rate"] = sampling_rate
        df["sta_lat"] = lat
        df["sta_lon"] = lon
        df["sta_elv"] = elv

        mseed_paths.append(df)
    df = pd.DataFrame(mseed_paths)
    return df

def make_PhaseNet_datalist(datadir, json_path, datalist_path, channel_list=[], 
                            filter_network=[], filter_station=[],**kwargs ):

    df = make_dataframe(datadir,json_path, filter_network,
                        filter_station,channel_list)
    datalist_dir = os.path.dirname(datalist_path)
    isfile(datalist_path)
    df.to_csv(datalist_path,index=False)

    return datalist_dir

def create_datalist(self,all_in_folder=False):
    tic = time.time()
    logger = logging.getLogger('PhaseNet: datalist')
    logger.info("Running to create datalist")
    
    if not os.path.exists(self.datalist_dir):
        os.makedirs(self.datalist_dir)

    if all_in_folder:
        datalist = os.path.join(self.datalist_dir,'fname.csv')
        datalist = make_PhaseNet_datalist(datadir=self.all_mseed, 
                                datalist_path=datalist, 
                                groupby='{network}.{station}')

    else:

        stations = [ sta for sta in os.listdir(self.mseed_storage)]

        for sta in stations:
            datadir = os.path.join(self.mseed_storage,sta)
            datalist_path = os.path.join(self.datalist_dir,sta,'fname.csv')
            datalist = make_PhaseNet_datalist(datadir=datadir, 
                                    datalist_path=datalist_path, 
                                    groupby='{network}.{station}.{location}')
    
    toc = time.time()
    exetime = timedelta(seconds=toc-tic)
    logger.info(f'Total time of execution: {exetime.total_seconds()} seconds')

def phasenet_from_console(pnet_obj,msg_author):
    phasenet_dir = os.path.dirname(os.path.dirname(pnet_obj.model_dir))
    run = os.path.join(phasenet_dir,"run.py")

    ### nt for windows
    if os.name == 'nt':
        command = f"python \"{run}\" --mode=pred \
            --model_dir=\"{pnet_obj.model_dir}\" --data_dir=\"{pnet_obj.data_dir}\" \
            --data_list=\"{pnet_obj.data_list}\" --output_dir=\"{pnet_obj.output_dir}\"\
            --batch_size={pnet_obj.batch_size} --tp_prob={pnet_obj.tp_prob}\
            --ts_prob={pnet_obj.ts_prob} --input_mseed"

    else:
        command = f"python {run} --mode=pred \
            --model_dir={pnet_obj.model_dir} --data_dir={pnet_obj.data_dir} \
            --data_list={pnet_obj.data_list} --output_dir={pnet_obj.output_dir}\
            --batch_size={pnet_obj.batch_size} --tp_prob={pnet_obj.tp_prob}\
            --ts_prob={pnet_obj.ts_prob} --resampling={pnet_obj.one_single_sampling_rate} --input_mseed"

    if pnet_obj.plot_figure == True:
        command += ' ' + '--plot_figure' +' '+'' 
    if pnet_obj.save_result == True:
        command += ' ' + '--save_result'
    # print(command)
    printlog("info",msg_author,"Console-> "+command)
    printlog("info",msg_author,"Running...")
    os.system(command)

def get_pickframe(path,csv_filename):
    """
    Parameters:
    -----------
    path: str
        Path of the directory that contain the stations folders with
        the picks.
    csv_filename: str
        Name of the file that contains the picks in each station folder
    """
    pick_paths = []
    for dp, dn, filenames in os.walk(path):
        for f in filenames:
            if f == csv_filename:
                pick_path = os.path.join(dp, f)
                pick_paths.append(pick_path)
    li = []
    for pick_path in pick_paths:
        df = pd.read_csv(pick_path, index_col=None, header=0)
        li.append(df)

    if len(li) == 1:
        pickframe = li[0]
    else:
        pickframe = pd.concat(li, axis=0, ignore_index=True)
    return pickframe

def merge_picks(path, csv_filename, output_path, sort=None):
    """
    Parameters:
    -----------
    path: str
        Path of the directory that contain the stations folders with
        the picks.
    csv_filename: str
        Name of the file that contains the picks in each station folder
    """
    pickframe = get_pickframe(path,csv_filename)
    if sort != None:
        pickframe = pickframe.sort_values(sort)
    
    if os.path.isdir(os.path.dirname(output_path)) == False:
        os.makedirs(os.path.dirname(output_path))

    pickframe.to_csv(output_path,index=False)

def rm_phasenet_duplicate_picks(path,output_path):
    
    def select_pick(df):
        logger = logging.getLogger('PhaseNet: Select_picks')

        def average(df):
            probs = df['probability'].to_numpy()
            return np.average(probs)

        df_overlap = df[ df['segment_type'] == 'overlap']
        df_single = df[ df['segment_type'] == 'single']

        if df_overlap.empty == True:
            logger.debug(f'Single selected')
            return df_single
        elif df_single.empty == True:
            logger.debug(f'Overlap selected')
            return df_overlap
        else:
            single_prob_avg = average(df_single)
            overlap_prob_avg = average(df_overlap)

            if single_prob_avg >= overlap_prob_avg:
                logger.debug(f'Single:{single_prob_avg} >= Overlap:{overlap_prob_avg} | Single selected')
                return df_single
            else:
                logger.debug(f'Overlap:{overlap_prob_avg} > Single:{single_prob_avg} | Overlap selected')
                return df_overlap
    
    df = pd.read_csv(path,dtype={"location":str})
    if df.empty == True:
        return df

    # df = df[ df['time'].between(starttime,endtime) ]
    cols = ['mseed_start_time','mseed_end_time','sampling_rate','segment','sample','segment_type']
    strftime = "%Y-%m-%d %H:%M:%S.%f"
    timecols = ['arrival_time','creation_time','mseed_start_time','mseed_end_time']
    df[timecols] = df[timecols].apply(pd.to_datetime)


    mygroup = ['mseed_start_time','mseed_end_time','network',\
                'station','instrument_type']
    # mygroup = ['mseed_start_time','mseed_end_time','network',\
    #             'station','instrument_type','author','snr',
    #             'station_lat','station_lon','station_elv']

    mseed_groups = df.groupby(mygroup,as_index=False)
    indexes = list(mseed_groups.indices.keys())

    picks_ok = []
    for ind in indexes:
        mseed_df = mseed_groups.get_group(ind)
        mseed_starttime = mseed_df['mseed_start_time'].iloc[0] 
        mseed_endtime = mseed_df['mseed_end_time'].iloc[0] 
        mseed_sr = mseed_df['sampling_rate'].iloc[0]

        # 1500 samples (phasenet) / sampling_rate of the mseed file
        overlap_in_sec =  1500 / mseed_sr 

        overlap_times = get_chunktimes(mseed_starttime,mseed_endtime,
                        overlap_in_sec,0)

        for ovl_time in overlap_times:
            ovl_starttime, ovl_endtime = ovl_time
            ovl_df = mseed_df[ (mseed_df['arrival_time'] <= ovl_endtime) & \
                               (mseed_df['arrival_time'] >= ovl_starttime)  ]
            if ovl_df.empty == True:
                pass
            else:
                pick_ok = select_pick(ovl_df)
                picks_ok.append(pick_ok)
                
    df_picks_ok = pd.concat(picks_ok,ignore_index=True,sort=False)
    df_picks_ok = df_picks_ok.sort_values('arrival_time')
    df_picks_ok = df_picks_ok.reset_index(drop=True)
    df_picks_ok.to_csv(output_path,index=False)
    
    return df_picks_ok

def picks2df(picks):
    df = []
    for pick in picks:
        seismonitor = {}
        
        seismonitor["pick_id"] = pick.resource_id.id
        seismonitor["arrival_time"] = pick.time.strftime("%Y-%m-%d %H:%M:%S.%f")
        seismonitor["probability"] = pick.time_errors.uncertainty
        seismonitor["phasehint"] = pick.phase_hint 
        seismonitor["network"] = pick.waveform_id.network_code
        seismonitor["station"] = pick.waveform_id.station_code
        seismonitor["instrument_type"] = pick.waveform_id.channel_code[0:2]
        seismonitor["author"] = pick.method_id.id
        seismonitor["creation_time"] = pick.creation_info.creation_time.strftime("%Y-%m-%d %H:%M:%S.%f")
        
        if isinstance(pick.waveform_id.location_code,int):
            if math.isnan(pick.waveform_id.location_code):
                seismonitor["location"] = ""
        elif pick.waveform_id.location_code == None:
            seismonitor["location"] = ""
        else:
            seismonitor["location"] = "{:02d}".format(int(pick.waveform_id.location_code))

        for key,val in pick.comments[0].items():
            seismonitor[key] = val

        df.append(seismonitor)
    df = pd.DataFrame(df)
    df["location"] = df["location"].astype(str)
    # print(df["location"])
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time",
                "mseed_start_time","mseed_end_time"]
    df[date_cols] = df[date_cols].apply(pd.to_datetime)
    return df

def get_picks(datapicks, datalist, 
             min_p_prob=0.3,
             min_s_prob=0.3,
             one_single_sampling_rate = -1,
            mode='df_obj',
            export=None):   #mode: pick_obj ; df_obj
    """Read phaseNet picks and returns list of Pick objects
    Parameters
    ----------
    datapicks : str
        Path to generated PhaseNet picks csv
    min_prob : float
        Minimum probability to consider a pick.
    
    Returns
    -------
    list
        List of Pick objects
    """

    df = pd.read_csv(datapicks)
    df = df.astype({'itp':'string','tp_prob':'string',
                    'its':'string','ts_prob':'string'})
    df = df.replace("[]",np.nan)
    df = df.dropna(subset=['itp','tp_prob','its','ts_prob'],how="all")

    picks = []
    def _get_picks(irow):
        i,row = irow

        wf_name = row["fname"]
        if not pd.isna(row["itp"]):
            picks_p = row["itp"].strip('[]').strip().split() 
            prob_p = row["tp_prob"].strip('[]').strip().split() 
            P_picks = pick_constructor(datalist,picks_p, 
                                prob_p, wf_name, 'P', min_p_prob,
                                one_single_sampling_rate )
        else:
            P_picks = []
        
        if not pd.isna(row["its"]):
            picks_s = row["its"].strip('[]').strip().split() 
            prob_s = row["ts_prob"].strip('[]').strip().split() 
            S_picks = pick_constructor(datalist,picks_s, 
                            prob_s, wf_name, 'S', min_s_prob,
                            one_single_sampling_rate )
        else:
            S_picks = []

        picks.append(P_picks+S_picks)

    # for irow in df.iterrows():
    #     _get_picks(irow)
    with cf.ThreadPoolExecutor() as executor:
        executor.map(_get_picks,df.iterrows())
    
    picks = [x for n in picks for x in n]

    if mode == 'pick_obj':
        pass
    if (mode ==  'df_obj') or (export != None):
        if not picks:
            return pd.DataFrame()
        picks= picks2df(picks)

        logger = logging.getLogger(f'PhaseNet: picks2df')
        if isinstance(picks, pd.DataFrame) == False:
            logger.warning('There are no picks to convert.')
            return pd.DataFrame()
        else:
            pnet_cols = ["mseed_start_time","mseed_end_time","sampling_rate","sample",
                    "segment","segment_type"]
            picks = picks[SeisMonitor_columns+pnet_cols]
            picks = picks.sort_values('arrival_time')
            picks = picks.reset_index(drop=True)
            logger.info('Conversion pick_sample in pick_time ok.'+
                        ' See your results in pick_df.csv.')
            if export != None:
                picks.to_csv(export,index=False,
                            date_format="%Y-%m-%d %H:%M:%S.%f")
    return picks

def pick_constructor(datalist,picks, prob, wf_name, ph_type, min_prob,
                    one_single_sampling_rate  = -1):
    """Construct Pick objects
    Parameters
    ----------
    picks : list
        List of phase pick sample point in a specific waveform
    prob : list
        List of probabilities associated with the picks
    wf_name : str
        Waveform name
    ph_type : str
        Pick phase type. Can be P or S
    min_prob : str
        Minimum probability to consider a pick
    
    Returns
    -------
    list
        List of picks objects
    """
    ## Add json parameter because need sampling rate of the stations.
    if os.path.isfile(datalist) == False:
        logger = logging.getLogger("Datalist")
        logger.error("Error in phasenet_pick2date()")
        raise Exception("No datalist ")

    else:
        datalist = pd.read_csv(datalist,dtype={"location":str})
        datalist[['mseed_start_time','mseed_end_time']] = datalist[['mseed_start_time','mseed_end_time']].apply(pd.to_datetime)
        datalist = datalist.set_index("fname")


    segment = 0.0
    picks_list = []
    if picks != ['']:
        for pick, prob in zip(picks, prob):
            prob = float(prob)
            if prob >= min_prob:
                

                #----------
                # se obtienen los parámetros para la creación del objeto pick
                #----------
                # algunos datos se obtienen del nombre de la forma de onda
                # net, station, loc, ch, df, *to_segment = wf_name.split('_')
                # datajson[sta]
                wfpath,segment = wf_name.split(".mseed_")
                wfpath = wfpath +".mseed"
                segment = int(segment)

                row = datalist.loc[wfpath]
                net = row["network"]
                sta = row["station"]
                loc = row["location"]
                # print(loc,type(loc))
                ch = row["instrument_type"]
                
                if one_single_sampling_rate == -1:
                    sampling_rate = row["sampling_rate"]
                else:
                    sampling_rate = one_single_sampling_rate

                mseed_starttime = row['mseed_start_time']
                mseed_endtime = row['mseed_end_time']
                sta_lat = row["sta_lat"]
                sta_lon = row["sta_lon"]
                sta_elv = row["sta_elv"]

                dtt = (mseed_endtime - mseed_starttime).seconds

                # se transforma las cuentas asociadas al pick en tiempo
                pick_time, creation_time,obs = sample2time(pick, mseed_starttime, sampling_rate, segment, dtt)
                # se crea el Id usando el tiempo del pick
                ID = id_maker(pick_time, net, sta, loc, ch, ph_type)    

                # Se evalua si la probabilidad es lo suficientemente buena 
                # como para considerarlo manual
                evaluation_mode = 'automatic'
                if prob >= 0.95:
                    evaluation_mode = 'manual'

                # Se crea el objeto Pick
                pick_obj = Pick(resource_id=ResourceIdentifier( id= ID),
                            time=pick_time,
                            time_errors=QuantityError(uncertainty=prob,
                                                    confidence_level=prob*100),
                            waveform_id=WaveformStreamID(network_code=net,
                                                            station_code=sta,
                                                            location_code=loc,
                                                            channel_code=ch,
                                                            resource_uri= ResourceIdentifier( id= ID)),
                            phase_hint = ph_type,
                            evaluation_mode = evaluation_mode,
                            creation_info= CreationInfo(author="PhaseNet_picker",
                                                        creation_time=UTCDateTime.now()),
                            method_id="PhaseNet",
                            comments= [{"event_start_time":None,
                                        "event_end_time":None,
                                        "detection_probability":None,
                                        "snr":None,
                                        "station_lat":sta_lat,
                                        "station_lon":sta_lon,
                                        "station_elv":sta_elv,
                                        "sampling_rate":sampling_rate,
                                        "file_name":wf_name,
                                        "mseed_start_time":mseed_starttime.strftime("%Y-%m-%d %H:%M:%S.%f"),
                                        "mseed_end_time":mseed_endtime.strftime("%Y-%m-%d %H:%M:%S.%f"),
                                        "sample":pick,
                                        "segment":segment,
                                        "segment_type":obs,
                                        }]
                            # comments= [f'{to}to{endtime} sr{df} segment{segment}_sample{pick} {obs}',coords]
                            )
                # Se agrega cada pick a la lista de picks
                picks_list.append(pick_obj)
    return picks_list

def id_maker(pick_time, net, station, loc, ch, phaseHint):
    """Creates the seiscomp Pick PublicID
    
    Parameters
    ----------
    pick_time : Obspy UTCDateTime object
        Time for phase pick
    
    Returns
    ------
    str
       Seiscomp pick PublicID 
    """
    dateID = pick_time.strftime('%Y%m%d.%H%M%S.%f')[:-4]
    if phaseHint == 'P':
        # publicID = dateID+f'-AIC-{net}.{station}.{loc}.{ch}'
        publicID = dateID+f'-P-{net}.{station}.{loc}.{ch}'
    elif phaseHint == 'S':
        # publicID = dateID+f'-S-L2-{net}.{station}.{loc}.{ch}'
        publicID = dateID+f'-S-{net}.{station}.{loc}.{ch}'
    return publicID

def sample2time(sample, to, df, segment, dtt):
    """Transforma las cuentas de un pick de PhaseNet en fecha
    
    Parameters
    ----------
    sample : str
        Sample point of PhaseNet pick
    to : str
        Initial time of the waveform that contains the pick in format
        YYYYmmddHHmmssff. Whit ff as deciseconds
    df : str
        Sampling rate
    
    Returns
    -------
    pick_time : Obspy UTCDateTime object
        Time for phase pick
    creation_time : UTCDateTime
        Time for creation time
    """

    init_time = to
    # if segment is different to 0 which implies that we are using
    # pred_mseed mode
    obs = 'single'
    if segment is not 0:
        
        # if segment is bigger than the lenght of the waveform; then,
        # the segment is an overlapping one, and then we need to
        # include the 1500 samples (15 s) of shiftfing 
        if segment >= dtt*df:
            segment = segment - dtt*df
            init_time -= datetime.timedelta(seconds=1500/df)
            obs = 'overlap'

    add_seconds = (segment + float(sample) )/df
    pick_time = init_time + dt.timedelta( seconds=add_seconds )
    creation_time = UTCDateTime()
    return pick_time, creation_time, obs
