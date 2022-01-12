# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:49:23
#  * @modify date 2021-12-22 09:49:23
#  * @desc [description]
#  */


"""
Utils file to use EQTransformer Monitor
"""


import sys
import json
import obspy
import time
import os
import shutil
import logging
import numpy as np
import pandas as pd
import datetime as dt
from SeisMonitor.utils import printlog
from SeisMonitor.downloader.utils import get_chunktimes
from datetime import timedelta
from obspy.clients.fdsn import Client as FDSN_Client
from obspy import UTCDateTime
from obspy.core.inventory.inventory import (Inventory,read_inventory)

import csv
import datetime
from obspy.core.event.origin import Pick
from obspy.core.event.base import (QuantityError,
                                WaveformStreamID,
                                CreationInfo,
                                Comment)
from obspy.core.event import ResourceIdentifier

class EQTobj(object):
    def __init__(self,model_path,chunk_size=3600,
                n_processor=2,overlap=0.3,
                detection_threshold=0.1, P_threshold=0.1,
                S_threshold=0.1,number_of_plots=1,
                batch_size=1,
                plot_mode=1,
                overwrite=False):

        """
        EQTransformer parameters
        """

        self.model_path = model_path
        self.chunk_size = chunk_size
        self.n_processor = n_processor
        self.overlap = overlap
        self.detection_threshold = detection_threshold
        self.P_threshold = P_threshold
        self.S_threshold = S_threshold
        self.number_of_plots = number_of_plots
        self.batch_size = batch_size
        self.plot_mode = plot_mode
        self.overwrite = overwrite
        self.name = "EQTransformer"

class PhaseNetobj(object):
    def __init__(self, model_path,
                 mode='pred', P_threshold=0.3, S_threshold=0.3,
                batch_size=2, plot = False, save_result=False,
                epochs = 100,learning_rate= 0.01,decay_step = -1,
                decay_rate = 0.9,momentum = 0.9,filters_root = 8,
                depth = 5, kernel_size = [7, 1], pool_size = [4, 1],
                drop_rate = 0, dilation_rate = [1, 1],loss_type = "cross_entropy",
                weight_decay = 0, optimizer = "adam",   summary = True,
                class_weights = [1, 1, 1],log_dir = "log",  num_plots = 10,
                input_length = None,input_mseed = True, filename_picks = "picks",
                data_dir = "./dataset/waveform_pred/",
                data_list = "./dataset/waveform.csv",
                train_dir = "./dataset/waveform_train/",
                train_list ="./dataset/waveform.csv",
                valid_dir = None, valid_list = None,
                output_dir = None ):
        """
        PhaseNet parameters
        """

        self.model_dir = model_path
        self.mode = mode
        self.tp_prob = P_threshold
        self.ts_prob = S_threshold
        self.batch_size = batch_size
        self.plot_figure = plot
        self.save_result = save_result
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.decay_step = decay_step
        self.decay_rate = decay_rate
        self.momentum = momentum
        self.filters_root = filters_root
        self.depth = depth
        self.kernel_size = kernel_size
        self.pool_size = pool_size
        self.drop_rate = drop_rate
        self.dilation_rate = dilation_rate
        self.loss_type = loss_type
        self.weight_decay = weight_decay
        self.optimizer = optimizer
        self.summary = summary
        self.class_weights = class_weights
        self.log_dir = log_dir
        self.num_plots = num_plots
        self.input_length = input_length
        self.input_mseed = input_mseed
        self.fpred = filename_picks
        self.data_dir = data_dir
        self.data_list = data_list
        self.train_dir = train_dir
        self.train_list = train_list
        self.valid_dir = valid_dir
        self.valid_list = valid_list
        self.output_dir = output_dir
        self.name = "PhaseNet"

def filter_inventory(inv,network,station,location,
                    channel,starttime,endtime):
    """
    Filter inventory according to parameters

    Parameters
    ----------
    inv:Inventory
        Inventory that will be filtered
    network: str
        Select one or more network codes. 
        Multiple codes are comma-separated (e.g. "IU,TA"). 
        Wildcards are allowed.
    station: str
        Select one or more SEED station codes. 
        Multiple codes are comma-separated (e.g. "ANMO,PFO"). 
        Wildcards are allowed.
    location: str
        Select one or more SEED location identifiers. 
        Multiple identifiers are comma-separated (e.g. "00,01"). 
        Wildcards are allowed.
    channel: str
        Select one or more SEED channel codes. 
        Multiple codes are comma-separated (e.g. "BHZ,HHZ").
    starttime: obspy.UTCDateTime
        Limit results to time series samples on or 
        after the specified start time.
    endtime: obspy.UTCDateTime
        Limit results to time series samples on or 
        before the specified end time.

    Returns:
    --------
    inventory: Inventory
        Filtered Inventory
    """
    networks = network.split(',')
    stations = station.split(',')
    locations = location.split(',')
    channels = channel.split(',')

    cha0 = channels[0]
    inventory = inv.select(network=networks[0],
                            station=stations[0],
                            location=locations[0],
                            channel=cha0,
                            starttime=starttime,
                            endtime=endtime)
    for net in networks:
        for sta in stations:
            for loc in locations:
                for cha in channels:
                    one_inv = inv.select(network=net,station=sta,
                                            location=loc,channel=cha,
                                            starttime=starttime,
                                            endtime=endtime)
                    inventory = inventory.__add__(one_inv) 
    return inventory

def makeJSON(json_path,providers, restrictions, from_xml=None,
                    channel_list=[], filter_network=[], filter_station=[],**kwargs ):


    """
    
    Uses fdsn to find available stations in a specific geographical location and time period.  
    Parameters
    ----------
    json_path: str
        Path of the json file that will be returned
    providers:list
        list of Client object
    restrictions: DownloadRestrictions or MassDownloaderRestrictions
        Restrictions to download mseed
    from_xml : str
        Path of xml file
    **kwargs: str
        get_stations kwargs


    Returns
    ----------
    stations_list.json: A dictionary containing information for the available stations.      
        
    """  
    station_list = {}

    if from_xml != None:
        inv = read_inventory(from_xml)
        # inv = read_inventory(from_xml, format="STATIONXML")
        inventory = filter_inventory(inv=inv ,network=restrictions.network,
                                    station=restrictions.station,
                                    location=restrictions.location,
                                    channel=restrictions.channel,
                                    starttime=restrictions.starttime, 
                                    endtime=restrictions.endtime)

    else:
        # inventory = Inventory()
        inventory = Inventory()
        for provider in providers:
            one_inv = provider.get_stations(network=restrictions.network,
                                                    station=restrictions.station,
                                                    location=restrictions.location,
                                                    channel=restrictions.channel,
                                                    starttime=restrictions.starttime, 
                                                    endtime=restrictions.endtime, 
                                                    level='channel',**kwargs)
            inventory = inventory.__add__(one_inv) 
    
    logger = logging.getLogger('json') 
    toprint = []
    for ev in inventory:
        net = ev.code
        if net not in filter_network:
            for st in ev:
                station = st.code
                msg = f'{net}-{station}'
                logger.debug(msg)

                if station not in filter_station:

                    elv = st.elevation
                    lat = st.latitude
                    lon = st.longitude

                    new_chan = [ch.code for ch in st.channels]
                    if len(channel_list) > 0:
                        chan_priority=[ch[:2] for ch in channel_list]

                        for chnn in chan_priority:
                            if chnn in [ch[:2] for ch in new_chan]:
                                new_chan = [ch for ch in new_chan if ch[:2] == chnn]     
                

                    if len(new_chan) > 0 and (station not in station_list):
                        channels = list(set(new_chan))
                        sample_rates_gen = lambda x: st.select(channel=x).channels[0].sample_rate
                        sample_rates = list(map(sample_rates_gen,channels))
                        station_list[str(station)] ={"network": net,
                                                "channels": list(set(new_chan)),
                                                "coords": [lat, lon, elv],
                                                "sampling_rate": sample_rates
                                                }
                    toprint.append(msg)
    logger.info(str(toprint) + ' ok')

    json_dir = os.path.dirname(json_path)
    if not os.path.exists(json_dir):
        os.makedirs(json_dir)
    with open(json_path, 'w') as fp:
        json.dump(station_list, fp)

    return json_dir  





## PhaseNet implementation
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

def get_one_stream(input_path,output_path,**kwargs):
    """
    Parameters:
    -----------
    input_path: str
        Where is stored the waveform files
    output_path: str
        Path to move the waveform files in only
        into a single 3D waveform file.

    kwargs: corresponding to get_filenames() parameters

    Returns:
    --------
        Save one stream with all available channels
        in the mseed folder.
    """
    strftime = "%Y%m%dT%H%M%SZ"

    if os.path.isdir(output_path) == False:
        os.makedirs(output_path)

    for station_path in os.listdir(input_path):
        if station_path != 'all_mseed':
            path = os.path.join(input_path,station_path)
            filenames = get_filenames(path,**kwargs)

            st = obspy.core.Stream()
            for onefile in filenames:
                onest = obspy.read(os.path.join(path,onefile))
                st += onest
            if len(st) != 0:
                tr = st[0]
                stats = tr.stats
                network,station,\
                    channel,location,\
                        starttime, endtime = (stats.network,stats.station,\
                                            stats.channel,stats.location,\
                                                stats.starttime.strftime(strftime),\
                                                stats.endtime.strftime(strftime))
                fmt = f"{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
                printlog("info","mv_onefolder: ",f"Done: {network}-{station}")
                streampath = os.path.join(output_path,fmt)
                st.write(streampath)

def get_stname(str_id,network=None,station=None,
                location=None,channel=None):
    if ("{network}" in str_id) and ("{station}" in str_id) and \
            ("{location}" in str_id) and ("{channel}" in str_id):
        
        str_id = str_id.format(network=network,station=station,
                            location=location,channel=channel)

    elif ("{network}" in str_id) and ("{station}" in str_id) and \
            ("{location}" in str_id):
        str_id = str_id.format(network=network,station=station,
                            location=location)
    elif ("{network}" in str_id) and ("{station}" in str_id):
        str_id = str_id.format(network=network,station=station)

    return str_id

def make_dataframe(mseed,groupby, filter_net=[],
                filter_sta=[],filter_cha=[]):
    filename = []
    E = []
    N = []
    Z = []
    
    filenames = get_filenames(mseed, filter_net,
                filter_sta,filter_cha)
    logger = logging.getLogger('PhaseNet: datalist')
    for onefile in filenames:
        C = []
        net,sta,loc,cha = onefile.split('__')[0].split('.')
        path = os.path.join(mseed,onefile)
        st = obspy.read(path)
        st = st._groupby(groupby)
        st_name = get_stname(groupby,net,sta,loc,cha)
        st = st[st_name] #st_by_groupby

        st = st.merge()
        for tr in st:
            channel = tr.stats.channel
            C.append(channel)
        existence = list(map(lambda x: x[-1] \
                            if x[-1] in ('E','N','Z') else '',C))
        if 'E' in existence:    
            e = existence.index('E')
            if 'Z' not in existence:
                Z.append("None")
            if 'N' not in existence:
                N.append("None")
            E.append(C[e])
        if 'N' in existence:    
            n = existence.index('N')
            if 'E' not in existence:
                E.append("None")
            if 'Z' not in existence:
                Z.append("None")
            N.append(C[n])
        if 'Z' in existence:    
            z = existence.index('Z')
            if 'E' not in existence:
                E.append("None")
            if 'N' not in existence:
                N.append("None")
            Z.append(C[z])
    
        logger.info(f'{onefile} ok')
        filename.append(onefile)
    data = {'fname':filename,'E':E,'N':N,'Z':Z}
    df = pd.DataFrame.from_dict(data)
    return df

def make_PhaseNet_datalist(datadir, datalist_path, groupby, channel_list=[], 
                            filter_network=[], filter_station=[],**kwargs ):

    df = make_dataframe(datadir,groupby, filter_network,
                        filter_station,channel_list)
    datalist_dir = os.path.dirname(datalist_path)
    if not os.path.exists(datalist_dir):
        os.makedirs(datalist_dir)
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

def phasenet_from_console(pnet_obj):
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
            --ts_prob={pnet_obj.ts_prob} --input_mseed"

    if pnet_obj.plot_figure == True:
        command += ' ' + '--plot_figure' +' '+'' 
    if pnet_obj.save_result == True:
        command += ' ' + '--save_result'
    # print(command)
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

def rm_phasenet_duplicate_picks(path,output_path,group_by_loc=True):
    def text2col(comments_text):
        text = comments_text.strip("['").strip("']")
        dates,sr,pick,obs = text.split(' ')
        segment,sample = pick.split('_')

        starttime,endtime = dates.split('to')
        starttime = starttime
        endtime = endtime
        try:
            sr = int(float(sr.strip('sr')))
        except:
            sr = None
        segment = int(segment.strip('segment'))
        sample = int(sample.strip('sample'))
        
        return pd.Series({'mseed_starttime':starttime,
                    'mseed_endtime':endtime,
                    'sampling_rate':sr,
                    'segment':segment,
                    'sample':sample,
                    'obs':obs})

    def select_pick(df):
        logger = logging.getLogger('PhaseNet: Select_picks')

        def average(df):
            probs = df['uncertainty'].to_numpy()
            return np.average(probs)

        df_overlap = df[ df['obs'] == 'overlap']
        df_single = df[ df['obs'] == 'single']

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
    
    df = pd.read_csv(path)
    if df.empty == True:
        return df

    # df = df[ df['time'].between(starttime,endtime) ]
    cols = ['mseed_starttime','mseed_endtime','sampling_rate','segment','sample','obs']
    df[cols] = df['comments_text'].apply(lambda x: text2col(x))
    df['mseed_starttime'] = pd.to_datetime(df['mseed_starttime'])
    df['mseed_endtime'] = pd.to_datetime(df['mseed_endtime'])
    df['time'] = pd.to_datetime(df['time'])

    if group_by_loc:
        mygroup = ['mseed_starttime','mseed_endtime','network_code',\
                    'station_code','channel_code','location_code']
    else:
        mygroup = ['mseed_starttime','mseed_endtime','network_code',\
                    'station_code','channel_code']
    mseed_groups = df.groupby(mygroup,as_index=False)
    indexes = list(mseed_groups.indices.keys())
    picks_ok = []
    for ind in indexes:
        mseed_df = mseed_groups.get_group(ind)
        # + 999.9 miliseconds because the mseed filename doesn't have it.
        mseed_starttime = mseed_df['mseed_starttime'].iloc[0] +\
                            dt.timedelta(milliseconds=999.9)
        mseed_endtime = mseed_df['mseed_endtime'].iloc[0] +\
                            dt.timedelta(milliseconds=999.9) 
        mseed_sr = mseed_df['sampling_rate'].iloc[0]

        # 15000 samples (phasenet) / sampling_rate of the mseed file
        overlap_in_sec =  1500 / mseed_sr 

        overlap_times = get_chunktimes(mseed_starttime,mseed_endtime,
                        overlap_in_sec,0)

        for ovl_time in overlap_times:
            ovl_starttime, ovl_endtime = ovl_time
            ovl_df = mseed_df[ (mseed_df['time'] <= ovl_endtime) & \
                               (mseed_df['time'] >= ovl_starttime)  ]
            if ovl_df.empty == True:
                pass
            else:
                pick_ok = select_pick(ovl_df)
                picks_ok.append(pick_ok)

    df_picks_ok = pd.concat(picks_ok,ignore_index=True,sort=False)
    df_picks_ok = df_picks_ok.sort_values('time')
    df_picks_ok = df_picks_ok.drop(['Unnamed: 0'],axis=1)
    df_picks_ok = df.reset_index(drop=True)
    df_picks_ok.to_csv(output_path,index=False)
    
    return df_picks_ok

## pick2date
def picks2df(picks):
    appended_picks = []
    for pick in picks:
        appended_items = []
        for pick_key,pick_value in pick.items():
            if pick_key in ('time_errors','waveform_id',
                            'horizontal_slowness_errors',
                            'backazimuth_errors',
                            'creation_info'):

                if pick_value != None:
                    for new_key, new_item in pick_value.items():
                        appended_items.append((new_key, new_item))

            elif pick_key in ('comments'):
                if not pick_value:
                    appended_items.append((pick_key,None))
                else:
                    appended_items.append((pick_key+'_text',pick_value))
            else:
                appended_items.append((pick_key,pick_value))
            
        pick_dict = dict(appended_items)
        pick_df = pd.DataFrame( pick_dict.values(), 
                                index= pick_dict.keys() ).T
        appended_picks.append(pick_df)

    if not appended_picks:
        return None
    df_picks = pd.concat(appended_picks,ignore_index=True,sort=False)

    return df_picks

def get_picks(phaseNet_picks, jsonfile, dt,
             min_prob=0.3, mode='df_obj', export='csv'):   #mode: pick_obj ; df_obj
    '''Read phaseNet picks and returns list of Pick objects
    Parameters
    ----------
    phaseNet_picks : str
        Path to generated PhaseNet picks csv
    min_prob : float
        Minimum probability to consider a pick.
    
    Returns
    -------
    list
        List of Pick objects
    '''
    """
    # script directory for phaseNet.inp searching
    main_dir = os.path.dirname(os.path.abspath(__file__))
    par_fn = 'phaseNet.inp'
    rel_par_path = os.path.join('../', par_fn)
    main_par_path = os.path.join(main_dir, par_fn)
    # verifying if phaseNet.inp exist in any of the following 3 paths. 
    check_inp_dirs = [par_fn, rel_par_path, main_par_path]
    for path in check_inp_dirs:
        if os.path.isfile(path):
            print(f'Reading params for: {path} \n')
            params = read_params(path)
            break
    
    mode = params['mode']
    """
    export_dirname = os.path.dirname(phaseNet_picks)

    picks = []
    with open(phaseNet_picks, newline='') as csvfile: 
        reader = csv.reader(csvfile, delimiter=',') 
        for i, row in enumerate(reader): 
            if i != 0: 
                wf_name = row[0]
                picks_p = row[1].strip('[]').strip().split() 
                prob_p = row[2].strip('[]').strip().split() 
                picks_s = row[3].strip('[]').strip().split() 
                prob_s = row[4].strip('[]').strip().split() 
                P_picks = pick_constructor(jsonfile,picks_p, prob_p, wf_name, 'P', min_prob, dt)
                S_picks = pick_constructor(jsonfile,picks_s, prob_s, wf_name, 'S', min_prob, dt)
                picks += P_picks + S_picks


    if mode == 'pick_obj':
        pass
    if (mode ==  'df_obj') or (export != None):
        picks= picks2df(picks)

        logger = logging.getLogger(f'PhaseNet: picks2df')
        if isinstance(picks, pd.DataFrame) == False:
            logger.warning('There are no picks to convert.')
            return None
        else:
            picks = picks.sort_values('time')
            picks = picks.reset_index(drop=True)
            logger.info('Conversion pick_sample in pick_time ok.'+
                        ' See your results in pick_df.csv.')
            if export != None:
                export_file = os.path.join(export_dirname,'picks_df')
                if export == 'csv':
                    picks.to_csv(f'{export_file}.csv')
                elif export == 'excel':
                    picks.to_excel(f'{export_file}.xlsx')
                else:
                    logger.error('Format not supported!; only csv or excel')
    return picks

def pick_constructor(jsonfile,picks, prob, wf_name, ph_type, min_prob, dt):
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
    if os.path.isfile(jsonfile) == False:
        logger = logging.getLogger("JSON")
        logger.error("Error in phasenet_pick2date()")
        logger.error("You need to create json_file to get the sampling rate of the stations")
        raise Exception("Check makeJSON in AIpicker.picker.utils2pick ")

    with open(jsonfile, "r") as jfile:
        datajson = json.load(jfile)

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
                # print(wf_name)
                str_id,to,end_name = wf_name.split('__')
                endtime,segment = end_name.split('_')
                endtime = endtime.split('.mseed')[0]
                segment = int(segment)
                net,station,loc,ch = str_id.split('.')
                # print(str_id)
                my_cha = datajson[station]['channels']
                my_samplingrate = datajson[station]['sampling_rate']
                # print(ch)
                idxs = [i for i,x in enumerate(my_cha) if x == ch]
                if not idxs:
                    continue
                else:
                    df = my_samplingrate[idxs[0]]


    #             to = to_segment[0].split('.')[0]
                # if len(to_segment)==2:
                #     segment = int(to_segment[1])

                # se transforma las cuentas asociadas al pick en tiempo
                pick_time, creation_time,obs = sample2time(pick, to, df, segment, dt)
                # se crea el Id usando el tiempo del pick
                ID = id_maker(pick_time, net, station, loc, ch, ph_type)    

                # Se evalua si la probabilidad es lo suficientemente buena 
                # como para considerarlo manual
                evaluation_mode = 'automatic'
                if prob >= 0.95:
                    evaluation_mode = 'manual'

                # print(to)
                # Se crea el objeto Pick
                pick_obj = Pick(resource_id=ResourceIdentifier( id= ID),
                                time=pick_time,
                                time_errors=QuantityError(uncertainty=prob,
                                                        confidence_level=prob*100),
                                waveform_id=WaveformStreamID(network_code=net,
                                                             station_code=station,
                                                             location_code=loc,
                                                             channel_code=ch,
                                                             resource_uri= ResourceIdentifier( id= ID)),
                                phase_hint = ph_type,
                                evaluation_mode = evaluation_mode,
                                creation_info= CreationInfo(author="PhaseNet_picker",
                                                            creation_time=UTCDateTime.now()),
                                method_id="PhaseNet",
                                comments= [f'{to}to{endtime} sr{df} segment{segment}_sample{pick} {obs}']
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
        publicID = dateID+f'-AIC-{net}.{station}.{loc}.{ch}'
    elif phaseHint == 'S':
        publicID = dateID+f'-S-L2-{net}.{station}.{loc}.{ch}'
    return publicID

def sample2time(sample, to, df, segment, dt):
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
    df = float(df)

    # init_time = UTCDateTime(to[:-2]+'.'+to[-2:])
    init_time = UTCDateTime(to)
    # if segment is different to 0 which implies that we are using
    # pred_mseed mode
    obs = 'single'
    if segment is not 0:
        
        # if segment is bigger than the lenght of the waveform; then,
        # the segment is an overlapping one, and then we need to
        # include the 1500 samples (15 s) of shiftfing 
        if segment >= dt*df:
            if df == 0:
                df = 100
            segment = segment - dt*df
            init_time -= datetime.timedelta(seconds=1500/df)
            obs = 'overlap'

    # +100 because in the mseed name, the start time is -1 second
    pick_time = init_time + (segment + float(sample) )/df 
    creation_time = UTCDateTime()
    return pick_time, creation_time, obs
