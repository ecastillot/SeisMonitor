"""
 * @author Emmanuel David Castillo Taborda
 * @email ecastillot@unal.edu.co
 * @create date 2021-03-03 23:51:21
 * @modify date 2021-06-02 04:31:09
 * @desc [description]
"""

## not available yet

import datetime as dt
import time
import os
import json
import pandas as pd
from .utils import Pick,prepare_eqt,picks2xml
import concurrent.futures
from obspy import UTCDateTime
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.catalog import read_events
from obspy.core.event import Catalog

def merge_csv(storage,ai_picker,csv_type="event",sort=None,export=None):
    """
    Parameters:
    -----------
    storage: str
        Path of the storage folder
    ai_picker: str
        "phasenet" or "eqt"
    csv_type: str
        "event" or "pick"
    sort: str
        sort by a key word.
    export: str
        path of the csv file

    Returns:
    --------
        Dataframe with all information merged.
    """
    if ai_picker in ("eqt","EQT","EQTransformer"):
        ai_picker = "eqt"
        events_csv_name = "eqt_events.csv"
        # picks_csv_name = "eqt_picks.csv"
        picks_csv_name = "X_prediction_results.csv"

    elif ai_picker in ("pnet","phasenet","PhaseNet"):
        ai_picker = "phasenet"
        events_csv_name = "phasenet_events.csv"
        picks_csv_name = "phasenet_picks.csv"
    
    if csv_type == "event":
        search = events_csv_name 
    elif csv_type == "pick":
        search = picks_csv_name

    events = []
    for dp, dn, filenames in os.walk(storage):
        for f in filenames:
            if f == search:
                search_path = os.path.join(dp, f)
                print(search_path)
                df = pd.read_csv(search_path, index_col=0)
                events.append(df)

    events_df = pd.concat(events,ignore_index=True)
    if sort != None:   
        events_df = events_df.sort_values(by=sort,ascending=True,ignore_index=True)
    
    events_df.index+=1
    events_df.index.name= 'No'

    if export != None:
        if os.path.isdir(os.path.dirname(export)) == False:
            os.makedirs(os.path.dirname(export))
        events_df.to_csv(export)
        print(f"Events_csv_file: {export}")
    return events_df

def filter_phasenet_probs(df,phase_hint,json_path):
    """
    Parameters:
    -----------
    df : DataFrame
        PhaseNet dataframe
    phase_hint: str
        "P" or "S"
    json_path: str
        path where is located the filter probabilities.

    Returns:
    --------
    df: DataFrame
        Filter DataFrame according to the filters located in the json_path.
    """
    with open(json_path) as f:
        probs_dict = json.load(f)

    dfs = []
    for station,probs in probs_dict.items():
        station_df = df[ (df["uncertainty"] >= probs[f'{phase_hint.upper()}'] ) & (df['station_code']==station)]
        dfs.append(station_df)
    
    df = pd.concat(dfs)
    return df

def phasenetDF2Pick(df,min_Pprob=None,min_Sprob=None,no_stations=[]):
    """
    Parameters:
    -----------
    df: DataFrame
        AIPicker-PhaseNet output csv file

    Returns:
    --------
    """
    Pdf = df[df["phase_hint"]=="P"]
    Sdf = df[df["phase_hint"]=="S"]


    ###### filter probability ######df,phase_hint,json_path
    json_path = "/home/ecastillo/tesis/catalog/aipicker/filters/phasenet_filter_prob.json"
    Pdf = filter_phasenet_probs(Pdf,'p',json_path)
    Sdf = filter_phasenet_probs(Sdf,'s',json_path)


    if min_Pprob != None:
        Pdf = Pdf[ (Pdf["uncertainty"] >= min_Pprob) ]
    if min_Sprob != None:
        Sdf = Sdf[ (Sdf["uncertainty"] >= min_Sprob) ]
        # print(min_Sprob)
        # print(Sdf)
    df = pd.concat([Pdf,Sdf])

    df["station_code"] = df["station_code"].map(lambda x: x.strip())
    if no_stations:
        df = df [ ~df["station_code"].isin(no_stations)]

    df = df.sort_values(by="time",ascending=True,ignore_index=True)
    picklist = []
    for index,row in df.iterrows():
        picktime = dt.datetime.fromisoformat(row["time"])
        p = Pick(publicID=row["resource_id"],
                pick_time=UTCDateTime(picktime),
                net=row["network_code"],
                station=row["station_code"],
                loc="{:02d}".format(row["location_code"]) ,
                ch=row["channel_code"],
                prob=row["uncertainty"],
                phaseHint=row["phase_hint"],
                creation_time=row["creation_time"],
                evaluation=row["evaluation_mode"],
                author=row["author"],
                snr=None,
                ev_prob =None,
                )
        picklist.append(p)
    return picklist

def eqtDF2Pick(df,min_Pprob=None,min_Sprob=None,no_stations=[]):
    """
    Parameters:
    -----------
    df: DataFrame
        AIPicker-EQTransformer output csv file

    Returns:
    --------
    """

    if min_Pprob != None:
        df = df[ (df["p_probability"] >= min_Pprob) | (df["p_probability"].isnull() == True) ]
    if min_Sprob != None:
        df = df[ (df["s_probability"] >= min_Sprob)  | (df["s_probability"].isnull() == True) ]

    df["station"] = df["station"].map(lambda x: x.strip())
    if no_stations:
        df = df [ ~df["station"].isin(no_stations)]
        
    df = df.sort_values(by="event_start_time",ascending=True,ignore_index=True)

    picklist = prepare_eqt(df)
    return picklist

def write_xml(pick_list,output_file=None):
    """
    Parameters:
    -----------
    pick_list: list
        List of Pick objects from utils
    output_file: str
        Path of the xml file

    """

    if os.path.isdir(os.path.dirname(output_file)) == False:
        os.makedirs(os.path.dirname(output_file))
    xml_text = picks2xml(pick_list)
    # Writting in the output file
    with open(output_file, 'w') as f:
        f.write(xml_text)
    
    print(f'\nOutput file: {output_file}')

def get_path(file_path,output_storage, new_filename):  
    """
    Parameters:
    -----------
    file_path: str
        julianday path
    output_storage: str
        path where will be located the new file
    new_filename: str
        name filename

    Returns:
    --------
    Path of the new filename

    example: >> file_path: /home/ecastillo/tesis/picksOK/CM/2020/001/eqt_picks.csv
                output_storage: /home/ecastillo/tesis/xml
                new_filename: eqt_picks.xml
            
            R>> output_path: /home/ecastillo/tesis/xml/CM/2020/001/eqt_picks.xml

    """
    f = os.path.basename(file_path)
    dp =os.path.dirname(file_path)

    julday = os.path.basename(dp)
    year_path = os.path.dirname(dp)
    year = os.path.basename(year_path)
    net = os.path.basename(os.path.dirname(year_path))

    # xml_name = os.path.splitext(f)[0]
    output_file = os.path.join(output_storage,net,year,julday, new_filename)
    return  output_file

def get_PickObjects(csv_file,ai_picker,min_Pprob=None,min_Sprob=None,no_stations=[]):
    """
    Parameters:
    -----------
    csv_file: str
        Path of the picker output
    ai_picker: str
        'eqt' or 'phasenet' or 'both'
    min_Pprob: float
        Minimum probability to select P phases.
    min_Sprob: float
        Minimum probability to select S phases.

    Returns:
    --------
    picklist : list
        List of pick objects
    """
    df = pd.read_csv(csv_file)
    if ai_picker in ("eqt","EQT","EQTransformer"):
        picklist = eqtDF2Pick(df,min_Pprob,min_Sprob,no_stations)
    elif ai_picker in ("pnet","phasenet","PhaseNet"):
        picklist = phasenetDF2Pick(df,min_Pprob,min_Sprob,no_stations)
    return picklist

def get_xml_Picks(pickobjects_list,output_file):
    """
    Write pickobjects in a picks xml file.

    Parameters:
    -----------
    pickobjects_list: list
        List of PickObjects
    output_file: str
        Path where will be located the xml picks file.

    """
    write_xml(pickobjects_list,output_file)

def get_xml_Origins(xml_picks,output_file,locator_type="LOCSAT",
                    locator_profile="iasp91",
                    db="sysop:sysopp@10.100.100.13/seiscomp3"):
    """
    Write origins in an origins xml file.

    Parameters:
    -----------
    xml_picks: str
        path where is lcoated the picks xml file
    output_file: str
        path where will be lcoated the origins xml file
    locator_type: str
        "LOCSAT" or "Hypo71"
    locator_profile: str
        "iaspei91" for LOCSAT and "RSNC" for Hypo71
    db: str
        "sysop:sysopp@10.100.100.13/seiscomp3" 
    run: Bolean
        True to run the scanloc message automatically


    Returns:
    --------
    origin xml file
    """

    if os.path.isdir(os.path.dirname(output_file)) == False:
        os.makedirs(os.path.dirname(output_file))

    scanloc_cmd = f'scanloc -u playback --locator-type {locator_type} --locator-profile {locator_profile} '
    scanloc_cmd += '--ep %s -d %s  > %s'%(xml_picks, db, output_file)

    
    if  locator_type == "Hypo71":

        tmp_file = os.path.join(os.path.dirname(output_file),
                                os.path.basename(output_file).split('.')[0]+'_tmp.xml')
        mv_msg = f"mv {output_file} {tmp_file} "
        # print(mv_msg)
        # os.system(mv_msg)

        key_word = '<?xml version="1.0" encoding="UTF-8"?>'
        clean_msg = f"sed -n '/{key_word}/,$p' {tmp_file} > {output_file}"
        # print(clean_msg)
        # os.system(clean_msg)

        rm_msg = f"rm {tmp_file}"
        # print(rm_msg)

        scanloc_cmd = ';'.join([scanloc_cmd,mv_msg,clean_msg,rm_msg])
        # os.system(scanloc_cmd)
    return scanloc_cmd

def get_xml_Amplitudes(xml_origins,output_file,
                db="sysop:sysopp@10.100.100.13/seiscomp3"):
    """
    Write amplitudes in an amplitudes xml file.

    Parameters:
    -----------
    xml_origins: str
        path where is located the origins xml file
    output_file: str
        path where will be lcocated the amplitudes xml file
    db: str
        "sysop:sysopp@10.100.100.13/seiscomp3" 

    """

    if os.path.isdir(os.path.dirname(output_file)) == False:
        os.makedirs(os.path.dirname(output_file))

    scamp_cmd = 'scamp -u playback --ep %s -d %s  > %s'%(xml_origins,
                                                            db, output_file)
    
    return scamp_cmd

def get_xml_Magnitudes(xml_amplitudes,output_file,
                db="sysop:sysopp@10.100.100.13/seiscomp3"):
    """
    Write magnitudes in an magnitudes xml file.

    Parameters:
    -----------
    xml_amplitudes: str
        path where is located amplitudes xml file
    output_file: str
        path where will be located the magnitudes xml file
    db: str
        "sysop:sysopp@10.100.100.13/seiscomp3" 

    """

    if os.path.isdir(os.path.dirname(output_file)) == False:
        os.makedirs(os.path.dirname(output_file))

    scmag_cmd = 'scmag -u playback --ep %s -d %s  > %s'%(xml_amplitudes, db, output_file)
    
    return scmag_cmd

def get_xml_Events(xml_magnitudes,output_file,
                    db="sysop:sysopp@10.100.100.13/seiscomp3"):
    """
    Write events in an magnitudes xml file.

    Parameters:
    -----------
    xml_magnitudes: str
        path where is located magnitudes xml file
    output_file: str
        path where will be located the events xml file
    db: str
        "sysop:sysopp@10.100.100.13/seiscomp3" 

    """

    if os.path.isdir(os.path.dirname(output_file)) == False:
        os.makedirs(os.path.dirname(output_file))

    scevent_cmd = 'scevent -u playback --ep %s -d %s  > %s'%(xml_magnitudes, db,
                                                                output_file)
    
    return scevent_cmd

def write_pref_origin(seiscomp_file,pref_origin_file,version="0.9",):
    def change_xml_version(ev_file,new_version="0.9"):
        lines = open(ev_file).readlines()
        lines[1] = f'<seiscomp xmlns="http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/{new_version}" version="{new_version}">\n'
        with open(ev_file, 'w') as f:
            f.write(''.join(lines))

    if version not in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.10"]:
        change_xml_version(seiscomp_file,new_version="0.10")


    catalog = read_events(seiscomp_file,format = "SC3ML" )
    event_list = catalog.events 

    events = []
    for event in event_list:
        pref_origin = event.preferred_origin() 
        del event.origins
        event.origins = [pref_origin]
        events.append(event)
    catalog.events = events
    # print(type(catalog))
    catalog.write(pref_origin_file,"SC3ML")

def get_pick_counts(df):
    df["P_Count"] = 1
    df["S_Count"] = 1
    Pcols = ["id","pick_p"]
    Scols = ["id","pick_s"]
    P_df = df.groupby(Pcols).P_Count.count().reset_index()
    S_df = df.groupby(Scols).S_Count.count().reset_index()
    counts_df = pd.merge(P_df,S_df,on="id")
    counts_df = counts_df[["id","P_Count","S_Count"]]

    counts = counts_df.set_index('id').T.to_dict('list')
    counts_df = df['id'].map(counts)
    df ["P_Count"] = counts_df.apply(lambda x: x[0])
    df ["S_Count"] = counts_df.apply(lambda x: x[1])
    df = df.rename(columns={"P_Count":"#P_picks",
                        "S_Count":"#S_picks"})
    return df

def get_csv_events(seiscomp_file, version="0.9", with_magnitude=True, 
                    picker=None,sort='time_event',export = None, 
                    from_format="SC3ML",pick_counts=False, inside_polygon=False ):
    """
    parameters
    ----------
    seiscom_file : str
        path of the sc3ml seiscomp file that contains only one event
    version: str
        String of the xml version number.
    with_magnitude: Bolean
        If True return the magnitude and magnitude information
    picker: None
        In 'eqt' or 'phasenet' doesn't matter the channel where is located the pick.
        While in others matters the channel. None is select to  have importance in the channel
    export: str (deault : None)
        Path to export in a csv file. None don't create any file, then only returns.

    returns
    -------
    appended_events : list
        It's a list of Pandas DataFrames where each dataframe contains all information about the picks from the respective event
    """
    def change_xml_version(ev_file,new_version="0.9"):
        lines = open(ev_file).readlines()
        lines[1] = f'<seiscomp xmlns="http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/{new_version}" version="{new_version}">\n'
        with open(ev_file, 'w') as f:
            f.write(''.join(lines))

    if from_format == "SC3ML":
        if version not in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.10"]:
            change_xml_version(seiscomp_file,new_version="0.10")


    seiscomp_file_path = os.path.splitext(seiscomp_file)[0]
    catalog = read_events(seiscomp_file,format = from_format)
    event_list = catalog.events 

    events = []
    for event in event_list:
        loc_id = os.path.basename(str(event.resource_id))
        ev_type = event.event_type
        agency = event.creation_info.agency_id
        region = event.event_descriptions[0].text

        ## Preferred Origin
        pref_origin = event.preferred_origin() 
        time = pref_origin.time
        latitude = pref_origin.latitude
        latitude_error = pref_origin.latitude_errors.uncertainty
        longitude = pref_origin.longitude
        longitude_error = pref_origin.longitude_errors.uncertainty

        depth = pref_origin.depth
        depth_error = pref_origin.depth_errors
        rms = pref_origin.quality.standard_error
        method = os.path.basename(str(pref_origin.method_id))
        earth_model = os.path.basename(str(pref_origin.earth_model_id))
        evaluation_mode = pref_origin.evaluation_mode

        if depth != None:
            depth = float(depth)/1000 #in km
        if depth_error != None:
            if depth_error.uncertainty != None:
                depth_error = depth_error.uncertainty/1000 #in km
            else:
                depth_error = None
        else:
            depth_error = None
        ## Preferred Magnitude
        if with_magnitude:
            pref_magnitude = event.preferred_magnitude()
            # print("aca",pref_magnitude)
            if pref_magnitude != None:
                magnitude = pref_magnitude.mag
                magnitude_type = pref_magnitude.magnitude_type
                if magnitude != None:
                    magnitude = round(magnitude,2)
            else:
                magnitude = None
                magnitude_type = None

        else:
            magnitude = None
            magnitude_type = None
        
        ## Dictionary with the picks information
        picks = {}
        for pick in event.picks:
            if pick.resource_id.id not in picks.keys():
                
                ## conditions to get the filter_id
                if pick.filter_id == None:
                    prob = None
                    snr = None
                    ev_prob = None
                    pick.filter_id = ResourceIdentifier(id=None)
                else:
                    if "Probability" in str(pick.filter_id):
                        prob,snr,ev_prob = str(pick.filter_id).split("+")
                        prob = float(str(prob).split("_")[-1])
                        snr = str(snr).split("_")[-1]
                        ev_prob = str(ev_prob).split("_")[-1]
                        if snr != 'None':
                            snr = float(snr)
                        if ev_prob != 'None':
                            ev_prob = float(ev_prob)
                        pick.filter_id= ResourceIdentifier(id=None)
                    else:
                        prob = None
                        snr = None
                        ev_prob = None

                if pick.creation_info == None:
                    _author = None
                else:
                    _author = pick.creation_info.author

                picks[pick.resource_id.id] = {"network_code":pick.waveform_id.network_code,
                                            "station_code":pick.waveform_id.station_code,
                                            "location_code":pick.waveform_id.location_code,
                                            "channel_code":pick.waveform_id.channel_code,
                                            "phase_hint":pick.phase_hint,
                                            "time":pick.time,
                                            "author":_author,
                                            "probability": prob,
                                            "snr":snr,
                                            "ev_prob":ev_prob,
                                            "time_errors":pick.time_errors,
                                            "filter_id":pick.filter_id,
                                            "method_id":pick.method_id,
                                            "polarity":pick.polarity,
                                            "evaluation_mode":pick.evaluation_mode,
                                            "evaluation_status":pick.evaluation_status } 
        
        ## Preferred picks
        columns = ["agency","id","time_event","latitude","latitude_uncertainty",
                    "longitude","longitude_uncertainty","depth","depth_uncertainty",
                    "rms","region","method","earth_model","event_type","magnitude",
                    "magnitude_type","picker","network","station","location","channel",
                    "pick","time_pick","probability","snr","detection_probability"]
        data = []
        fmt = "%Y-%m-%d %H:%M:%S.%f" #2020-01-01 00:02:38
        for i,arrival in enumerate(pref_origin.arrivals):

            pick = picks[arrival.pick_id.id]

            line = [agency,loc_id,time.strftime(fmt),latitude,latitude_error,\
                    longitude,longitude_error,depth,depth_error,rms,region,\
                    method, earth_model,ev_type,  magnitude, magnitude_type,\
                    pick["author"],pick["network_code"],pick["station_code"],\
                    pick["location_code"],pick["channel_code"],pick["phase_hint"],\
                    pick["time"].strftime(fmt),pick["probability"],pick["snr"],pick["ev_prob"]]
            data.append(line)

        picks_df = pd.DataFrame(data,columns=columns)
        picks_p = picks_df[ picks_df["pick"] == 'P']
        picks_s = picks_df[ picks_df["pick"] == 'S']
        if picker in ("eqt","EQTransformer","eqtransformer",
                    "phasenet","PhaseNet"):
            merge_on = columns[:-5]
        else:
            merge_on = columns[:-6]
        picks_df = pd.merge(picks_p,picks_s,how='left',suffixes=("_p", "_s"),on=merge_on)
        events.append(picks_df)

    events_df = pd.concat(events)

    events_df.dropna(subset=['latitude','longitude'],inplace=True)

    events_df.index+=1
    events_df.index.name= 'No'

    if sort != None:   
        events_df = events_df.sort_values(by=sort,ascending=True,ignore_index=True)
    
    if pick_counts:
        events_df = get_pick_counts(events_df)

    if export != None:
        if os.path.isdir(os.path.dirname(export)) == False:
            os.makedirs(os.path.dirname(export))
        events_df.to_csv(export)
        print(f"Events_csv_file: {export}")
    return events_df

def merge_xml(xmls,output_file):
    """
    Parameters:
    -----------
    xmls: list
        List of paths of the xml files
    output_file: str
        Path of the merged xml
    """
    if len(xmls) <2 :
        raise Exception("Need two or more xml files to merge")

    xmls = " ".join(xmls)

    msg = f"scxmlmerge {xmls} > {output_file}"

    # os.system(msg)
    return msg

class SeiscompAssociator(object):
    def __init__(self,picks_csv_file,ai_picker,picks_xml_file,origins_file,
                events_xml_file,events_csv_file,
                amplitudes_file=False,magnitudes_file=False,
                origins_loc={"LOCSAT":"iasp91"},
                db="sysop:sysopp@10.100.100.13/seiscomp3",
                min_Pprob=None,min_Sprob=None,no_stations=[]):

        self.picks_csv_file = picks_csv_file
        self.picks_xml_file = picks_xml_file
        self.origins_file = origins_file
        self.amplitudes_file = amplitudes_file
        self.magnitudes_file = magnitudes_file
        self.events_xml_file = events_xml_file
        self.events_csv_file = events_csv_file
        self.location_types = list(origins_loc.keys())
        self.location_profiles = list(origins_loc.values())
        self.db = db
        self.min_Pprob = min_Pprob
        self.min_Sprob = min_Sprob
        self.no_stations = no_stations

        if ai_picker in ("eqt","EQT","EQTransformer"):
            self.ai_picker = "eqt"
        elif ai_picker in ("pnet","phasenet","PhaseNet"):
            self.ai_picker = "phasenet"
        elif ai_picker in ("sgc","SGC"):
            self.ai_picker = "sgc"

    def run_2(self):
        tic = time.time()
        picklist = get_PickObjects(self.picks_csv_file,self.ai_picker,self.min_Pprob,self.min_Sprob,self.no_stations)
        get_xml_Picks(picklist,self.picks_xml_file)
        

        def get_xml_path(xml_path,location_type,location_profile):
            name = os.path.basename(xml_path).split('.')[0]
            extra_name = f"_{location_type}_{location_profile}.xml"
            xml_path = os.path.join(os.path.dirname(xml_path),name + extra_name)
            return  xml_path

        msgs = []
        events = []
        for i in range(0,len(self.location_types)):
            origin_path = get_xml_path(self.origins_file,self.location_types[i],
                                        self.location_profiles[i])
            amplitude_path = get_xml_path(self.amplitudes_file,self.location_types[i],
                                        self.location_profiles[i])
            magnitude_path = get_xml_path(self.magnitudes_file,self.location_types[i],
                                        self.location_profiles[i]) 
            event_path = get_xml_path(self.events_xml_file,self.location_types[i],
                                        self.location_profiles[i]) 

            msg_orig = get_xml_Origins(self.picks_xml_file,origin_path,
                                            self.location_types[i],self.location_profiles[i],
                                            self.db)

            if self.amplitudes_file and self.magnitudes_file:
                msg_amp = get_xml_Amplitudes(origin_path,amplitude_path,self.db)
                msg_mag =get_xml_Magnitudes(amplitude_path,magnitude_path,self.db)
                msg_ev =get_xml_Events(magnitude_path,event_path,self.db)

                with_magnitude = True
                msg = ';'.join([msg_orig,msg_amp,msg_mag,msg_ev])
                # msg = ';'.join([msg_amp,msg_mag,msg_ev])

            else:
                msg_ev = get_xml_Events(origin_path,event_path,self.db)
                with_magnitude = False
                msg = ';'.join([msg_orig,msg_ev])

            msgs.append(msg)
            events.append(event_path)

        msg = ";".join(msgs)
        msg += ";"+ merge_xml(events,self.events_xml_file)

        print(msg + "\n")
        os.system(msg)

        get_csv_events(self.events_xml_file,version=0.11,with_magnitude=with_magnitude,
                        picker=self.ai_picker,export=self.events_csv_file)
        toc = time.time()
        print("{0:>15}".format(f'total time: {toc-tic:.2f}s'))

    def run(self):
        tic = time.time()
        # picklist = get_PickObjects(self.picks_csv_file,self.ai_picker,self.min_Pprob,self.min_Sprob,self.no_stations)
        # get_xml_Picks(picklist,self.picks_xml_file)
        
        filenames = []
        for i in range(0,len(self.location_types)):
            name = os.path.basename(self.origins_file).split('.')[0]
            extra_name = f"_{self.location_types[i]}_{self.location_profiles[i]}.xml"
            filename = os.path.join(os.path.dirname(self.origins_file),name + extra_name)
            merge_filename = os.path.join(os.path.dirname(self.origins_file), f"merge_{i}.xml")

            msg_orig = get_xml_Origins(self.picks_xml_file,filename,
                                            self.location_types[i],self.location_profiles[i],
                                            self.db)
            print(msg_orig)

            if i>0:
                msg_orig += ";" + merge_xml([filenames[i-1],filename],merge_filename)
                filename = merge_filename
            else:
                first_orig = msg_orig

            filenames.append(filename)

        msg_orig = first_orig +";"+ msg_orig

        if len(self.location_types) == 1:
            msg_orig += ";" + f"mv {filenames[0]} {self.origins_file}"
        else:
            msg_orig += ";" + f"mv {filenames[i]} {self.origins_file}"

        if self.amplitudes_file and self.magnitudes_file:
            msg_amp = get_xml_Amplitudes(self.origins_file,self.amplitudes_file,self.db)
            msg_mag =get_xml_Magnitudes(self.amplitudes_file,self.magnitudes_file,self.db)
            msg_ev =get_xml_Events(self.magnitudes_file,self.events_xml_file,self.db)

            with_magnitude = True
            msg = ';'.join([msg_orig,msg_amp,msg_mag,msg_ev])
            # msg = ';'.join([msg_amp,msg_mag,msg_ev])

        else:
            msg_ev = get_xml_Events(self.origins_file,self.events_xml_file,self.db)
            with_magnitude = False
            msg = ';'.join([msg_orig,msg_ev])

        print(msg + "\n")
        os.system(msg)

        get_csv_events(self.events_xml_file,version=0.11,with_magnitude=with_magnitude,
                        picker=self.ai_picker,export=self.events_csv_file)
        toc = time.time()
        print("{0:>15}".format(f'total time: {toc-tic:.2f}s'))


if __name__ == "__main__":

    ######################### write pref origin ###################
    # event = "/home/ecastillo/e1.xml"
    # pref_event = "/home/ecastillo/E1.xml"
    # write_pref_origin(event,pref_event,version=0.11)

    ################################ adding locator to origins ##########
    # xml_picks = "/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/xml/CM/2019/335/eqt_picks.xml"
    # output_file = "/home/ecastillo/o1.xml"
    # locator_type = "Hypo71"
    # locator_profile = "RSNC"
    # get_xml_Origins(xml_picks,output_file,locator_type,locator_profile)


    ############# sum xml ##############33
    # xml = '/home/ecastillo/eqt_origins.xml'
    # sum_origins(xml,xml, version="0.11")


    ### prove picks each 2 hours ###

    # picklist = get_PickObjects("/home/ecastillo/tesis/catalog/aipicker/picks/picks_2h/CM/2020/002/02__04/phasenet_stapicks_df.csv",
    #                         "phasenet")
    # get_xml_Picks(picklist,"/home/ecastillo/picks.xml")
    ############


    ################# filter phasenet probs ###############################
    # df = pd.read_csv("/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/csv/CM/2019/335/phasenet_picks.csv")
    # json_path =  "/home/ecastillo/tesis/catalog/aipicker/filters/phasenet_filter_prob_.json"
    # filter_phasenet_probs(df,'p',json_path)

    ######################################################################



    # get_julianpicks("/home/ecastillo/picksCM_2020-01",
    #                 "/home/ecastillo/picksCM_2020-01_cleaned",
    #                 "pnet")

    # store = "/home/ecastillo/tesis/picksOK"
    # out_store = "/home/ecastillo/tesis/picksxml"
    # get_xml_picks(store,out_store,"pnet",None,None)
    
    # xml_storage = "/home/ecastillo/tesis/picksxml"
    # output_storage = "/home/ecastillo/tesis/origins"
    # psa = PrepareSeiscompAssociator(xml_storage,output_storage,"eqt","sysop:sysopp@10.100.100.13/seiscomp3")
    # psa.get_origins("/home/ecastillo/tesis/picksxml/CM/2020/001/eqt_picks.xml","eqt_origins.xml")

    
    # csv_file = "/home/ecastillo/test/catalog/aipicker/picks/picks_1d/csv/CM/2019/335/eqt_picks.csv"
    # xml_file = "/home/ecastillo/test/catalog/aipicker/picks/picks_1d/xml/CM/2019/335/eqt_picks.xml"
    # origins_file = "/home/ecastillo/test/catalog/aipicker/origins/origins_1d/CM/2019/335/eqt_origins.xml"
    # amplitudes_file = "/home/ecastillo/test/catalog/aipicker/amplitudes/amplitudes_1d/CM/2019/335/eqt_amplitudes.xml"
    # magnitudes_file = "/home/ecastillo/test/catalog/aipicker/magnitudes/magnitudes_1d/CM/2019/335/eqt_magnitudes.xml"
    # events_csv_file = "/home/ecastillo/test/catalog/aipicker/events/events_1d/csv/CM/2019/335/eqt_events.csv"
    # events_file = "/home/ecastillo/test/catalog/aipicker/events/events_1d/xml/CM/2019/335/eqt_events.xml"
    # sc_associator = SeiscompAssociator(csv_file,"eqt",xml_file,
    #                                 origins_file, events_file,events_csv_file,amplitudes_file,
    #                                 magnitudes_file,{"LOCSAT":"iasp91","Hypo71":"RSNC"},min_Pprob=0,min_Sprob=0)
    # sc_associator.run()


    year = "2020"
    picker = "sgc"
    network = "CM_auto"  ###LO TENGO PARA Q NO VUELVA A CREAR XML DE PICKS
    picksdir_path = "/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d"
    originsdir_path = "/home/ecastillo/tesis/catalog/aipicker/origins/origins_1d"
    amplitudesdir_path = "/home/ecastillo/tesis/catalog/aipicker/amplitudes/amplitudes_1d"
    magnitudesdir_path = "/home/ecastillo/tesis/catalog/aipicker/magnitudes/magnitudes_1d"
    eventsdir_path = "/home/ecastillo/tesis/catalog/aipicker/events/events_1d"
    no_stations = ["BOG","NIZA","MEDEC","PAS2","MAN1C"]

    picks_file = f"{picker}_picks.csv"
    network_path = os.path.join(picksdir_path,"csv",network,year)
    association_files = []
    for dp, dn, filenames in os.walk(network_path):
        for f in filenames:
            if f == picks_file :
                csv_file = os.path.join(dp, f)
                xml_file = get_path(csv_file,os.path.join(picksdir_path,"xml"),f"{picker}_picks.xml")
                origins_file = get_path(csv_file,originsdir_path,f"{picker}_origins.xml")
                amplitudes_file = get_path(csv_file,amplitudesdir_path,f"{picker}_amplitudes.xml")
                magnitudes_file = get_path(csv_file,magnitudesdir_path,f"{picker}_magnitudes.xml")
                events_xml_file = get_path(csv_file,os.path.join(eventsdir_path,"xml"),f"{picker}_events.xml")
                events_csv_file = get_path(csv_file,os.path.join(eventsdir_path,"csv"),f"{picker}_events.csv")
                
                association_dict = {"csv_file":csv_file,"xml_file":xml_file,"origins_file":origins_file,
                                "amplitudes_file":amplitudes_file,"magnitudes_file":magnitudes_file,
                                "events_xml_file":events_xml_file,"events_csv_file":events_csv_file}
                association_files.append(association_dict)


    def SCassociator(assodict):
        print(assodict["csv_file"])
        sc_associator = SeiscompAssociator(assodict["csv_file"],picker,assodict["xml_file"],
                                            assodict["origins_file"], assodict["events_xml_file"],
                                            assodict["events_csv_file"],assodict["amplitudes_file"],
                                            assodict["magnitudes_file"],{"LOCSAT":"iasp91","Hypo71":"RSNC"},
                                            "sysop:sysopp@10.100.100.13/seiscomp3",
                                            # "sysop:sysop@localhost/seiscomp3", ##95
                                            None,None,no_stations)
        sc_associator.run()
        print("END:",assodict["csv_file"])


    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        executor.map(SCassociator,association_files)







    # merge_csv("/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/csv/VMM","eqt","pick","event_start_time","/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/csv/VMM/VMM_eqt_picks__20160101__20200901.csv")
    # merge_csv("/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/csv/CARMA","eqt","pick","event_start_time","/home/ecastillo/tesis/catalog/aipicker/picks/picks_1d/csv/CARMA/CARMA_eqt_picks__20160101__20170301.csv")
    # merge_csv("/home/ecastillo/tesis/catalog/aipicker/events/events_1d/csv/VMM","eqt","event","time_event","/home/ecastillo/tesis/catalog/aipicker/events/events_1d/csv/VMM/VMM_eqt_events__20160101__20200901.csv")
    # merge_csv("/home/ecastillo/tesis/catalog/aipicker/events/events_1d/csv/CARMA","eqt","event","time_event","/home/ecastillo/tesis/catalog/aipicker/events/events_1d/csv/VMM/CARMA_eqt_events__20160101__20170301.csv")


    # xml_prove = "/home/ecastillo/tesis/events/CM/2020/001/phasenet_events.xml"
    # csv_events = get_csv_events( xml_prove,"0.11",True,"phasenet",'time_event',"/home/ecastillo/tesis/events/CM/2020/001/phasenet_events.csv")
    # write_pref_origin(xml_prove,"/home/ecastillo/tesis/events/CM/2020/001/pref_eqt_events.xml",version="0.9")

    # # print(csv_events)
    # get_csv_events(self.events_xml_file,version=0.11,with_magnitude=with_magnitude,
    #                     picker=self.ai_picker,export=self.events_csv_file)
