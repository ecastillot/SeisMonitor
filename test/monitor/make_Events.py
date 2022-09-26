import os
import glob

import pathlib
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.catalog import read_events
from obspy.core.event import Catalog
import pandas as pd
import concurrent.futures as cf
def get_pick_counts(df):
    Pcols = ["id","pick_p"]
    Scols = ["id","pick_s"]
    P_df = df[Pcols]
    S_df = df[Scols]

    P_df = P_df.dropna()
    S_df = S_df.dropna()
    df ["#P_picks"] = len(P_df)
    df ["#S_picks"] = len(S_df)
    return df

    # df["P_Count"] = 1
    # df["S_Count"] = 1
    # Pcols = ["id","pick_p"]
    # Scols = ["id","pick_s"]
    # P_df = df.groupby(Pcols).P_Count.count().reset_index()
    # S_df = df.groupby(Scols).S_Count.count().reset_index()

    # if P_df.empty :
    #     df["P_Count"] = 0
    # if S_df.empty :
    #     df["P_Count"] = 0

    # print(P_df,S_df )
    # counts_df = pd.merge(P_df,S_df,on="id")
    # counts_df = counts_df[["id","P_Count","S_Count"]]
    # print(counts_df[["id","P_Count","S_Count"]])
    # counts = counts_df.set_index('id').T.to_dict('list')
    # counts_df = df['id'].map(counts)
    # print(counts_df)
    # df ["P_Count"] = counts_df.apply(lambda x: x[0])
    # df ["S_Count"] = counts_df.apply(lambda x: x[1])
    # df = df.rename(columns={"P_Count":"#P_picks",
    #                     "S_Count":"#S_picks"})
    # return df

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

    # events = []
    # for event in event_list:
    def make_event(event):
        # print(event)
        loc_id = os.path.basename(str(event.resource_id))
        ev_type = event.event_type
        agency = event.creation_info.agency_id
        # print(event.event_descriptions)
        if not event.event_descriptions:
            region = None
        else:
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
        # events.append(picks_df)
        return picks_df

    with cf.ThreadPoolExecutor() as executor:
        events = list(executor.map(make_event,event_list))

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
        events_df.to_csv(export,index=False)
        print(f"Events_csv_file: {export}")
    return events_df

def merge_csv(storage,sort=None,export=None):
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


    events = []
    for dp, dn, filenames in os.walk(storage):
        for f in filenames:
            if pathlib.Path(f).suffix ==".csv":
                search_path = os.path.join(dp, f)
                print(search_path)
                # df = pd.read_csv(search_path, index_col=0)
                df = pd.read_csv(search_path)
                # df["depth"] = df["depth"] *1e3
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





# file = "Ml_magnitude.xml"
# out_folder = r"/media/emmanuel/TOSHIBA EXT/Events"
# storage = r"/media/emmanuel/TOSHIBA EXT/ColSeismicity"
# print(os.path.join(storage,"**",file))
# x = glob.glob(os.path.join(storage,"**",file),recursive=True)


# def make_evs(path):

#     ev_path = path.split("/")[6]
#     ev_path = os.path.join(out_folder,ev_path+".csv")
#     get_csv_events(path,pick_counts=False,export=ev_path)

# for i in x:
#     make_evs(i)


# ### collect
# storage = r"/media/emmanuel/TOSHIBA EXT/Events"
# export="/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901.csv"
# merge_csv(storage,sort="time_event",export=export)


### good_events
# events = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901.csv"
# events_ok = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901_ok.csv"

# evs=[]
# events = pd.read_csv(events)
# gp_ev = events.groupby("id")
# for id,ev in gp_ev.__iter__():
#     # print(ev)
#     ev = get_pick_counts(ev)
#     # print("\n")
#     ev = ev[(ev["#P_picks"] >=4) & (ev["#S_picks"] >=2)]
#     if ev.empty:
#         continue
#     evs.append(ev)
# events_df = pd.concat(evs,ignore_index=True)
# events_df = events_df.sort_values(by="time_event",ascending=True,ignore_index=True)
# if os.path.isdir(os.path.dirname(events_ok)) == False:
#     os.makedirs(os.path.dirname(events_ok))
# events_df.to_csv(events_ok,index=False)
# print(f"Events_csv_file: {events_ok}")



###other
ev_path= "/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901_ok.csv"
events = pd.read_csv(ev_path)
print(events)
events["No"] = events.index
events.to_csv(ev_path,index=False)