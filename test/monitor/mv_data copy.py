import os 
import pandas as pd
import datetime as dt
from obspy import UTCDateTime
# startdate = dt.datetime(2022,12,1)
# enddate = dt.datetime(2022,12,1)
# storage = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/2019"
# search = "eqt_picks.csv"

SeisMonitor_columns = ["pick_id","arrival_time","probability","phasehint",
                    "network","station","location","instrument_type","author",
                    "creation_time","event_start_time","event_end_time",
                    "detection_probability","snr","station_lat","station_lon",
                    "station_elv","file_name"]

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

def eqt_picks_2_seismonitor_fmt(eqt_folder,out_path):
    #eqt
    date_cols = ["p_arrival_time","s_arrival_time",
                "event_start_time","event_end_time"]

    dfs = []
    for dp, dn, filenames in os.walk(eqt_folder):
        for f in filenames:
            if f == "eqt_picks.csv" :
                search_path = os.path.join(dp, f)
                print(search_path)
                df = pd.read_csv(search_path)
                if not df.empty:
                    df[date_cols] = df[date_cols].apply(pd.to_datetime)
                    dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs,ignore_index=True)

    # df = pd.read_csv("/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/CM_eqt_picks__20191201__20210101.csv")
    # df = pd.read_csv("/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/VMM/VMM_eqt_picks__20160101__20200901.csv")
    # df = pd.read_csv("/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CARMA/CARMA_eqt_picks__20160101__20170301.csv")
    df[date_cols] = df[date_cols].apply(pd.to_datetime)
    # df = df [df["event_start_time"] <= dt.datetime(2017,1,1)]
    df["station"] = df["station"].apply(lambda x: x.strip())
    df = df.sort_values(by="p_arrival_time",ignore_index=True)

    # x_results_path = os.path.join(os.path.dirname(out_path),
    #                             "X_prediction_results_merge.csv")
    # df.to_csv(x_results_path,index=False,
    #             date_format="%Y-%m-%d %H:%M:%S.%f")

    #seismonitor
    date_cols = ["arrival_time","creation_time",
                "event_start_time","event_end_time"]

    get_loc = lambda x: x.file_name.split(".")[2]
    df["location"] = df.apply(get_loc,axis=1)
    df["author"] = "EQTransformer"

    # filename_func = lambda x: os.path.join(mseed_folder,x)
    df["file_name"] = pd.Series(dtype='int')

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


def get_chunktimes(starttime,endtime,chunklength_in_sec, overlap_in_sec=0):
	"""
	Make a list that contains the chunktimes according to 
	chunklength_in_sec and overlap_in_sec parameters.

	Parameters:
	-----------
	starttime: obspy.UTCDateTime object
		Start time
	endtime: obspy.UTCDateTime object
		End time
	chunklength_in_sec: None or int
		The length of one chunk in seconds. 
		The time between starttime and endtime will be divided 
		into segments of chunklength_in_sec seconds.
	overlap_in_sec: None or int
		For more than one chunk, each segment will have overlapping seconds

	Returns:
	--------
	times: list
		List of tuples, each tuple has startime and endtime of one chunk.
	"""

	if chunklength_in_sec == 0:
		raise Exception("chunklength_in_sec must be different than 0")
	elif chunklength_in_sec == None:
		return [(starttime,endtime)]

	if overlap_in_sec == None:
		overlap_in_sec = 0

	deltat = starttime
	dtt = dt.timedelta(seconds=chunklength_in_sec)
	overlap_dt = dt.timedelta(seconds=overlap_in_sec)

	times = []
	while deltat < endtime:
		# chunklength can't be greater than (endtime-startime)
		if deltat + dtt > endtime:
			break
		else:
			times.append((deltat,deltat+dtt))
			deltat += dtt - overlap_dt

	if deltat < endtime:	
		times.append((deltat,endtime))
	# print(times)
	return times

## to seismonitor format
# eqt_folder = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM"
# out_path = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/eqt_picks_20191201_20210101.csv"
# eqt_picks_2_seismonitor_fmt(eqt_folder,out_path)

# print("begin")
# eqt_folder = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/VMM/2016"
# out_path = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/VMM/VMM_seismonitor_eqt_picks__20160101__20170101.csv"
# eqt_picks_2_seismonitor_fmt(eqt_folder,out_path)

# eqt_folder = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CARMA/2016"
# out_path = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CARMA/CARMA_seismonitor_eqt_picks__20160101__20170101.csv"
# eqt_picks_2_seismonitor_fmt(eqt_folder,out_path)
# exit()

# ################################################### mer VMM and YU 2016
# vmm = pd.read_csv("/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/VMM/VMM_seismonitor_eqt_picks__20160101__20170101.csv")
# yu = pd.read_csv("/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CARMA/CARMA_seismonitor_eqt_picks__20160101__20170101.csv")
# df = pd.concat([vmm,yu])
# df["arrival_time"] = pd.to_datetime(df["arrival_time"])
# df = df.sort_values("arrival_time")
# df.to_csv("/home/emmanuel/Tesis/seismonitor_eqt_picks__20160101__20170101.csv",index=False)
# exit()

############################################### to folders
# seismonitor_path = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/eqt_picks_20191201_20210101.csv"
# out_folder = "/home/emmanuel/Results"
# starttime = UTCDateTime("2019-12-01 00:00:00")
# endtime = UTCDateTime("2021-01-01 00:00:00")

seismonitor_path = "/home/emmanuel/CSE_082022/picks/seismonitor.csv"
out_folder = "/home/emmanuel/G-Ecopetrol/ecastillo/Avances/2022/Procesamiento_CSE/CSE_202208"
starttime = UTCDateTime("2022-08-01 00:00:00")
endtime = UTCDateTime("2022-09-01 00:00:00")

chunktimes = get_chunktimes(starttime,endtime,86400)
df = pd.read_csv(seismonitor_path)
df["arrival_time"] = pd.to_datetime(df["arrival_time"])

# print(df)
# exit()
for start,end in chunktimes:
    year = str(start.year)
    start_fmt, end_fmt = start.strftime("%Y%m%dT%H%M%S"),end.strftime("%Y%m%dT%H%M%S")
    folder_name = start_fmt+"__"+end_fmt
    detections = os.path.join(out_folder,folder_name,"detections","EQTransformer","results")
    csv_detections = os.path.join(detections,"seismonitor_picks.csv")

    if not os.path.isdir(detections):
        os.makedirs(detections)

    print(detections)
    time_df = df[(df["arrival_time"] >= start.datetime) & (df["arrival_time"] <= end.datetime)]
    time_df.to_csv(csv_detections,index=False)

    # print(time_df)


# events = []
# for dp, dn, filenames in os.walk(storage):
#     for f in filenames:
#         if f == search:
#             search_path = os.path.join(dp, f)
#             # print(search_path)
#             df = pd.read_csv(search_path, index_col=0)
#             events.append(df)

