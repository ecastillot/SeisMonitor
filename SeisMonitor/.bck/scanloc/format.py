import pandas as pd
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.event import Event

date_FMT = "%Y%m%d %H:%M:%S.%f"

eqt_FMT = ["file_name","network","station","instrument_type",
            "station_lat","station_lon","station_elv",
            "event_start_time","event_end_time","detection_probability",
            "detection_uncertainty","p_arrival_time","p_probability",
            "p_uncertainty","p_snr","s_arrival_time","s_probability",
            "s_uncertainty","s_snr"]

eqt_DATE_FORMAT = ["event_start_time","event_end_time","p_arrival_time","s_arrival_time"]
P_eqt_FMT = [ x for x in eqt_FMT if x not in ["s_arrival_time","s_probability","s_uncertainty","s_snr"]]
S_eqt_FMT = [ x for x in eqt_FMT if x not in ["p_arrival_time","p_probability","p_uncertainty","p_snr"]]

def make_resource_id_column(picktime,network,
                            station,channel):

    picktime_fmt = picktime.strftime("%Y-%m-%d %H:%M:%S.%f")
    id_name = f"{picktime_fmt}.{network}.{station}.{channel}"

    prefix = f"SeisMonitor:pick"

    return ResourceIdentifier(id =id_name,prefix=prefix)
    
def split_eqt_phases(df,sortby="arrival_time"):
    p_df = df[P_eqt_FMT]
    s_df = df[S_eqt_FMT]

    dfs = []
    for phasehint,onedf in [("p",p_df),("s",s_df)]:
        keys = [f"{phasehint.lower()}_arrival_time",f"{phasehint.lower()}_probability",
                f"{phasehint.lower()}_uncertainty",f"{phasehint.lower()}_snr"]
        removing_phasehint = lambda x:x.split(f"{phasehint.lower()}_")[-1] 
        new_keys = list(map(removing_phasehint,keys))
        columns = dict(zip(keys, new_keys))
        onedf = onedf.rename(columns=columns)
        onedf = onedf.dropna(subset=["arrival_time"])
        onedf["phasehint"] = phasehint.upper()
        onedf["event_box_id"] = df["network"] + "_" +df["station"] +\
                                 "_" + df["event_start_time"].astype(str) + \
                                 "_" + df["event_end_time"].astype(str)
        dfs.append(onedf)

    df = pd.concat(dfs)
    df = df.sort_values(sortby,ignore_index=True)
    return df



def eqt2seismonitor(csv):
    df = pd.read_csv(csv)

    df[eqt_DATE_FORMAT] = df[eqt_DATE_FORMAT].apply(pd.to_datetime)
    df = split_eqt_phases(df)


    print(df)

    # df[eqt_DATE_FORMAT] = df[eqt_DATE_FORMAT].apply(.dt.tz_localize)
    
    # print(df)

    # print(df)

if __name__ == "__main__":
    csv = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/CM_eqt_picks__20191201__20210101.csv"
    
    eqt2seismonitor(csv)
    # df = pd.read_csv(csv)

    # print(df)