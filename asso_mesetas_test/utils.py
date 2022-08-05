import pandas as pd
import datetime as dt

# eqt = "/home/emmanuel/Tesis/manual/picks/SGCmanual_picks_20191201000000_20210101000000.csv"
eqt = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/csv/CM/CM_eqt_events__20191201__20210101.csv"
# df = pd.read_csv(eqt,header=1)
df = pd.read_csv(eqt)
df["time_event"] = pd.to_datetime(df["time_event"])
df = df[(df["time_event"] >= dt.datetime(2019,12,24,19)) & \
        (df["time_event"] <= dt.datetime(2019,12,24,23))]
df = df.reset_index(drop=True)
df = df.drop_duplicates(subset="id")
print(len(df))
