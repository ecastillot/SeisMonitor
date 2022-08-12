import pandas as pd
import os 
import datetime as dt

manual = "/home/emmanuel/Tesis/manual/picks/SGCmanual_picks_20191201000000_20210101000000.csv"
df = pd.read_csv(manual,header=1)
df["time_event"] = df["time_event"].apply(pd.to_datetime)
df = df[df["time_event"]>=dt.datetime(2019,12,24,19)]
df = df[df["time_event"]<=dt.datetime(2019,12,25,0)]
df = df.drop_duplicates(subset="id",ignore_index=True)
print(df)
print("manual",len(df))


dfs = []
folder = "/home/emmanuel/associations_result/mesetas_seccionado"
for dp, dn, filenames in os.walk(folder):
    for f in filenames:
        if f == "catalog.csv" :
            search_path = os.path.join(dp, f)
            df = pd.read_csv(search_path)
            if not df.empty:
                df["time"] = df["time"].apply(pd.to_datetime)
                dfs.append(df)
df = pd.concat(dfs,ignore_index=True)
df = df.sort_values(by="time",ignore_index=True)
print(df)
print("seccionado",len(df))


manual = "/home/emmanuel/associations_result/mesetas_all/20191224T190000__20191225T000000/associations/GaMMA2EQTransformer/associations/catalog.csv"
df = pd.read_csv(manual)
df["time"] = df["time"].apply(pd.to_datetime)
df = df[df["time"]>=dt.datetime(2019,12,24,19)]
df = df[df["time"]<=dt.datetime(2019,12,25,0)]
print(df)
print("todo_auto",len(df))

manual = "/home/emmanuel/associations_result/mesetas/20191224T190000__20191225T000000/associations/GaMMA2EQTransformer/associations/catalog.csv"
df = pd.read_csv(manual)
df["time"] = df["time"].apply(pd.to_datetime)
df = df[df["time"]>=dt.datetime(2019,12,24,19)]
df = df[df["time"]<=dt.datetime(2019,12,25,0)]
print(df)
print("5stations_auto",len(df))