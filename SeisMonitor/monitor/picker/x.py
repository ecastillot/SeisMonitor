

r = get_max_allowed_batch_size(120,60,0.3)
print(r)

print(int(2.8))
# import pandas as pd
# import numpy as np
# import concurrent.futures as cf
# picks = "/home/emmanuel/sss/20191224T190000__20191224T200000/detections/PhaseNet/results/picks.csv"


# df = pd.read_csv(picks)
# df = df.astype({'itp':'string','tp_prob':'string',
#                 'its':'string','ts_prob':'string'})
# df = df.replace("[]",np.nan)
# df = df.dropna(subset=['itp','tp_prob','its','ts_prob'],how="all")
# print(df)



# picks = []

# # for i,row in df.iterrows():
# def _get_picks(irow):
#     i,row = irow
#     wf_name = row["fname"]
#     picks_p = row["itp"].strip('[]').strip().split() 
#     prob_p = row["tp_prob"].strip('[]').strip().split() 
#     picks_s = row["its"].strip('[]').strip().split() 
#     prob_s = row["ts_prob"].strip('[]').strip().split() 
#     picks.append((wf_name,picks_p,picks_s,prob_s))

# with cf.ThreadPoolExecutor() as executor:
#     executor.map(_get_picks,df.iterrows())
# print(picks)

# print(df)
# print(df.info())
# df.replace("[]")
# df['itp'] = df['itp'].str.strip('[]').str.split(',')
# df['tp_prob'] = df['tp_prob'].str.strip('[]').str.split(',')
# df['its'] = df['its'].str.strip('[]').str.split(',')
# df['ts_prob'] = df['ts_prob'].str.strip('[]').str.split(',')
# # df [["itp","tp_prob","its","ts_prob"]] = df[["itp","tp_prob","its","ts_prob"]].apply(lambda x: eval(str(x)),axis=0)
# # df = df.explode(['itp','tp_prob','its','ts_prob'])
# print(df)
# print(df.info())