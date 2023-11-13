import os
from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
from SeisMonitor.monitor.associator import utils as asut
import matplotlib.pyplot as plt

monitor_path = "/home/emmanuel/XO_monitor_results"

####### NO MODIFY THE FOLLOWING

region = [-162, -150,50, 60,0, 60]
gc = GaMMAObj(region,"EPSG:3338",
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False,
                method="BGMM",
                min_picks_per_eq=5,
                oversample_factor=1,
                max_sigma11=2.0,
                vel = {"p": 7.2, "s": 7.2/ 1.70})

inv = os.path.join(monitor_path,"stations","inv.xml")
picks = os.path.join(monitor_path,"picks","eqt","seismonitor_picks.csv")
out_dir = os.path.join(monitor_path,"gamma_asso","eqt")

g = GaMMA(gc)
obspy_catalog, df_catalog,df_picks = g.associate(picks,inv,out_dir)
print(obspy_catalog)


## just use this once to review the earthquake locations
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.set_aspect("equal")
ax.scatter(df_catalog["x(km)"], df_catalog["y(km)"])
ax.set_xlabel("x(km)")
ax.set_ylabel("y(km)")
plt.show()