import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.monitor.associator import utils as ut
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
import matplotlib.pyplot as plt
# eqt_f = "/home/emmanuel/EDCT/test/picks/picks/eqt_fmt"
# # ut.seismonitor_picks_to_eqt_fmt(sp,eqt_f,
# #                             tt=30,st=2,et=4)

picks = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
resp = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM_201601_202206.xml"
out = "/home/emmanuel/EDCT/test"

# region = [-84.798, -66.546,-1.628, 15.445,0, 41]
region = [-76.729, -72.315,1.55, 5.314,0, 150]

gc = GaMMAObj(resp,region,"EPSG:3116",
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False)
g = GaMMA(picks,resp,out)
catalog = g.associator(gc)
print(catalog)


exit()
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.set_aspect("equal")
cb = ax.scatter(catalog["x(km)"], catalog["y(km)"], c=catalog["z(km)"], s=8, cmap="viridis")
cbar = fig.colorbar(cb)
cbar.ax.set_ylim(cbar.ax.get_ylim()[::-1])
cbar.set_label("Depth[km]")

ax.plot(station_df["x(km)"], station_df["y(km)"], "r^", ms=10, mew=1, mec="k")
ax.set_xlabel("Easting [km]")
ax.set_ylabel("Northing [km]")
plt.show()
# print(x)