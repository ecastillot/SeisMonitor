import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

# test = "/home/emmanuel/EDCT/SeisMonitor/examples/test"
test = "/home/emmanuel/EDCT/SeisMonitor/examples/test"

from SeisMonitor.monitor.associator import utils as ut
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
import matplotlib.pyplot as plt
import os

picks = os.path.join(test,"seismonitor_picks.csv")
resp = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM_201601_202206.xml"

region = [-76.729, -72.315,1.55, 5.314,0, 150]
gc = GaMMAObj(resp,region,"EPSG:3116",
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False,
                oversample_factor=10,
                min_picks_per_eq=5,max_sigma11=10,
                max_sigma22=10,max_sigma12=20,
                vel = {"p": 7.0, 
                "s": 7.0 / 1.78})
g = GaMMA(picks,resp,test)
catalog = g.associator(gc)
print(catalog)
for event in catalog:
    print(event)