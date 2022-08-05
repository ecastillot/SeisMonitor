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

resp = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM_201601_202206.xml"
picks = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/test/seismonitor_picks.csv"
out = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/test"
# picks = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/picks_in_mesetas.csv"
# out = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/test2"

# region = [-84.798, -66.546,-1.628, 15.445,0, 41]
region = [-76.729, -72.315,1.55, 5.314,0, 150]

gc = GaMMAObj(resp,region,"EPSG:3116",
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False,
                max_sigma11=20,
                max_sigma22=2,max_sigma12=2)
g = GaMMA(picks,resp,out)
catalog = g.associator(gc)
print(catalog)
for event in catalog:
    print(event)
print(len(catalog))