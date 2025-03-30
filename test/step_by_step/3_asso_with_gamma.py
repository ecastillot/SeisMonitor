import sys
repository_path = r"/home/edc240000/SeisMonitor"  
sys.path.insert(0,repository_path)

import os
from SeisMonitor.monitor.associator import utils as ut
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj

# archive = "./out/download/fdsn"
archive = "/home/edc240000/SeisMonitor/test/others"

resp = os.path.join(archive,"stations","inv.xml")
picks = os.path.join(archive,"picks","eqt_seismonitor_picks.csv")
out_dir = os.path.join(archive,"asso","gamma")

# region = [-84.798, -66.546,-1.628, 15.445,0, 150]
region = [-76.729, -72.315,1.55, 5.314,0, 150]

gc = GaMMAObj(region,"EPSG:3116",
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False)
g = GaMMA(gc)
obspy_catalog, df_catalog,df_picks = g.associate(picks,resp,out_dir)
print(obspy_catalog)