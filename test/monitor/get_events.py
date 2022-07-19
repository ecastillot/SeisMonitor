import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.core import utils as ut


xml_file = "/home/emmanuel/EDCT/test/magnitude/Ml_magnitude.xml"
out= "/home/emmanuel/EDCT/test/csv_events/events.csv"
df = ut.get_csv_events(xml_file,from_format="SC3ML",export=out)
print(df)