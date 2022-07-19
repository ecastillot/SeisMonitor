import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.monitor.locator.hypoDD.core import HypoDD

xml_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
vel_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"
out = "/home/emmanuel/EDCT/test"

catalog = "/home/emmanuel/EDCT/test/magnitude/Ml_magnitude.xml"
hyp = HypoDD(catalog,xml_path,vel_path,out)
hyp.locate()