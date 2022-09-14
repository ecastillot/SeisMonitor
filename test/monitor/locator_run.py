import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
import logging
# logging.basicConfig(level=logging.DEBUG,
#                    format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',
#                    datefmt='%m-%d %H:%M') 

catalog = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/select.out"
station_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM_201601_202206.xml"
vel_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"
grid_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/model/layer"
time_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/time/layer"
loc_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/loc/SeisMonitor"
p_control_file_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/p_test.in"
s_control_file_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/s_test.in"

vel_model = lut.VelModel(vel_path)
stations = lut.Stations(station_path)
nlloc = NLLoc(region = [-84,-62,-5,15,-5,200],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 1,
        tmp_folder="/home/emmanuel/EDCT/test2_nlloc"
        )
nlloc.download()
nlloc.compute_travel_times()
# nlloc_folder="/home/emmanuel/EDCT/test_nlloc/loc_o"
# nlloc.locate(catalog,nlloc_folder)
exit()