import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.monitor.locator.nlloc.utils import *
import logging
# logging.basicConfig(level=logging.DEBUG,
#                    format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',
#                    datefmt='%m-%d %H:%M') 

catalog = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/select.out"
station_path = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/station.dat"
vel_path = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/model.dat"
grid_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/model/layer"
time_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/time/layer"
loc_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/loc/SeisMonitor"
p_control_file_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/p_test.in"
s_control_file_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/TEST4/s_test.in"

gen_control = GenericControlStatement(trans=["SIMPLE",
                                            5.0,-73.0,0.0])
vel2grid = Vel2Grid(vel_path,
                    grid_folder_out,
                    grid = [2,583,70,
                            -891.0,-891.0,-5.0,
                            3.0, 3.0, 3.0,
                            "SLOW_LEN"])
p_grid2time = Grid2Time(station_path,
                        grid_folder_out,
                        time_folder_out,
                        phase="P")
s_grid2time = Grid2Time(station_path,
                        grid_folder_out,
                        time_folder_out,
                        phase="S")
time2loc = Time2Loc(catalog=[catalog,"SEISAN"],
                    grid = [374,583,70,
                            -891.0,-891.0,-5.0,
                            3.0,3.0,3.0,
                            "PROB_DENSITY","SAVE"],
                    time_folder_out=time_folder_out,
                    loc_folder_out=loc_folder_out)
p_nlloc = NLLocControlFile(gen_control,vel2grid,
                p_grid2time,time2loc)
p_nlloc.write(p_control_file_out)
s_nlloc = NLLocControlFile(gen_control,vel2grid,
                s_grid2time,time2loc)
s_nlloc.write(s_control_file_out)