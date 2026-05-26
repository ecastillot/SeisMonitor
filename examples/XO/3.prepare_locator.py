import os
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from SeisMonitor.monitor.locator.nlloc import utils as nlloc_utils
from obspy import read_events
from obspy.core.event.catalog import Catalog

monitor_path = "/home/emmanuel/XO_monitor_results"
nlloc_core_path = "/home/emmanuel/NLLoc"
vel_path = "/home/emmanuel/XO_monitor/vel_model.csv"

# download nlloc #this only works for ubuntu, if it doesn't work, you need to install NLLOC by yourself
nlloc_utils.download_nlloc(nlloc_core_path) ##ONLY ONCE, IT'S TO INSTALL NLLOC



####### NO MODIFY THE FOLLOWING

out_dir = os.path.join(monitor_path,"loc","nlloc")
asso_dir = os.path.join(monitor_path,"gamma_asso","eqt")
inv = os.path.join(monitor_path,"stations","inv.xml")
nlloc_grid = os.path.join(out_dir,"nlloc_grid")

vel_model = lut.VelModel(vel_path,model_name="vel_model")
stations = lut.Stations(inv)

nlloc = NLLoc(
        core_path = nlloc_core_path,
        agency="SeisMonitor",
        region = [-162, -150,50, 60,-2, 60],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2,
        tmp_folder=nlloc_grid ### TAKE IN MIND THAT IN THIS FOLDER YOU WILL DOWNLOAD YOUR TTs, SO IT CONSUMES A LOT OF SPACE IN YOUR DISK
        )

nlloc.compute_travel_times() ### ONLY ONCE, to compute the travel times, after the first run, this line it's not necessary