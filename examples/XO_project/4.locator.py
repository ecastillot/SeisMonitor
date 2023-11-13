import os
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from SeisMonitor.monitor.locator.nlloc import utils as nlloc_utils
from obspy import read_events
from obspy.core.event.catalog import Catalog

nlloc_core_path = "/home/emmanuel/NLLoc"
monitor_path = "/home/emmanuel/XO_monitor_results"
vel_path = "/home/emmanuel/XO_monitor/vel_model.csv"


####### NO MODIFY THE FOLLOWING

out_dir = os.path.join(monitor_path,"loc","nlloc")
asso_dir = os.path.join(monitor_path,"gamma_asso","eqt")
gamma__catalog = os.path.join(asso_dir,"associations.xml")
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

eqt_nlloc_catalog = nlloc.locate(catalog=gamma__catalog,
                            nlloc_out_folder= out_dir,
                            out_filename = "nlloc_loc.xml",
                            out_format="SC3ML" )
print(eqt_nlloc_catalog)