import os
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from SeisMonitor.monitor.locator.nlloc import utils as nlloc_utils
from obspy import read_events
from obspy.core.event.catalog import Catalog

dataset = os.path.join("/home/emmanuel/SeisMonitor","data")

archive = "/home/emmanuel/SeisMonitor/out"
out_dir = os.path.join(archive,"loc","nlloc")
asso_dir = os.path.join(archive,"asso","gamma")
vel_path = os.path.join(dataset,"velmodel","vel1d_col.csv")
inv = os.path.join(dataset,"stations","inv.xml")

vel_model = lut.VelModel(vel_path,model_name="Ojeda&Havskov(2004)")
stations = lut.Stations(inv)

nlloc = NLLoc(
        agency="SeisMonitor",
        region = [-85, -68,0, 15,-5, 205],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2.5,
        tmp_folder="/home/emmanuel/NLLoc_grid/NLLoc_grid" ### CHANGE PATH TO YOUR OWN PATH AND ALSO TAKE IN MIND THAT CONSUME DISK
        )
nlloc.download()
# nlloc.compute_travel_times()
# nlloc.locate(catalog=os.path.join(asso_dir,"associations.xml"),
#              nlloc_out_folder= out_dir,
#               out_filename = "LOC.xml",
#               out_format="SC3ML" )