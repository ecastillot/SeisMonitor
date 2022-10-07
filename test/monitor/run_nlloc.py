import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from obspy import read_events
from obspy.core.event.catalog import Catalog
from SeisMonitor.monitor.locator.nlloc import utils as nlloc_utils
# catalog = read_events("/home/emmanuel/nlloc_iter/out/nlloc/SeisMonitor.20191225.144718.grid0.loc.hyp",
#                     format="NLLOC_HYP")
# catalog = nlloc_utils.write_pref_origin_removing_phaselocinfo(catalog)
# # for ev in catalog:
# #     picks = ev.picks
# #     print(picks)

# catalog.write("x.input",format="NORDIC")
# exit()

# catalog = read_events("/home/emmanuel/nlloc_iter/associations.xml")
# catalog.write("/home/emmanuel/nlloc_iter/events.xml",format="SC3ML")
# print(catalog)

# catalog = read_events("/home/emmanuel/nlloc_iter/events.xml")
# events = []
# for ev in catalog:
#     if ev.resource_id.id == "smi:local/20191225.144719.897000":
#         events.append(ev)
# catalog = Catalog(events)
# catalog.write("/home/emmanuel/nlloc_iter/2.xml",format="SC3ML")
# exit()

vel_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"
inv = "/home/emmanuel/NLLoc_grid/NLLoc_grid/time_grid/CM_YU.xml"
vel_model = lut.VelModel(vel_path,model_name="Ojeda&Havskov(2004)")
stations = lut.Stations(inv)

# print(stations.stations)
# exit()

nlloc = NLLoc(
        agency="SeisMonitor",
        region = [-85, -68,0, 15,-5, 205],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2.5,
        tmp_folder="/home/emmanuel/NLLoc_grid/NLLoc_grid"
        )
nlloc.iterlocate("/home/emmanuel/nlloc_iter/2.xml",
              "/home/emmanuel/nlloc_iter/out" ,
              out_filename = "LOC.xml",
              out_format="SC3ML" )
