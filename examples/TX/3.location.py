import os
import obsplus # To convert into pandas dataframe .to_df()
from pathlib import Path
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut

main_nlloc_path = "/opt/ohpc/pub/apps/nonlinloc"
vel_path = "/groups/igonin/ecastillo/SeisMonitor/examples/TX/sm/vel_model/DB_model.csv"

# Make sure the region covers all the stations, the expected earthquake locations and the velocity model.
# Specially the elevation part and the stations, 0 respect to sea level, 
# negative means below sea level, positive means above sea level.
region = [-104.17816, -103.80355, 31.48921, 31.72133,-2, 12] 
delta_in_km = 0.5 # grid spacing for nlloc, smaller means more precise but also more computationally expensive.

monitor_path = Path(__file__).parent / "sm2"
stations_path = monitor_path / "stations"
picks_path = monitor_path / "picks"
asso_path = monitor_path / "associations"
loc_path = monitor_path / "locations"

tt_loc_folder = monitor_path / "tt" #travel time folder for nlloc, can be anywhere but it could consume significant disk space.

####### default configuration

nlloc_output_path = loc_path / "nlloc"
xml_stations_path = monitor_path / "stations" / "inv.xml"
xml_asso_path = asso_path / "associations.xml"
xml_nlloc_path = nlloc_output_path / "nlloc_catalog.xml"

vel_model = lut.VelModel(str(vel_path),model_name="Velmodel1D")
stations = lut.Stations(str(xml_stations_path))

nlloc = NLLoc(
        core_path = str(main_nlloc_path), ### type your NLLoc path, 
        agency="SeisMonitor",
        region = region,
        vel_model = vel_model,
        stations = stations,
        delta_in_km = delta_in_km,
        tmp_folder=str(tt_loc_folder) 
        )

# Use this to compute travel times only once, and then you can reuse the travel time files for future locations.
nlloc.compute_travel_times()
print(nlloc.tmp_folder)
print("nlloc.tmp_folder dir: ",nlloc.tmp_folder)
print("Important folders in nlloc.tmp_folder directory",os.listdir(nlloc.tmp_folder))

eqt_nlloc_catalog = nlloc.locate(catalog=str(xml_asso_path),
                            nlloc_out_folder= str(nlloc_output_path),
                            out_filename = str(xml_nlloc_path.name),
                            out_format="QUAKEML" )
print("Catalog\n",eqt_nlloc_catalog )
print("Events\n",eqt_nlloc_catalog.to_df())
print("Picks\n",eqt_nlloc_catalog.picks_to_df())