import os
import obsplus # To convert into pandas dataframe .to_df()
from pathlib import Path
from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
from SeisMonitor.monitor.associator import utils as asut
import matplotlib.pyplot as plt

monitor_path = Path(__file__).parent / "sm2"
stations_path = monitor_path / "stations"
picks_path = monitor_path / "picks"
asso_path = monitor_path / "associations"

# region = [min_lon, max_lon, min_lat, max_lat, min_depth_km, max_depth_km]
# EPSG:3081 -> UTM zone 14N, suitable for Texas
region = [-104.17816, -103.80355, 31.48921, 31.72133,0, 12] 
epsg = "EPSG:3081" 


####### default configuration
sm_picks_path = picks_path / "seismonitor_picks.csv"
inv_stations_path = monitor_path / "stations" / "inv.xml"

gc = GaMMAObj(region,epsg,
                use_amplitude = False,
                use_dbscan=False,
                calculate_amp=False,
                method="BGMM",
                min_picks_per_eq=5,
                oversample_factor=1,
                max_sigma11=2.0,
                vel = {"p": 6, "s": 6/ 1.70})


g = GaMMA(gc)
obspy_catalog, df_catalog,df_picks = g.associate(picks_csv=sm_picks_path,
                                    xml_path=inv_stations_path,
                                    out_dir=asso_path)
print("Catalog\n",obspy_catalog)
print("Events\n",obspy_catalog.to_df())
print("Picks\n",obspy_catalog.picks_to_df())
