import os 
import math
from obspy import read_events
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event.origin import OriginUncertainty
from obspy.core.event.event import Event
from SeisMonitor.monitor.seismonitor import SeisMonitor
from SeisMonitor.monitor.picker import ai as ai_picker
from SeisMonitor.monitor.downloader import utils as dut
from SeisMonitor.monitor.associator import ai as ai_asso
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.utils import get_merged_inv_and_json

out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/test"

sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-25T01:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    filter_domain= [-83.101,-64.549,-2.229,14.945],
                    )
sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_provider = Provider(sgc_client,sgc_rest)

seismo = SeisMonitor(providers = [sgc_provider],
                    chunklength_in_sec=7200,
                    out_folder = out)

seismo.add_downloader(picker_args= {"batch_size":100,"overlap":0.3,"length":60})

dataset = os.path.join(os.path.dirname(os.path.dirname(__file__)),"data")
eqt_model = os.path.join(dataset,"models",'EqT_model.h5')
seismo.add_picker(
                  pickers={
                            "EQTransformer":ai_picker.EQTransformerObj(
                                            model_path = eqt_model,
                                            n_processor = 32,
                                            overlap = 0.3,
                                            detection_threshold =0.1,
                                            P_threshold = 0.01,
                                            S_threshold = 0.01,
                                            batch_size = 100,
                                            number_of_plots = 0,
                                            plot_mode = 1,
                                            rm_downloads=True ) 
                            }
                )
seismo.add_associator(input=["EQTransformer"],
                        associators={
                        "GaMMA":ai_asso.GaMMAObj(
                                            [-85, -68,-2, 15,0, 180],
                                            "EPSG:3116",
                                            use_amplitude = False,
                                            use_dbscan=False,
                                            max_sigma11=5.0,
                                            calculate_amp=False)
                        }
                        )

dataset = os.path.join(os.path.dirname(os.path.dirname(__file__)),"data")
vel_path = os.path.join(dataset,"metadata","vel1d_col.csv")
vel_model = lut.VelModel(vel_path)
inv,_,_,_ = dut.get_merged_inv_and_json(seismo.providers)
stations = lut.Stations(inv)

nlloc = NLLoc(
        agency="SeisMonitor",
        region = [-85, -68,0, 15,-5, 205],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2.5,
        tmp_folder="/home/emmanuel/NLLoc_grid/NLLoc_grid" ### CHANGE PATH TO YOUR OWN PATH AND ALSO TAKE IN MIND THAT CONSUME DISK
        )

seismo.add_locator(input={"associations":("GaMMA","EQTransformer")},
                    locators={
                        "NLLOC":nlloc}
                        )


ml_params = {"a":1.019,"b":0.0016,"r_ref":140} #ojeda
k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
                    ml_params["b"]* (ml_params["r_ref"]-100) +3
Ml = lambda ampl,epi_dist : math.log10(ampl * 1e3) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
                                ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k

seismo.add_magnitude(input={"locations":("NLLOC","GaMMA/EQTransformer")},
                    magnitudes={
                        "Ml":{"mag_type":Ml,
                                "trimmedtime":5,
                                "out_format":"SC3ML"}}
                                )
seismo.run()