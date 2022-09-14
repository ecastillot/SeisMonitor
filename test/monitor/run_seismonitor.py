import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from SeisMonitor.monitor.seismonitor import SeisMonitor
from SeisMonitor.monitor.picker import ai as ai_picker
from SeisMonitor.monitor.downloader import utils as dut
from SeisMonitor.monitor.associator import ai as ai_asso
from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
from SeisMonitor.monitor.locator import utils as lut
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
import os 
out = "/home/emmanuel/Results/2019"

sgc_client = FDSNClient('http://10.100.100.13:8091')
# sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    # starttime=UTCDateTime("2021-06-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-01-01T00:00:00.000000Z"),
                    starttime=UTCDateTime("2019-12-01T00:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-02T00:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml )

# carma_client = FDSNClient(base_url="IRIS")
# carma_rest = WaveformRestrictions(network="YU",
#                     station="*",
#                     location="*",
#                     channel="H*",
#                     starttime=UTCDateTime("2016-01-01T00:00:00.000000Z"),
#                     endtime=UTCDateTime("2017-01-01T00:00:00.000000Z"),
#                     location_preferences=["","00","20","10","40"],
#                     channel_preferences=["HH","BH","EH","HN","HL"],
#                     filter_networks=[], 
#                     filter_stations=[],
#                     # filter_domain= [-83.101,-64.549,-2.229,14.945],
#                     # filter_domain= [-76.536,-71.168,7.758,11.823],
#                     )
# carma_provider = Provider(carma_client,carma_rest)

seismo = SeisMonitor(providers = [sgc_provider],chunklength_in_sec=86400,
# seismo = SeisMonitor(providers = [sgc_provider,carma_provider],chunklength_in_sec=86400,
# seismo = SeisMonitor(providers = [sgc_provider],
                    out_folder = out)
# seismo.add_downloader()
# seismo.add_picker(
#                   pickers={
#                             "EQTransformer":ai_picker.EQTransformerObj(
#                                             model_path = ai_picker.EQTransformer_model_path,
#                                             n_processor = 32,
#                                             overlap = 0.3,
#                                             detection_threshold =0.1,
#                                             P_threshold = 0.01,
#                                             S_threshold = 0.01,
#                                             batch_size = 100,
#                                             number_of_plots = 0,
#                                             plot_mode = 1,
#                                             rm_downloads=True ) ,
#                             # "PhaseNet":ai_picker.PhaseNetObj(model_path = ai_picker.PhaseNet_model_path,
#                             #                         mode='pred',
#                             #                         P_threshold=0.75,
#                             #                         S_threshold=0.75,
#                             #                         batch_size=100, 
#                             #                         one_single_sampling_rate=100,
#                             #                         plot=False, 
#                             #                         save_result=False,
#                             #                         rm_downloads=True) 
#                             }
#                 )
# print(seismo.process["download"])
# seismo.add_associator(input=["EQTransformer"],
#                         associators={
#                         "GaMMA":ai_asso.GaMMAObj(
#                                             [-85, -68,-2, 15,0, 180],
# #                                             [-76.729, -72.315,1.55, 5.314,0, 150],
#                                             "EPSG:3116",
#                                             use_amplitude = False,
#                                             use_dbscan=False,
#                                             max_sigma11=5.0,
#                                             calculate_amp=False)

#                         }
#                         )



vel_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"
vel_model = lut.VelModel(vel_path)
inv,_,_,_ = dut.get_merged_inv_and_json(seismo.providers.copy())
stations = lut.Stations(inv)
# nlloc = NLLoc(region = [-84,-62,-5,15,-5,200],
# nlloc = NLLoc(region = [-85, -68,0, 15,-5, 180],
nlloc = NLLoc(region = [-81, -65,-3, 13,-5, 180],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 3.5,
        tmp_folder="/home/emmanuel/EDCT/test2_nlloc"
        )
nlloc.download()
nlloc.compute_travel_times()

seismo.add_locator(input={"associations":("GaMMA","EQTransformer")},
                    locators={
                        "NLLOC":nlloc}
                        )
seismo.run()