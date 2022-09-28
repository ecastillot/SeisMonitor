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
from SeisMonitor.monitor.downloader.utils import get_merged_inv_and_json
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
import os 
from obspy import read_events
from obspy.core.event.origin import OriginUncertainty
from obspy.core.event.event import Event

# out = "/home/emmanuel/inventories"
# out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2018"
# out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2020"
# out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2021"
out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2019"
# out = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2022"

sgc_client = FDSNClient('http://10.100.100.13:8091')
sgc2_client = FDSNClient('http://10.100.100.232:8091')
# # sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    # starttime=UTCDateTime("2020-01-05T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2020-03-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2020-03-06T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2020-05-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2020-05-19T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2020-07-01T00:00:00.000000Z"),

                    # starttime=UTCDateTime("2019-01-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2019-03-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2019-03-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2019-05-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2019-05-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2019-07-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2019-07-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2019-09-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2019-09-19T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2019-11-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2019-11-14T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2020-01-01T00:00:00.000000Z"),

                    # starttime=UTCDateTime("2021-01-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2021-03-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2021-03-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2021-05-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2021-05-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2021-07-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2021-07-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2021-09-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2021-09-19T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2021-11-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2021-11-14T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-01-01T00:00:00.000000Z"),

                    # starttime=UTCDateTime("2022-01-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-03-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2022-03-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-05-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2022-05-03T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-07-01T00:00:00.000000Z"),
                    # starttime=UTCDateTime("2022-07-29T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-09-01T00:00:00.000000Z"),


                #     starttime=UTCDateTime("2018-05-09T00:00:00.000000Z"),
                #     endtime=UTCDateTime("2019-01-01T00:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml )
sgc2_provider = Provider(sgc2_client,sgc_rest,xml=sgc_xml )

# carma_client = FDSNClient(base_url="IRIS")
# carma_rest = WaveformRestrictions(network="YU",
#                     station="*",
#                     location="*",
#                     channel="H*",
#                     starttime=UTCDateTime("2018-01-31T00:00:00.000000Z"),
#                     endtime=UTCDateTime("2018-01-01T00:00:00.000000Z"),
#                     location_preferences=["","00","20","10","40"],
#                     channel_preferences=["HH","BH","EH","HN","HL"],
#                     filter_networks=[], 
#                     filter_stations=[],
#                     # filter_domain= [-83.101,-64.549,-2.229,14.945],
#                     # filter_domain= [-76.536,-71.168,7.758,11.823],
#                     )
# carma_provider = Provider(carma_client,carma_rest)
# providers = [carma_provider,sgc_provider]
providers = [sgc_provider,sgc2_provider]

# inventory,json_info,updated_providers,stations_outside_domains = get_merged_inv_and_json(providers)
# inventory.write(r"/media/emmanuel/TOSHIBA EXT/ColSeismicity/NLLoc_grid/time_grid/CM_YU.xml",
#                 format="STATIONXML")  
# exit()
seismo = SeisMonitor(providers = providers,chunklength_in_sec=86400,
# seismo = SeisMonitor(providers = [sgc_provider,carma_provider],chunklength_in_sec=86400,
# # seismo = SeisMonitor(providers = [sgc_provider],
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
inv = "/home/emmanuel/NLLoc_grid/NLLoc_grid/time_grid/CM_YU.xml"
# inv,_,_,_ = dut.get_merged_inv_and_json(seismo.providers.copy())
stations = lut.Stations(inv)
# nlloc = NLLoc(region = [-84,-62,-5,15,-5,200],
nlloc = NLLoc(region = [-85, -68,0, 15,-5, 205],
# nlloc = NLLoc(region = [-81, -65,-3, 13,-5, 205],
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2.5,
        tmp_folder="/home/emmanuel/NLLoc_grid/NLLoc_grid"
        )
# nlloc.download()
# nlloc.compute_travel_times()
# nlloc.locate("/home/emmanuel/NLLoc_grid/NLLoc_grid/t.xml",
#               "/home/emmanuel/NLLoc_grid/NLLoc_grid/out" ,
#               "SC3ML" )

# catalog = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/xml/CM/2019/335/eqt_events.xml"
# catalog = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/xml/CM/2019/336/eqt_events.xml"
# catalog = read_events(catalog,format="SC3ML")
# for event in catalog:
#     origins = event.origins
#     for origin in origins:
#         origin_unc = origin.origin_uncertainty

#         if origin_unc == None:
#             origin_unc = OriginUncertainty(horizontal_uncertainty=0,
#                            min_horizontal_uncertainty=0,
#                             max_horizontal_uncertainty=0    )
#             origin.origin_uncertainty = origin_unc
# # pref origin
# # events = []
# # for event in catalog.events:
# #     pref_origin = event.preferred_origin() 
# #     origin_unc = pref_origin.origin_uncertainty
# #     pref_magnitude = event.preferred_magnitude() 

# #     if origin_unc == None:
# #         continue

# #     picks = {}
# #     for pick in event.picks:
# #         picks[pick.resource_id.id] = pick

# #     new_picks = []
# #     for i,arrival in enumerate(pref_origin.arrivals):
# #         pick = picks[arrival.pick_id.id]
# #         new_picks.append(pick)

# #     event.picks = new_picks
# #     ev = Event(origins = [pref_origin],magnitudes = [pref_magnitude],picks = new_picks)
# #     events.append(ev)
# # catalog.events = events

# for event in catalog.events:
#     origins = event.origins
#     print(len(origins))
        # print(origin_unc)
# print(catalog)
# nlloc.locate(catalog,"/home/emmanuel/1")

seismo.add_locator(input={"associations":("GaMMA","EQTransformer")},
                    locators={
                        "NLLOC":nlloc}
                        )
seismo.add_magnitude(input={"locations":("NLLOC","GaMMA/EQTransformer")},
                    magnitudes={
                        "Ml":{"mag_type":"RSNC",
                                "trimmedtime":5,
                                "out_format":"SC3ML"}}
                                )
seismo.run()