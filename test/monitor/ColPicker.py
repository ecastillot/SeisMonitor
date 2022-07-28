import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

import os
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut

# sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_client = FDSNClient('http://10.100.100.232:8091')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-24T20:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"

carma_client = FDSNClient(base_url="IRIS")
carma_rest = WaveformRestrictions(network="YU",
                    station="*",
                    location="*",
                    channel="H*",
                    starttime=UTCDateTime("2016-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2016-12-24T20:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )

# print(rest.domain)
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml)
carma_provider = Provider(carma_client,carma_rest)
md = MseedDownloader(providers=[sgc_provider])
json_path = "/home/emmanuel/EDCT/test/json/test.json"
mseed_storage = ("/home/emmanuel/EDCT/test/downloads/"
                  "{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
md.make_inv_and_json(json_path)
md.download(mseed_storage,chunklength_in_sec=3600,n_processor=None)

eqt_model = os.path.join("/home/emmanuel/test/seismo",'EQTransformer/ModelsAndSampleData/EqT_model.h5')
eqtobj = EQTransformerObj(model_path = eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 20,
            number_of_plots = 0,
            plot_mode = 1 )  
            
mseed_storage = "/home/emmanuel/EDCT/test/downloads"
json_path = "/home/emmanuel/EDCT/test/json/test.json"
out_dir = "/home/emmanuel/EDCT/test/picks/eqt"
eqt = EQTransformer(mseed_storage,json_path,out_dir)
eqt.picker(eqtobj)



eqt_folder = "/home/emmanuel/EDCT/test/picks/eqt"
result = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
piut.eqt_picks_2_seismonitor_fmt(eqt_folder,mseed_storage,result)