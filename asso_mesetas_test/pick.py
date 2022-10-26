import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

test = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/test"

import os
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-24T23:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM_201601_202206.xml"
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml)

json_path = os.path.join(test,"json/test.json")
mseed_storage = os.path.join(test,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

# md = MseedDownloader(providers=[sgc_provider])
# md.make_inv_and_json(json_path)
# md.download(mseed_storage,chunklength_in_sec=3600,n_processor=None)

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
mseed_storage = os.path.join(test,"downloads")
json_path = os.path.join(test,"json","test.json")
out_dir = os.path.join(test,"picks","eqt")
result = os.path.join(test,"seismonitor_picks.csv")

eqt = EQTransformer(mseed_storage,json_path,out_dir)
eqt.pick(eqtobj)
piut.eqt_picks_2_seismonitor_fmt(out_dir,mseed_storage,result)