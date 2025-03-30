# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-10-26 14:13:40
#  * @modify date 2022-10-26 14:13:40
#  * @desc [description]
#  */
import sys
repository_path = r"/home/edc240000/SeisMonitor"  
sys.path.insert(0,repository_path)

import os
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

out_download_folder = "/home/edc240000/SeisMonitor/test/out/download/fdsn"

sgc_rest = WaveformRestrictions(network="YU",
                    station="CS*,FC*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2017-01-02T19:00:00.000000Z"),
                    endtime=UTCDateTime("2017-01-03T01:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    #filter_domain= [-83.101,-64.549,-2.229,14.945],
                    )
sgc_client = FDSNClient('IRIS')
sgc_provider = Provider(sgc_client,sgc_rest)

json_path = os.path.join(out_download_folder,"stations")


mseed_storage = os.path.join(out_download_folder,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

md = MseedDownloader(providers=[sgc_provider])
md.make_inv_and_json(json_path)
md.download(mseed_storage,picker_args={"batch_size":100,"overlap":0.3,"length":60},
            chunklength_in_sec=7200,n_processor=None)