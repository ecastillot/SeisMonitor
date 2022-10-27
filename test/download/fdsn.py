# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-10-26 14:13:40
#  * @modify date 2022-10-26 14:13:40
#  * @desc [description]
#  */

import os
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

out_download_folder = "./out/download/fdsn"

carma_rest = WaveformRestrictions(network="YU",
                    station="GJ*,CS*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2017-12-24T00:00:00.000000Z"),
                    endtime=UTCDateTime("2017-12-24T06:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    filter_domain= [-83.101,-64.549,-2.229,14.945],
                    )
carma_client = FDSNClient(base_url="IRIS")
carma_provider = Provider(carma_client,carma_rest)
json_path = os.path.join(out_download_folder,"json","test.json")


mseed_storage = os.path.join(out_download_folder,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

md = MseedDownloader(providers=[carma_provider])
md.make_inv_and_json(json_path)
md.download(mseed_storage,
            chunklength_in_sec=3600,n_processor=16)