# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-10-26 14:13:40
#  * @modify date 2022-10-26 14:13:40
#  * @desc [description]
#  */

import os
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

out_folder = "./out/download/local"

carma_rest = WaveformRestrictions(network="YU",
                    station="GJ*",
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

archive = "/home/emmanuel/Descargas/SeisMonitor_dataset/archive/mseed"
my_local_fmt = os.path.join("{year}-{month:02d}", 
                "{year}-{month:02d}-{day:02d}", 
                "{network}.{station}.{location}.{channel}.{year}.{julday:03d}")
carma_client = LocalClient(archive,my_local_fmt)
# st = carma_client.get_waveforms(network="YU",
#                     # station="GJ*,CS*",
#                     station="GJ*",
#                     location="*",
#                     channel="*",
#                     starttime=UTCDateTime("2017-12-24T00:00:00.000000Z"),
#                     endtime=UTCDateTime("2017-12-24T00:10:00.000000Z"))
# print(st)

xml_path = "/home/emmanuel/Descargas/SeisMonitor_dataset/archive/dataless/YU.xml"
carma_provider = Provider(carma_client,carma_rest,xml=xml_path)

json_path = os.path.join(out_folder,"stations")



dld_fmt = "downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
mseed_storage = os.path.join(out_folder,dld_fmt)

md = MseedDownloader(providers=[carma_provider])
md.make_inv_and_json(json_path)
md.download(mseed_storage,
            chunklength_in_sec=3600,n_processor=16)