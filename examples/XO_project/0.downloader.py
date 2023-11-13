# /**
#  * @author Emmanuel Castillo
#  * @email ecastillot@unal.edu.co / castillo.280997@gmail.com
#  * @create date 2023-08-05 21:03:30
#  * @modify date 2023-08-05 21:03:30
#  * @desc [description]
#  */

import os
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

monitor_path = r"/home/emmanuel/XO_monitor_results"

sgc_rest = WaveformRestrictions(network="XO",
                    station="EP*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2018-07-01T00:00:00.000000Z"),
                    endtime=UTCDateTime("2018-07-01T04:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    )


####### NO MODIFY THE FOLLOWING


sgc_client = FDSNClient('IRIS')
sgc_provider = Provider(sgc_client,sgc_rest)
md = MseedDownloader(providers=[sgc_provider])

json_path = os.path.join(monitor_path,"stations")
inv,json = md.make_inv_and_json(json_path)

mseed_storage = os.path.join(monitor_path,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
md.download(mseed_storage,
            picker_args={"batch_size":100,"overlap":0.3,"length":60},
            chunklength_in_sec=7200,n_processor=None)