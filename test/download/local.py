# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-10-26 14:13:40
#  * @modify date 2022-10-26 14:13:40
#  * @desc [description]
#  */


import os
# import wget
import requests
# from mega import Mega
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

from zipfile import (
    BadZipFile,   
    ZipFile,
)
from io import BytesIO,StringIO
import requests




out_folder = "./out/download/local"
if not os.path.isdir(out_folder):
    os.makedirs(out_folder)


SeisMonitor_dataset = '738dbcaed1b6fbb8ef72dd4885ae448ab4f6bfda'
# SeisMonitor_dataset = 'https://mega.nz/file/bl4TETbb#IfaVZvYce3zN6Ek7LTwckdt991NOSSPqz789EfDcYuk'

# # response = wget.download(SeisMonitor_dataset, os.path.join(out_folder,"SeisMonitor_dataset"))
r =  requests.get(SeisMonitor_dataset, allow_redirects=True, stream = True) 
z = ZipFile(BytesIO(r.content))
z.extractall(out_folder)



# # out_download_folder = "./out/download/fdsn"
# out_download_folder = "./data/archive"

# carma_rest = WaveformRestrictions(network="YU",
#                     station="GJ*,CS*",
#                     location="*",
#                     channel="*",
#                     starttime=UTCDateTime("2017-12-24T00:00:00.000000Z"),
#                     endtime=UTCDateTime("2017-12-24T06:00:00.000000Z"),
#                     location_preferences=["","00","20","10","40"],
#                     channel_preferences=["HH","BH","EH","HN","HL"],
#                     filter_networks=[], 
#                     filter_stations=[],
#                     filter_domain= [-83.101,-64.549,-2.229,14.945],
#                     )
# carma_client = FDSNClient(base_url="IRIS")
# carma_provider = Provider(carma_client,carma_rest)
# json_path = os.path.join(out_download_folder,"json","test.json")

# my_local_fmt = os.path.join("{year}-{month:02d}", 
#                 "{year}-{month:02d}-{day:02d}", 
#                 "{network}.{station}.{location}.{channel}.{year}.{julday:03d}")


# mseed_storage = os.path.join(out_download_folder,my_local_fmt)
# # mseed_storage = os.path.join(out_download_folder,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

# md = MseedDownloader(providers=[carma_provider])
# md.make_inv_and_json(json_path)
# md.download(mseed_storage,
#             chunklength_in_sec=3600,n_processor=16)