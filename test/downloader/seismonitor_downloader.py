import sys
Seismonitor_path = "/home/emmanuel/Ecopetrol/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from obspy.clients.fdsn import Client as FDSN_Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.utils import DownloadRestrictions
from SeisMonitor.downloader.seismonitor import MseedDownloader

client = FDSN_Client('http://sismo.sgc.gov.co:8080')
rest = DownloadRestrictions(network="CM",
                    station="AGCC,EZNC,BAR2,URMC,YOT,PDSC",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2020-09-22T00:00:00.000000Z"),
                    endtime=UTCDateTime("2020-09-23T00:00:00.000000Z"),
                    chunklength_in_sec=86400,
                    overlap_in_sec=None,
                    groupby='{network}.{station}.{location}',
                    threshold=60,
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    to_pick=(10,0.3))
ppc_rest = None
mseed_storage = ("/home/emmanuel/Ecopetrol/SeisMonitor/out/downloads/seismonitor_dld/"
                  "{network}/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
md = MseedDownloader([client])
md.download(rest,mseed_storage,
            ppc_restrictions=ppc_rest ,
            n_processor=16,
            concurrent_feature="thread")