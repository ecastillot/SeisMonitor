import sys
Seismonitor_path = "/home/emmanuel/Ecopetrol/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from obspy.clients.fdsn import Client as FDSN_Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.obspy import MseedDownloader
from obspy.clients.fdsn.mass_downloader import Restrictions
from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain

# client = FDSN_Client(base_url="IRIS")
# # client = FDSN_Client(base_url="IRIS")
# restrictions = Restrictions(starttime=UTCDateTime("2016-12-24T19:00:00.000000Z"),
#                             endtime=UTCDateTime("2016-12-24T20:00:00.000000Z"),
#                             network="YU", 
#                             station="*",
#                             # location="*", channel="*",
#                             chunklength_in_sec=3600,
#                             reject_channels_with_gaps=False,
#                             minimum_length=0,
#                             minimum_interstation_distance_in_m=0.0,
#                             channel_priorities=["HH[ZNE]", "BH[ZNE]","EH[ZNE]","HN[ZNE]"],
#                             location_priorities=["", "00", "20", "10"])

client = FDSN_Client(base_url="http://sismo.sgc.gov.co:8080")
# client = FDSN_Client(base_url="IRIS")
restrictions = Restrictions(starttime=UTCDateTime("2017-12-24T19:00:00.000000Z"),
                            endtime=UTCDateTime("2017-12-24T20:00:00.000000Z"),
                            network="CM", 
                            station="*",
                            # location="*", channel="*",
                            chunklength_in_sec=3600,
                            reject_channels_with_gaps=False,
                            minimum_length=0,
                            minimum_interstation_distance_in_m=0.0,
                            channel_priorities=["HH[ZNE]", "BH[ZNE]","EH[ZNE]","HN[ZNE]"],
                            location_priorities=["", "00", "20", "10"])

# domain = RectangularDomain(minlatitude=35.50, maxlatitude=35.60,
#                            minlongitude=-117.80, maxlongitude=-117.40)
domain = RectangularDomain(minlatitude=7.758, maxlatitude=11.823,
                           minlongitude=-76.536, maxlongitude=-71.168)

mseed_storage = "/home/emmanuel/EDCT/test/downloads/"
stationxml_storage = "/home/emmanuel/EDCT/test/stations"

mseed_dl = MseedDownloader([client])
mseed_dl.download(domain=domain,
                restrictions=restrictions, 
                mseed_storage=mseed_storage,
                stationxml_storage=stationxml_storage,
                workers=1,parallel_mode="thread")