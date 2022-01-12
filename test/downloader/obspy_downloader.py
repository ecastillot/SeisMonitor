import sys
Seismonitor_path = "/home/emmanuel/Ecopetrol/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from obspy.clients.fdsn import Client as FDSN_Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.obspy import MseedDownloader
from obspy.clients.fdsn.mass_downloader import Restrictions
from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain

client = FDSN_Client(base_url="IRIS")
restrictions = Restrictions(starttime=UTCDateTime(2020, 9, 1),
                            endtime=UTCDateTime(2020, 9, 2),
                            network="GS", 
                            station="CA10",
                            # location="*", channel="*",
                            chunklength_in_sec=86400,
                            reject_channels_with_gaps=False,
                            minimum_length=0.0,
                            minimum_interstation_distance_in_m=0.0,
                            channel_priorities=["HH[ZNE]", "BH[ZNE]","EH[ZNE]","HN[ZNE]"],
                            location_priorities=["", "00", "20", "10"])

domain = RectangularDomain(minlatitude=35.50, maxlatitude=35.60,
                           minlongitude=-117.80, maxlongitude=-117.40)

mseed_storage = "/home/emmanuel/Ecopetrol/SeisMonitor/out/downloads/obspy_dld/waveforms"
stationxml_storage = "/home/emmanuel/Ecopetrol/SeisMonitor/out/downloads/obspy_dld/stations"

mseed_dl = MseedDownloader([client])
mseed_dl.download(domain=domain,
                restrictions=restrictions, 
                mseed_storage=mseed_storage,
                stationxml_storage=stationxml_storage,
                workers=1,parallel_mode="thread")