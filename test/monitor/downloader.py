import sys
Seismonitor_path = "/home/emmanuel/EDCT/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
# sgc_client = FDSNClient('http://10.100.100.232:8091')
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
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-24T20:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )

root = "/home/emmanuel/RSNC-2016"
fmt = '{year}/{network}/{station}/{channel}.{sds_type}/{network}.{station}.{location}.{channel}.{sds_type}.{year}.{doy:03d}'
local_client = LocalClient(root,fmt)

# print(rest.domain)
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml)
carma_provider = Provider(carma_client,carma_rest)

providers = [carma_provider,sgc_provider]

md = MseedDownloader(providers=[sgc_provider])
json_path = "/home/emmanuel/test_downloads/test_metadata"
mseed_storage = ("/home/emmanuel/test_downloads/"
                  "{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
md.make_inv_and_json(json_path)
md.download(mseed_storage,chunklength_in_sec=3600,n_processor=None)

# client2 = FDSNClient(base_url="IRIS")

# md = MseedDownloader([client],rest,ppc_restrictions=ppc_rest)
# md.download(mseed_storage,
#             n_processor=16,
#             concurrent_feature="thread")
# json_path = "/home/emmanuel/EDCT/test/json/test.json"
# md.make_json(json_path,from_xml=xml)