import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.magnitude.mag import Magnitude,MwPhysicalMagParams, MwProcessingMagParams
from SeisMonitor.core.objects import WaveformRestrictions,Provider

sgc_client = FDSNClient('http://10.100.100.13:8091')
# sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    # starttime=UTCDateTime("2021-06-01T00:00:00.000000Z"),
                    # endtime=UTCDateTime("2022-01-01T00:00:00.000000Z"),
                    starttime=UTCDateTime("2019-12-01T00:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-02T00:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    # filter_domain= [-83.101,-64.549,-2.229,14.945],
                    # filter_domain= [-76.536,-71.168,7.758,11.823],
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml )

catalog = "/home/emmanuel/EDCT/test/associations/associations.xml"
out = "/home/emmanuel/EDCT/test"

mag = Magnitude([sgc_provider],catalog,out) #catalog,providers,out
cat = mag.get_Ml(mag_type="RSNC",trimmedtime=5,out_format="SC3ML")


# physparams = MwPhysicalMagParams(vp=8200,
#                             p_radiation_pattern=0.05,
#                             )
# procparams = MwProcessingMagParams(
#                         time_before_pick = 2,
#                         time_after_pick = 15,
#                         only_proc_p_pick=True)
# cat = mag.get_Mw(physparams,procparams,out_format="SC3ML")
print(cat)