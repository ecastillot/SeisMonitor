import math
from pathlib import Path
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.monitor.magnitude.mag import Magnitude
from SeisMonitor.core.objects import WaveformRestrictions,Provider

monitor_path = Path(__file__).parent / "sm"

loc_path = monitor_path / "locations"
mag_path = monitor_path / "magnitude"

nlloc_output_path = loc_path / "nlloc"
xml_nlloc_path = nlloc_output_path / "nlloc_catalog.xml"

client = FDSNClient("http://rtserve.beg.utexas.edu")
rest = WaveformRestrictions(network="*",
                    station="PB13,PB17,PB20,PB23,PB34,PB58,PB40,PB24,PB25,PB26,WB03",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2022-11-16T21:00:00.000000Z"),
                    endtime=UTCDateTime("2022-11-16T23:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    )


provider = Provider(client,rest)
mag = Magnitude([provider],xml_nlloc_path,mag_path) #catalog,providers,out

# Use your own Ml formula, here is an example for local magnitude, but make sure to check if the built-in formula is suitable for your region and data.
# For Texas using Kavoura et al (2020)
Ml = lambda ampl,epi_dist : math.log10(ampl*1e3 ) + 1.54 * math.log10(epi_dist) + 0.0002*(epi_dist-100) - 0.08

cat = mag.get_Ml(mag_type=Ml ,
            trimmedtime=5, #seconds after pick S to trim the signal
            out_format="SC3ML")
print(cat)