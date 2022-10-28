import os
import math
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.magnitude.mag import Magnitude,MwPhysicalMagParams, MwProcessingMagParams
from SeisMonitor.core.objects import WaveformRestrictions,Provider

archive = "./out/download/fdsn"

sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_rest = WaveformRestrictions(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-25T01:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    filter_domain= [-83.101,-64.549,-2.229,14.945],
                    )
sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
sgc_provider = Provider(sgc_client,sgc_rest)

catalog = os.path.join(archive,"loc","nlloc","LOC.xml")
out_dir = os.path.join(archive,"mag","Ml")

mag = Magnitude([sgc_provider],catalog,out_dir) #catalog,providers,out


ml_params = {"a":1.019,"b":0.0016,"r_ref":140} #ojeda
k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
                    ml_params["b"]* (ml_params["r_ref"]-100) +3
Ml = lambda ampl,epi_dist : math.log10(ampl * 1e3) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
                                ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k

cat = mag.get_Ml(mag_type=Ml ,
            trimmedtime=5,
            out_format="SC3ML")