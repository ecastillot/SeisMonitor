import os
import math
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.magnitude.mag import Magnitude,MwPhysicalMagParams, MwProcessingMagParams
from SeisMonitor.core.objects import WaveformRestrictions,Provider

monitor_path = "/home/emmanuel/XO_monitor_results"


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
ml_params = {"a":1.019,"b":0.0016,"r_ref":140} #YOU WILL NEED TO FIND THE VALUES ACCORDING TO THE PAPERS
                                              ### SEE HERE THE EQUATION ANT THE MEANNING OF THE VALUES
                                              #https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/5.magnitude.ipynb#scrollTo=4xsALpmpdmAk

####### NO MODIFY THE FOLLOWING

sgc_client = FDSNClient('IRIS')
sgc_provider = Provider(sgc_client,sgc_rest)
nlloc_catalog_path = os.path.join(monitor_path,"loc","nlloc","nlloc_loc.xml")
out_dir = os.path.join(monitor_path,"magnitude","nlloc","Ml")
mag = Magnitude([sgc_provider],nlloc_catalog_path ,
out_dir) #catalog,providers,out


k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
                    ml_params["b"]* (ml_params["r_ref"]-100) +3

Ml = lambda ampl,epi_dist : math.log10(ampl * 1e3) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
                                ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k

cat = mag.get_Ml(mag_type=Ml ,
            trimmedtime=5, #seconds after pick S to trim the signal
            out_format="SC3ML")
print(cat)