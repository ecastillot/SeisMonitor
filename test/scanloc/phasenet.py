import sys
import os
seismopath = "/home/emmanuel/test/seismo"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)


from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.utils import DownloadRestrictions
from SeisMonitor.scanloc.utils import PhaseNetobj
from SeisMonitor.scanloc.monitor import Monitor

# route = "/home/emmanuel/archive/sds"
route = 'http://sismo.sgc.gov.co:8080'
resp = os.path.join(f"{seismonitor}/data","metadata/CM.xml")
out_folder = f"{seismonitor}/out/scanloc/test_CM"

# outs
info_dir = os.path.join(out_folder,'downloads')
picks_dir = os.path.join(out_folder,'detections')
events_dir = os.path.join(out_folder,'events')

## locator parameters (SEISAN software is mandatory to locate events)
vel_file = os.path.join("../data","vel_model/CM_vel1d.csv")
station0_file = os.path.join("../data","CM/STATION0.HYP")


client = Client(route)
restrictions = DownloadRestrictions(network="CM",
                          station="BAR2",
                          location="*",
                          channel="*",
                          starttime=UTCDateTime("2020-12-09T08:15:10.0"),
                          endtime=UTCDateTime("2020-12-09T08:16:10.0"),
                          chunklength_in_sec=43200,
                          overlap_in_sec=None,
                          groupby='{network}.{station}.{location}',
                          threshold=60,
                          location_preferences=["","00","20","10"],
                          channel_preferences=["HH","BH"],
                          to_pick=(1,0.3))      

pnet_model = os.path.join(seismopath,'PhaseNet/model/190703-214543')
pnetobj = PhaseNetobj(model_path = pnet_model ,
                mode='pred',
                P_threshold=0.1,
                S_threshold=0.1,
                batch_size=1, 
                plot=True, 
                save_result=False) 
monitor = Monitor(providers = [client],
                    restrictions=restrictions,
                    info_dir=info_dir,
                    picks_dir=picks_dir,
                    events_dir=events_dir)

monitor.make_json(from_xml=resp)
monitor.download()
monitor.mv_downloads2onefolder()
monitor.make_datalist()
monitor.phasenet_picker(pnetobj)