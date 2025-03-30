import sys
repository_path = r"/home/edc240000/SeisMonitor"  
sys.path.insert(0,repository_path)

import os
from SeisMonitor.monitor.picker.ai import PhaseNet,PhaseNetObj
from SeisMonitor.monitor.picker import utils as piut
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

archive = "/home/edc240000/SeisMonitor/test/out/download/fdsn"
mseed_storage = os.path.join(archive,"downloads","{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

dataset = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"data")
pnet_path = "/home/edc240000/PhaseNet"
pnet_model_path = "/home/edc240000/PhaseNet/model/190703-214543"
pnetobj = PhaseNetObj(
            pnet_path=pnet_path,
            model_path=pnet_model_path,
            P_threshold=0.7, S_threshold=0.6,
            batch_size=100
            ) 

mseed_storage = os.path.join(archive,"downloads")
json_dir = os.path.join(archive,"stations")
out_dir = os.path.join(archive,"picks","pnet")
result = os.path.join(archive,"picks","seismonitor_picks.csv")

pnet = PhaseNet(pnetobj)
pnet.pick(mseed_storage,json_dir,out_dir)