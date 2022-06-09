import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)


from SeisMonitor.monitor.picker.picker import PhaseNet,PhaseNetObj
from obspy.core.utcdatetime import UTCDateTime

pnet_model = os.path.join("/home/emmanuel/test/seismo",'PhaseNet/model/190703-214543')
pnetobj = PhaseNetObj(model_path = pnet_model ,
                mode='pred',
                P_threshold=0.75,
                S_threshold=0.75,
                batch_size=100, 
                plot=False, 
                save_result=False) 

mseed_storage = "/home/emmanuel/EDCT/test/downloads"
json_path = "/home/emmanuel/EDCT/test/json/test.json"
out_dir = "/home/emmanuel/EDCT/test/picks/pnet"
pnet = PhaseNet(mseed_storage,json_path,out_dir)
pnet.mv_downloads2onefolder()
pnet.make_datalist()
pnet.picker(pnetobj)