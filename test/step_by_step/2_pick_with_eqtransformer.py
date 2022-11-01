import os
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut
from obspy.clients.fdsn import Client as FDSNClient
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.core.client import LocalClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

archive = "./out/download/fdsn"
json_path = os.path.join(archive,"json/stations.json")
mseed_storage = os.path.join(archive,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")

dataset = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"data")
eqt_model = os.path.join(dataset,"models",'EqT_model.h5')
eqtobj = EQTransformerObj(model_path = eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 20,
            number_of_plots = 0,
            plot_mode = 1 ) 

mseed_storage = os.path.join(archive,"downloads")
json_dir = os.path.join(archive,"stations")
out_dir = os.path.join(archive,"picks","eqt")
result = os.path.join(archive,"picks","seismonitor_picks.csv")

# eqt = EQTransformer(mseed_storage,json_path,out_dir)
eqt = EQTransformer(eqtobj)
eqt.pick(mseed_storage,json_dir,out_dir)
piut.eqt_picks_2_seismonitor_fmt(out_dir,mseed_storage,result)