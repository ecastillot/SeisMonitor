import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)


from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.monitor.picker.picker import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut

  

eqt_model = os.path.join("/home/emmanuel/test/seismo",'EQTransformer/ModelsAndSampleData/EqT_model.h5')
eqtobj = EQTransformerObj(model_path = eqt_model,
            n_processor = 4,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 100,
            number_of_plots = 1,
            plot_mode = 1 )  
            
mseed_storage = "/home/emmanuel/EDCT/test/downloads/CM"
json_path = "/home/emmanuel/EDCT/test/json/test.json"
out_dir = "/home/emmanuel/EDCT/test/picks/eqt"
# eqt = EQTransformer(mseed_storage,json_path,out_dir)
# eqt.picker(eqtobj)



eqt_folder = "/home/emmanuel/EDCT/test/picks/eqt"
result = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
piut.eqt_picks_2_seismonitor_fmt(eqt_folder,result)