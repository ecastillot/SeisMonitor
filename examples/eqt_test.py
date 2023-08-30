import os
import zipfile
from SeisMonitor.utils4examples import clone_seismonitor_data
import os
from SeisMonitor.monitor.picker.ai import PhaseNet,PhaseNetObj
from SeisMonitor.monitor.picker import utils as piut

import os
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut

monitor_path = "/home/emmanuel/AIOBS/SeisMonitor-dataset/downloads"

downloads = os.path.join(monitor_path ,"downloads")
stations = os.path.join(monitor_path ,"stations")

eqt_model = "/home/emmanuel/AIOBS/SeisMonitor-models/EQTransformer_models/EqT_model.h5"

eqtobj = EQTransformerObj(model_path=eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 20,
            number_of_plots = 0,
            plot_mode = 1 ) 

out_dir = os.path.join(monitor_path ,"picks","eqt")
result = os.path.join(monitor_path ,"picks","eqt","seismonitor_picks.csv")

eqt = EQTransformer(eqtobj)
eqt.pick(downloads,stations,out_dir)
piut.eqt_picks_2_seismonitor_fmt(out_dir,downloads,result)