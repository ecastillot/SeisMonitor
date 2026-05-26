import os
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut

monitor_path = "/home/emmanuel/XO_monitor_results"
eqt_model = "/home/emmanuel/XO_monitor_results/picking_models/eqt/EqT_original_model.h5"
eqtobj = EQTransformerObj(model_path=eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 100,
            number_of_plots = 10,
            plot_mode = None ) 

####### NO MODIFY THE FOLLOWING

out_dir = os.path.join(monitor_path ,"picks","eqt")
result = os.path.join(monitor_path ,"picks","eqt","seismonitor_picks.csv")
downloads = os.path.join(monitor_path ,"downloads")
stations = os.path.join(monitor_path ,"stations")

eqt = EQTransformer(eqtobj)
eqt.pick(downloads,stations,out_dir)
piut.eqt_picks_2_seismonitor_fmt(out_dir,downloads,result)