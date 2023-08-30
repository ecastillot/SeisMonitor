import os
import zipfile
from SeisMonitor.utils4examples import clone_seismonitor_data
import os
from SeisMonitor.monitor.picker.ai import PhaseNet,PhaseNetObj
from SeisMonitor.monitor.picker import utils as piut

monitor_path = "/home/emmanuel/AIOBS/SeisMonitor-dataset/downloads"

downloads = os.path.join(monitor_path ,"downloads")
stations = os.path.join(monitor_path ,"stations")
models = "/home/emmanuel/AIOBS/SeisMonitor-models"
pnet_model = os.path.join(models,"PhaseNet_models","190703-214543")
pnet_path = os.path.join(monitor_path ,"PhaseNet")


# piut.clone_aipicker("PhaseNet",pnet_path)

pnetobj = PhaseNetObj(pnet_path=pnet_path,
            model_path=pnet_model,
            P_threshold=0.7, S_threshold=0.6,
            batch_size=100
            ) 

out_dir = os.path.join(monitor_path,"picks","pnet")
result = os.path.join(monitor_path,"picks","pnet","seismonitor_picks.csv")

pnet = PhaseNet(pnetobj)
pnet.pick(downloads,stations,out_dir)