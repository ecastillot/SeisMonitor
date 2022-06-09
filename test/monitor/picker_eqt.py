import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)


from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.monitor.picker.picker import EQTransformer,EQTransformerObj

  

eqt_model = os.path.join("/home/emmanuel/test/seismo",'EQTransformer/ModelsAndSampleData/EqT_model.h5')
eqtobj = EQTransformerObj(model_path = eqt_model,
            n_processor = 4,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 1,
            number_of_plots = 1,
            plot_mode = 1 )  
            
mseed_storage = "/home/emmanuel/EDCT/test/downloads/CM"
json_path = "/home/emmanuel/EDCT/test/json/test.json"
out_dir = "/home/emmanuel/EDCT/test/picks/eqt"
eqt = EQTransformer(mseed_storage,json_path,out_dir)
eqt.picker(eqtobj)




# monitor.make_station0(xml=resp,vel_model=vel_file)
# catalog = monitor.hypocenter_locator(station0_file)


# catalog = "/home/emmanuel/results/SeisMonitor/11/magnitude/events/events.xml"

# monitor.compute_magnitude(client,catalog,resp)

# Mw_out = "/home/emmanuel/results/SeisMonitor/11/magnitude/events/events_ml.xml"
# E_df = stats.quakexml2fmt(Mw_out)
# E_df.to_csv("/home/emmanuel/results/SeisMonitor/11/magnitude/events/events_ml.csv")
# print(E_df)



#### merge csv
# store = "/home/emmanuel/results/SeisMonitor/11/detections/eqt"
# out = "/home/emmanuel/analisis/automatic.csv"
# merge_csv(store,"eqt","pick","event_start_time",out)

