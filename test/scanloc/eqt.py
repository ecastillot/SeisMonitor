import sys
import os
seismopath = "/home/emmanuel/test/seismo"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)


from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.utils import DownloadRestrictions
from SeisMonitor.scanloc.utils import EQTobj
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

eqt_model = os.path.join(seismopath,'EQTransformer/ModelsAndSampleData/EqT_model.h5')
eqtobj = EQTobj(model_path = eqt_model,
            chunk_size = 3600,
            n_processor = 4,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 1,
            number_of_plots = 1,
            plot_mode = 1 )  
            
monitor = Monitor(providers = [client],
                    restrictions=restrictions,
                    info_dir=info_dir,
                    picks_dir=picks_dir,
                    events_dir=events_dir)

monitor.make_json(from_xml=resp)
monitor.download()
monitor.eqt_picker(eqtobj,rm_downloads=False)
monitor.eqt_associator()



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

