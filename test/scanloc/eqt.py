import sys

Seismonitor_path = "/home/emmanuel/Ecopetrol/SeisMonitor"
sys.path.insert(0,Seismonitor_path)

import os
from obspy.clients.filesystem.sds import Client as SDS_Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.downloader.utils import DownloadRestrictions
from SeisMonitor.scanloc.utils import EQTobj
from SeisMonitor.scanloc.monitor import Monitor
# from SeisMonitor.tools import stats
from SeisMonitor.locator.seiscomp import merge_csv

sds_archive = "/home/emmanuel/archive/sds"
eqt_model = '/home/emmanuel/EQTransformer/ModelsAndSampleData/EqT_model.h5'
data = "/home/emmanuel/Ecopetrol/SeisMonitor/data"
out_folder = "/home/emmanuel/Ecopetrol/SeisMonitor/out_test/scanloc/1"


resp = os.path.join(data,"metadata/RESP.EY")
# xml_file = os.path.join(data,"metadata/RESP.EY.xml")
vel_file = os.path.join(data,"vel_model/castilla_vel1d.csv")
station0_file = os.path.join(data,"castilla/STATION0.HYP")
info_dir = os.path.join(out_folder,'downloads')
picks_dir = os.path.join(out_folder,'detections')
events_dir = os.path.join(out_folder,'events')

client = SDS_Client(sds_archive,
                    sds_type='D', format='MSEED',)
restrictions = DownloadRestrictions(network="EY",
                        station="CA*",
                        location="00",
                        channel="HH*",
                        starttime=UTCDateTime("2020-11-04T03:20:00"),
                        endtime=UTCDateTime("2020-11-04T03:22:00"),
                        # starttime=UTCDateTime("2020-11-01T03:00:00."),
                        # endtime=UTCDateTime("2020-11-01T04:00:00.0"),
                        chunklength_in_sec=3600,
                        overlap_in_sec=None,
                        groupby='{network}.{station}.{channel}')

eqtobj = EQTobj(model_path = eqt_model,
            chunk_size = 3600,
            n_processor = 4,
            overlap = 0.5,
            detection_threshold =0.01,
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

