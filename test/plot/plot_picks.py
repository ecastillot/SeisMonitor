import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.plot.picks import Tracer,Streamer
from SeisMonitor.core.objects import WaveformRestrictions,Provider,Processing

# client = Client(base_url='http://sismo.sgc.gov.co:8080/')

# starttime = UTCDateTime("2020-09-22T00:00:00.000000Z")
# endtime = UTCDateTime("2020-09-22T02:00:00.000000Z")
# st = client.get_waveforms(network="CM",station="BAR2",
#                                 location="*",
#                                 channel="HHZ",
#                                 starttime=starttime,
#                                 endtime=endtime)
# eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
# sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
# phasenet_csv = "/home/emmanuel/EDCT/test/picks/pnet/results/seismonitor_picks.csv"
# # csvs = [eqt_csv,sgc_csv,phasenet_csv]
# csvs = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}
# plot_multiple_picker(st,csvs)

eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
phasenet_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
picks = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}

sgc_client = Client(base_url='http://sismo.sgc.gov.co:8080/')
sgc_rest = WaveformRestrictions(network="CM",
                                station="OCA",
                                location="*",
                                channel="HHZ",
                                starttime=UTCDateTime("2019-12-24T19:04:00.000000Z"),
                                endtime=UTCDateTime("2019-12-24T19:06:00.000000Z")
                    )
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
processing = Processing()
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml,processing=processing)


# tr = Tracer(sgc_provider,picks)
# fig = tr.plot()

# exit()

############################

sgc_rest = WaveformRestrictions(network="CM",
                                station="*",
                                location="*",
                                channel="HHZ",
                                starttime=UTCDateTime("2019-12-24T19:07:00.000000Z"),
                                endtime=UTCDateTime("2019-12-24T19:15:00.000000Z"),
                                location_preferences=["","00","20","10","40"],
                                channel_preferences=["HH","BH","EH","HN","HL"],
                                filter_networks=[], 
                                filter_stations=[],
                                filter_domain=[]
                                )
                                # mesetas -74.3719,3.2326
sgc_xml = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
processing = Processing()
sgc_provider = Provider(sgc_client,sgc_rest,xml=sgc_xml,processing=processing)
st = Streamer(providers=[sgc_provider],picks=picks)
st.plot("EQT",
        starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
        endtime=UTCDateTime("2019-12-24T20:00:00.000000Z"))