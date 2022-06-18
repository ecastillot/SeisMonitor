import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from SeisMonitor.plot.picks import plot_multiple_picker

client = Client(base_url='http://sismo.sgc.gov.co:8080/')

starttime = UTCDateTime("2020-09-22T00:00:00.000000Z")
endtime = UTCDateTime("2020-09-22T02:00:00.000000Z")
st = client.get_waveforms(network="CM",station="BAR2",
                                location="*",
                                channel="HHZ",
                                starttime=starttime,
                                endtime=endtime)
eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
phasenet_csv = "/home/emmanuel/EDCT/test/picks/pnet/results/seismonitor_picks.csv"
# csvs = [eqt_csv,sgc_csv,phasenet_csv]
csvs = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}
plot_multiple_picker(st,csvs)