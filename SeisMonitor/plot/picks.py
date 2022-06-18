from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy.realtime.signal import scale
# from obspy.core.stream import Stream

import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import datetime as dt
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
import matplotlib.dates as mdates
# from . import utils as ut
import utils as ut

# class PickFigure():
#     def __init__(self,):
#         pass


if __name__ == "__main__":
    # client = Client(base_url='http://10.100.100.232:8091')
    client = Client(base_url='http://sismo.sgc.gov.co:8080/')
    # starttime = UTCDateTime("20191224T185959")
    starttime = UTCDateTime("20191224T190600")
    endtime = UTCDateTime("20191224T191359")
    st = client.get_waveforms(network="CM",station="URMC",
                                    location="*",
                                    channel="HHZ",
                                    starttime=starttime,
                                    endtime=endtime)

    eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    phasenet_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    # csvs = [eqt_csv,sgc_csv,phasenet_csv]
    csvs = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}

    # df = get_picks(phasenet_csv,starttime,endtime)
    # print(df)
    fig = plot_multiple_picker(st,csvs)
    fig.savefig("/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/plot/mpt.png",dpi=300)
    # plt.show()