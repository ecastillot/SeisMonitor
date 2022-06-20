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
from . import utils as ut
# import utils as ut

class Tracer():
    def __init__(self,client,picks,
                xml=None):
        self.client = client
        self.picks = picks
        self.xml = xml
    
    def plot(self,network,station,
                            location,channel,
                            starttime,
                            endtime,
                            st_proc = {
                            "normalize":True,
                            "merge":{"fill_value":0},
                            "detrend":{"type":"demean"},
                            "taper":{"max_percentage":0.001, 
                                    "type":"cosine", 
                                    "max_length":2} , 
                            "filter":{"type":'bandpass', 
                                    "freqmin" : 1, 
                                    "freqmax" : 45, 
                                    "corners":2, 
                                    "zerophase":True},
                            },
                            show=True
                            ):
        st = self.client.get_waveforms(network,station,
                                    location,
                                    channel,
                                    starttime,
                                    endtime)
        fig = ut.plot_multiple_picker(st,self.picks,st_proc)
        if show:
            plt.show()
        return fig

    def plot_all_stations(self,picker,
                            starttime,endtime,
                            select_networks=[],
                            select_stations=[],
                            filter_networks=[],
                            filter_stations=[],
                            st_proc = {
                            "normalize":True,
                            "merge":{"fill_value":0},
                            "detrend":{"type":"demean"},
                            "taper":{"max_percentage":0.001, 
                                    "type":"cosine", 
                                    "max_length":2} , 
                            "filter":{"type":'bandpass', 
                                    "freqmin" : 1, 
                                    "freqmax" : 45, 
                                    "corners":2, 
                                    "zerophase":True},
                            },
                            show=True):
        pass

if __name__ == "__main__":
    client = Client(base_url='http://sismo.sgc.gov.co:8080/')

    eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    phasenet_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    picks = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}

    Trace(client,picks)
    # # client = Client(base_url='http://10.100.100.232:8091')
    # # starttime = UTCDateTime("20191224T185959")
    # starttime = UTCDateTime("20191224T190600")
    # endtime = UTCDateTime("20191224T191359")
    # st = client.get_waveforms(network="CM",station="URMC",
    #                                 location="*",
    #                                 channel="HHZ",
    #                                 starttime=starttime,
    #                                 endtime=endtime)

    # # csvs = [eqt_csv,sgc_csv,phasenet_csv]

    # # df = get_picks(phasenet_csv,starttime,endtime)
    # # print(df)
    # fig = plot_multiple_picker(st,csvs)
    # fig.savefig("/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/plot/mpt.png",dpi=300)
    # # plt.show()