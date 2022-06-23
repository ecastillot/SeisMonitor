from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy.realtime.signal import scale
# from obspy.core.stream import Stream

import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import datetime as dt
import pandas as pd
import numpy as np
import json
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
import matplotlib.dates as mdates
from . import utils as ut
# import utils as ut
from matplotlib.dates import DateFormatter
from matplotlib.transforms import blended_transform_factory

class Tracer():
    def __init__(self,provider,picks):
        self.provider = provider
        self.picks = picks
    
    def plot(self,show=True):
        waveform_restrictions = self.provider.waveform_restrictions
        processing = self.provider.processing
        st = self.provider.client.get_waveforms(waveform_restrictions.network,
                                waveform_restrictions.station,
                                waveform_restrictions.location,
                                waveform_restrictions.channel,
                                waveform_restrictions.starttime,
                                waveform_restrictions.endtime)
        fig = ut.plot_multiple_picker(st,self.picks,processing )
        if show:
            plt.show()
        return fig

class Streamer():
    def __init__(self,providers,picks):
        self.providers = providers
        self.picks = picks

    def plot(self,picker,starttime,endtime,
                order,
                fontsize=6,
                show=True):

        picker_csv = self.picks[picker]
        streams = ut.get_ordered_streams(self.providers,order,
                                    starttime,endtime)
        fig = ut.get_streamer_plot(streams,picker_csv,
                                    starttime,endtime,
                                    fontsize,show)
        return fig

    def plot_by_station(self,picker,netsta,n_picks=50,
                        align="S",
                        phase_second = 2,
                        window = 15,
                        show=True):

        picker_csv = self.picks[picker]
        df = ut.get_picks(picker_csv,
                    select_networks=[netsta[0]],
                    select_stations=[netsta[1]])
        fig = ut.get_plot_by_station(self.providers,netsta,
                            df,n_picks,align,phase_second,
                            window,
                            show)
        return fig
        
if __name__ == "__main__":
    client = Client(base_url='http://sismo.sgc.gov.co:8080/')

    eqt_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    sgc_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    phasenet_csv = "/home/emmanuel/EDCT/test/picks/eqt/seismonitor_picks.csv"
    picks = {"EQT":eqt_csv,"PNET":phasenet_csv,"SGC":sgc_csv}

    # Trace(client,picks)
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