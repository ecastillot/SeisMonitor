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
                order=("CM","URMC"),
                fontsize=6,
                show=True):
        picker_csv = self.picks[picker]

        streams = ut.get_ordered_streams(self.providers,order,
                                    starttime,endtime)
        
        n_traces = len(list(streams.keys()))

        # Pcmap = cm.get_cmap("Blues",10)
        # Scmap = cm.get_cmap("Reds",10)

        fig, ax = plt.subplots(n_traces, sharex=True,
                    gridspec_kw = {'wspace':0, 'hspace':0})
        
        for i,(strid,st) in enumerate(streams.items()):
            tr = ut.get_proc_tr(st)
            STARTTIME = mdates.date2num(starttime.datetime)

            x = tr.times("matplotlib")-STARTTIME
            ax[i].plot(x, tr.data, "k-")
            ax[i].xaxis_date()
            ax[i].set_xlim([min(x),max(x)])
            ymin, ymax = ax[i].get_ylim()
            ax[i].set_yticks([])
            ax[i].set_ylabel(strid,rotation=0,labelpad=24,fontsize=fontsize)

            df = ut.get_picks(picker_csv,starttime,endtime,
                        select_stations=[tr.stats.station])
            for j, (_, row) in enumerate(df.iterrows()):
                pick = mdates.date2num(row['arrival_time'])-STARTTIME
                if row['phasehint'].upper() == 'P':
                    # color = Pcmap(row['probability']-0.001) # -0.001 removes a bug in Cbar 
                    color = "blue" # -0.001 removes a bug in Cbar 
                elif row['phasehint'].upper() == 'S':
                    # color = Scmap(row['probability']-0.001)
                    color = "red"
                
                ax[i].vlines(pick, ymin, ymax, color=color, 
                            # linewidth=0.7, 
                            linewidth=5, 
                            label=row['phasehint'])
                ax[i].set_facecolor('lightgray')

            # text = ".".join((network,station,
                 
            #         location,channel))
            # text = text.rjust(len(text) + 18)
        text = "starttime: " + starttime.strftime("%Y-%m-%d %H:%M:%S.%f")+\
                "\nendtime : " + endtime.strftime("%Y-%m-%d %H:%M:%S.%f")

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # ax[0].text(.03, .95, text,
        ax[0].text(0.8, 5, text,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax[0].transAxes,
                    backgroundcolor= "white",
                    fontdict={"fontsize":9},
                    bbox=props)
        ax[i].set_xlabel(f"dt", size=16)
        # ax.set_xlabel(f"dt", size=16)
        fig.autofmt_xdate()
        plt.tight_layout()
        plt.show()
        # print(json_info)
        # print(provider[0].waveform_restrictions.__dict__)

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