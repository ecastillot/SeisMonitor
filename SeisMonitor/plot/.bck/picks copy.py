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

def get_multiple_picker_figure():
    fig = plt.figure(figsize=(12,12))
    gs1 = gs.GridSpec(1, 1, figure=fig)
    gs1.update(left=0.05, right=0.9, wspace=0.05)

    gss = gs.GridSpecFromSubplotSpec(4, 1,subplot_spec=gs1[0], hspace=0.0,height_ratios=[4,0.6,0.6,0.6])
    ax0 = fig.add_subplot(gss[0])
    ax0.tick_params(axis="x", labelbottom=0)
    ax0.set_yticks([])
    ax1 = fig.add_subplot(gss[1], sharex=ax0)
    ax1.tick_params(axis="x", labelbottom=0)
    ax1.set_yticks([])
    ax2 = fig.add_subplot(gss[2], sharex=ax0)
    ax2.tick_params(axis="x", labelbottom=0)
    ax2.set_yticks([])
    ax3 = fig.add_subplot(gss[3], sharex=ax0)
    ax3.set_yticks([])
    fig.autofmt_xdate()

    gs2 = gs.GridSpec(1, 1)
    gs2.update(left=0.92, right=0.96, hspace=0.2)
    gss2 = gs.GridSpecFromSubplotSpec(1, 2,subplot_spec=gs2[0], hspace=0.5,width_ratios=[1,1])
    ax4 = fig.add_subplot(gss2[0])
    ax4.set_yticks([])
    ax4.set_xticks([])
    ax5 = fig.add_subplot(gss2[1],sharey=ax4)
    # ax4 = plt.subplot(gs2[:, 0])
    # ax5 = plt.subplot(gs2[:, 1])
    ax5.set_yticks([])
    ax5.set_xticks([])
    ax5.text(1.1,0.5,"Probability", size=10,
                           verticalalignment='center', rotation=270)

    return fig, [ax0,ax1,ax2,ax3,ax4,ax5]

def plot_multiple_picker(st,info):

    fig,ax = get_multiple_picker_figure()

    st.merge(fill_value=0) 
    st.detrend('demean')
    st.taper(max_percentage=0.001, type='cosine', max_length=2) 
    st.filter(type='bandpass', freqmin = 1.0, freqmax = 45, corners=2, zerophase=True)
    st.trim(min([tr.stats.starttime for tr in st]), max([tr.stats.endtime for tr in st]), pad=True, fill_value=0)

    tr = st[0]

    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    xlabel_date = starttime.strftime("%Y/%m/%d %H:%M:%S")
    STARTTIME = mdates.date2num(starttime.datetime)

    # tr.data = scale(tr,factor=1)
    x = tr.times("matplotlib")-STARTTIME
    ax[0].plot(x, tr.data, "k-")
    # ax[0].plot(tr.times("matplotlib"), tr.data, "k-")
    ax[0].xaxis_date()
    # ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%M'))
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%S'))
    ax[0].set_xlim([min(x),max(x)])
    ax[3].set_xlabel(f"Seconds from {xlabel_date}", size=20)
    # ax[3].set_xlabel(f"Minutes from {xlabel_date}", size=20)
    ymin, ymax = ax[0].get_ylim()


    # print(df)
    # https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
    Pcmap = cm.winter.reversed()
    Scmap = cm.autumn.reversed()
    for i,(picker,csv) in enumerate(info.items(),1):
        df = ut.get_picks(csv,starttime,endtime,
                        station_list=[tr.stats.station])
        for j, (_, row) in enumerate(df.iterrows()):
            pick = mdates.date2num(row['arrival_time'])-STARTTIME
            # pick = row['arrival_time']  
            if row['phasehint'].upper() == 'P':
                color = Pcmap(row['probability']-0.001) # -0.001 por un bug en Cbar 
            elif row['phasehint'].upper() == 'S':
                color = Scmap(row['probability']-0.001)

            print(row['probability'])
            ax[i].vlines(pick, ymin, ymax, color=color, 
                            linewidth=0.5, 
                            label=row['phasehint'])
            ax[i].set_facecolor('lightgray')
            ax[i].set_yticks([])
            ax[i].set_ylabel(picker,rotation=0,labelpad=10,fontsize=8)

    # # fig.autofmt_xdate()


    Pclb = fig.colorbar(cm.ScalarMappable(cmap=Pcmap), cax=ax[4],orientation="vertical")
    Pclb.ax.tick_params(labelsize=0.7)
    Pclb.ax.set_title('P',fontsize=20)
    Sclb = fig.colorbar(cm.ScalarMappable(cmap=Scmap), cax=ax[5],orientation="vertical")
    Sclb.ax.set_title('S',fontsize=20)
    # Sclb.ax.set_title('Probability',fontsize=10)

    fig.autofmt_xdate()
    plt.show()
