import pandas as pd
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.inventory.inventory import (Inventory,read_inventory)
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import json
from matplotlib.dates import date2num,AutoDateLocator, num2date, datestr2num, HourLocator, MinuteLocator
def order_by_coord(coord,csv_id,dataless):
    lat,lon = coord
    df_id = pd.read_csv(csv_id)
    stations = df_id["station_code"].to_list()
    inv = read_inventory(dataless)
    stas= []
    lats = []
    lons = []
    for net in inv:
        for sta in net:
            if sta.code in stations:
                stas.append(sta.code)
                lats.append(sta.latitude)
                lons.append(sta.longitude)

    data = {"station":stas,"latitude":lats,"longitude":lons}
    df = pd.DataFrame.from_dict(data)

    gps2da_gen = lambda y: gps2dist_azimuth(y.latitude,y.longitude,lat,lon)
    gps2da = df.apply(gps2da_gen,axis=1).tolist()
    df[["r","az","baz"]] = pd.DataFrame(gps2da)
    df["r"] = df["r"]/1e3

    df = df.drop_duplicates(subset=["station"],
                                ignore_index=True)
    df = df.sort_values("r",ignore_index=True,ascending=False)
    df = df[["station","r"]]
    df['id'] = df.index
    # print(df)
    return df

def plot_asso_pnet(coord,csv_id,dataless,df,ev_df,starttime,endtime):
    df_id = order_by_coord(coord,csv_id,dataless)

    df = df[(df["time"]<= endtime) & (df["time"]>= starttime)]

    f = open('/home/emmanuel/Tesis/auto/aipicker/filters/phasenet_filter_prob.json')
    data = json.load(f)


    f_df = []
    nf_df = []
    for station,values in data.items():
        val_p = values["P"]
        val_s = values["S"]
        # val_p = 0.7
        # val_s = 0.6
        fdf_p = df[(df["phasehint"]=="P") & \
                        (df['uncertainty']>=val_p) & \
                        (df['station'].replace(' ', '')==station) ]
        fdf_s = df[(df["phasehint"]=="S") & \
                        (df['uncertainty']>=val_s) & \
                        (df['station'].replace(' ', '')==station) ]
        
        nfdf_p = df[(df["phasehint"]=="P") & \
                        ~(df['uncertainty']>=val_p) & \
                        (df['station'].replace(' ', '')==station) ]
        nfdf_s = df[(df["phasehint"]=="S") & \
                        ~(df['uncertainty']>=val_s) & \
                        (df['station'].replace(' ', '')==station) ]

        if  not fdf_p.empty:
            f_df.append(fdf_p)
        if  not fdf_s.empty:
            f_df.append(fdf_s)
        if  not nfdf_p.empty:
            nf_df.append(nfdf_p)
        if  not nfdf_s.empty:
            nf_df.append(nfdf_s)

    f_df = pd.concat(f_df,ignore_index=True)
    nf_df = pd.concat(nf_df,ignore_index=True)

    ev_df_p = ev_df[(ev_df["time_pick_p"]<= endtime) & (ev_df["time_pick_p"]>= starttime)]
    ev_df_s = ev_df[(ev_df["time_pick_s"]<= endtime) & (ev_df["time_pick_s"]>= starttime)]

    f_df = pd.merge(f_df,df_id)
    nf_df = pd.merge(nf_df,df_id)

    print(f_df)
    # # print(nf_df)
    # exit()

    ev_df_p = pd.merge(ev_df_p,df_id)
    ev_df_s = pd.merge(ev_df_s,df_id)


    fig, ax = plt.subplots(figsize=(9,7))
    plt.rcParams["font.family"]="arial"
    ax.set_axisbelow(True)

    dates = nf_df["time"].to_list()
    idxs =  nf_df["id"].to_list()
    pick = ax.scatter(dates, idxs, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                label="Deleted phases",
                c="lightgray",
                marker='o', 
                edgecolor='lightgray',
                zorder=3,
                alpha=0.7)

    dates = f_df["time"].to_list()
    idxs =  f_df["id"].to_list()
    pick = ax.scatter(dates, idxs, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                c="gray",
                label="Unassociated phases",
                marker='o', 
                edgecolor='k',
                zorder=3,
                alpha=0.7)


    dates_p = ev_df_p["time_pick_p"].to_list()
    idxs_p =  ev_df_p["id"].to_list()

    dates_s = ev_df_s["time_pick_s"].to_list()
    idxs_s =  ev_df_s["id"].to_list()

    pick_p = ax.scatter(dates_p, idxs_p, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                marker='o', 
                c = "blue",
                label= "Associated P-phase",
                edgecolor='k',
                zorder=3,
                alpha=0.7)

    pick_s = ax.scatter(dates_s, idxs_s, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                marker='o', 
                c = "red",
                label= "Associated S-phase",
                edgecolor='k',
                zorder=3,
                alpha=0.7)

    ax.set_xlabel('Time', size=14)
    ax.set_ylabel('Index station', size=14)

    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5)) 
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.set_xlim(starttime,endtime)

    for label in ax.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='right')
    
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.set_yticks(np.arange(0, max(idxs)+5, 5.0))


    ax.grid(linestyle='--',zorder=20,axis="both",visible=True)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 - box.height * 0.05,
                    box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2)                

    plt.tight_layout() 

    return fig
    # plt.savefig("./phase",dpi=300)
    # plt.show()
    # print(df)
    # print(df_id)

def plot_asso_eqt(coord,csv_id,dataless,df,ev_df,starttime,endtime):
    df_id = order_by_coord(coord,csv_id,dataless)

    df_p = df[(df["p_arrival_time"]<= endtime) & (df["p_arrival_time"]>= starttime)]
    df_s = df[(df["s_arrival_time"]<= endtime) & (df["s_arrival_time"]>= starttime)]

    ev_df_p = ev_df[(ev_df["time_pick_p"]<= endtime) & (ev_df["time_pick_p"]>= starttime)]
    ev_df_s = ev_df[(ev_df["time_pick_s"]<= endtime) & (ev_df["time_pick_s"]>= starttime)]

    df_p = pd.merge(df_p,df_id)
    df_s = pd.merge(df_s,df_id)
    ev_df_p = pd.merge(ev_df_p,df_id)
    ev_df_s = pd.merge(ev_df_s,df_id)

    dates_p = df_p["p_arrival_time"].to_list()
    idxs_p =  df_p["id"].to_list()

    dates_s = df_s["s_arrival_time"].to_list()
    idxs_s =  df_s["id"].to_list()

    fig, ax = plt.subplots(figsize=(9,7))
    plt.rcParams["font.family"]="arial"
    ax.set_axisbelow(True)

    pick_p = ax.scatter(dates_p, idxs_p, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                c="gray",
                marker='o', 
                edgecolor='k',
                label="Unassociated phases",
                zorder=3,
                alpha=0.7)

    pick_s = ax.scatter(dates_s, idxs_s, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                c="gray",
                marker='o', 
                edgecolor='k',
                zorder=3,
                alpha=0.7)

    dates_p = ev_df_p["time_pick_p"].to_list()
    idxs_p =  ev_df_p["id"].to_list()

    dates_s = ev_df_s["time_pick_s"].to_list()
    idxs_s =  ev_df_s["id"].to_list()

    pick_p = ax.scatter(dates_p, idxs_p, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                marker='o', 
                c = "blue",
                edgecolor='k',
                label= "Associated P-phase",
                zorder=3,
                alpha=0.7)

    pick_s = ax.scatter(dates_s, idxs_s, 
                # s=m/2, c=c, cmap=cmap, norm=norm
                marker='o', 
                c = "red",
                label= "Associated S-phase",
                edgecolor='k',
                zorder=3,
                alpha=0.7)

    ax.set_xlabel('Time', size=14)
    ax.set_ylabel('Index station', size=14)

    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5)) 
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.set_xlim(starttime,endtime)
    
    for label in ax.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='right')

    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.set_yticks(np.arange(0, max(idxs_p)+5, 5.0))

    ax.grid(linestyle='--',zorder=20,axis="both",visible=True)


    box = ax.get_position()
    ax.set_position([box.x0, box.y0 - box.height * 0.05,
                    box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=3)                

    plt.tight_layout() 
    # plt.show()
    return fig
    # print(df)
    # print(df_id)

if __name__ == "__main__":
    # coord = (6.80,-73.15) #santos
    coord = (3.45,-74.19) #santos
    csv_id = "/home/emmanuel/Tesis/data/association/id.csv"
    dataless = "/home/emmanuel/Tesis/metadata/CM.xml"

    starttime = dt.datetime.strptime('20191224T190000', '%Y%m%dT%H%M%S')
    endtime = dt.datetime.strptime('20191224T210000', '%Y%m%dT%H%M%S')

    #######33order by coord
    # order_by_coord(coord,csv_id,dataless)

    ######## plot asso

    # ## pnet
    # pnet = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/2019/358/phasenet_picks.csv"
    # pnet = pd.read_csv(pnet)
    # pnet["time"] = pd.to_datetime(pnet["time"]).dt.tz_localize(None)
    # pnet = pnet[["station_code","phase_hint","time","uncertainty"]]
    # pnet = pnet.rename(columns={"station_code":"station",
    #                             "phase_hint":"phasehint"})

    # pnet_ev = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/csv/CM/2019/358/phasenet_events.csv"
    # pnet_ev = pd.read_csv(pnet_ev)
    # pnet_ev["time_pick_p"] = pd.to_datetime(pnet_ev["time_pick_p"])
    # pnet_ev["time_pick_s"] = pd.to_datetime(pnet_ev["time_pick_s"])
    # pnet_ev = pnet_ev[["station","pick_p","time_pick_p",
    #                     "pick_s","time_pick_s"]]

    # fig = plot_asso_pnet(coord,csv_id,dataless,pnet,pnet_ev,starttime,endtime)
    # fig.savefig("/home/emmanuel/tesis_analisis/tesis_analisis/figures/asso/phasenet_asso.png",dpi=300)

    # # # ### eqt
    # eqt = "/home/emmanuel/Tesis/auto/aipicker/picks/picks_1d/csv/CM/2019/358/eqt_picks.csv"
    # eqt = pd.read_csv(eqt)
    # eqt["p_arrival_time"] = pd.to_datetime(eqt["p_arrival_time"])
    # eqt["s_arrival_time"] = pd.to_datetime(eqt["s_arrival_time"])
    # eqt = eqt[["station","p_arrival_time",
    #                     "s_arrival_time"]]

    # eqt_ev = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/csv/CM/2019/358/eqt_events.csv"
    # eqt_ev = pd.read_csv(eqt_ev)
    # eqt_ev["time_pick_p"] = pd.to_datetime(eqt_ev["time_pick_p"])
    # eqt_ev["time_pick_s"] = pd.to_datetime(eqt_ev["time_pick_s"])
    # eqt_ev = eqt_ev[["station","pick_p","time_pick_p",
    #                     "pick_s","time_pick_s"]]

    # fig = plot_asso_eqt(coord,csv_id,dataless,eqt,eqt_ev,starttime,endtime)
    # fig.savefig("/home/emmanuel/tesis_analisis/tesis_analisis/figures/asso/eqt_asso.png",dpi=300)

    # # # ### sgc
    # sgc = "/home/emmanuel/Tesis/auto/sgc/picks/picks_1d/csv/CM/2019/358/sgc_picks.csv"
    # sgc = pd.read_csv(sgc)
    # sgc["p_arrival_time"] = pd.to_datetime(sgc["p_arrival_time"])
    # sgc["s_arrival_time"] = pd.to_datetime(sgc["s_arrival_time"])
    # sgc = sgc[["station","p_arrival_time",
    #                     "s_arrival_time"]]

    # sgc_ev = "/home/emmanuel/Tesis/auto/sgc/events/events_1d/csv/CM/2019/358/sgc_events.csv"
    # sgc_ev = pd.read_csv(sgc_ev)
    # sgc_ev["time_pick_p"] = pd.to_datetime(sgc_ev["time_pick_p"])
    # sgc_ev["time_pick_s"] = pd.to_datetime(sgc_ev["time_pick_s"])
    # sgc_ev = sgc_ev[["station","pick_p","time_pick_p",
    #                     "pick_s","time_pick_s"]]

    # fig = plot_asso_eqt(coord,csv_id,dataless,sgc,sgc_ev,starttime,endtime)
    # fig.savefig("/home/emmanuel/tesis_analisis/tesis_analisis/figures/asso/sgc_asso.png",dpi=300)

# # ### eqt-nlloc
    eqt = "/home/emmanuel/EDCT/SeisMonitor/test/plot/asso/seismonitor_picks.csv"
    eqt = pd.read_csv(eqt)
    eqt["p_arrival_time"] = pd.to_datetime(eqt[eqt["phasehint"]=="P"]["arrival_time"])
    eqt["s_arrival_time"] = pd.to_datetime(eqt[eqt["phasehint"]=="S"]["arrival_time"])
    eqt = eqt[["station","p_arrival_time",
                        "s_arrival_time"]]

    eqt_ev = "/home/emmanuel/EDCT/SeisMonitor/test/plot/asso/20191224T000000__20191225T000000.csv"
    eqt_ev = pd.read_csv(eqt_ev)
    eqt_ev["time_pick_p"] = pd.to_datetime(eqt_ev["time_pick_p"])
    eqt_ev["time_pick_s"] = pd.to_datetime(eqt_ev["time_pick_s"])
    eqt_ev = eqt_ev[["station","pick_p","time_pick_p",
                        "pick_s","time_pick_s"]]

    fig = plot_asso_eqt(coord,csv_id,dataless,eqt,eqt_ev,starttime,endtime)
    fig.savefig("/home/emmanuel/EDCT/SeisMonitor/test/plot/asso/gamma_asso.png",dpi=300)