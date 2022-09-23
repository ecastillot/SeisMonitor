import pygmt
import os
import string
import pandas as pd
import numpy as np
import utils as ut
import datetime as dt
from obspy.geodetics.base import gps2dist_azimuth
pygmt.config(FORMAT_GEO_MAP="ddd.xx")
# events = "/home/emmanuel/Tesis/data/map/merge_eqt.csv"
import geopandas as gpd

def rmap_with_profile(catalog,points,
                    region,
                    width = [-15,15],
                    depth=[-3,200],
                    between=None,
                    color=None,
                    truncate=[0,200],
                    cpt_kwargs = {"cmap":'rainbow', "reverse":True,
                                "series":[0, 200] }
                    ):
    
    events = ut.get_df(catalog,
                    between=between,
                    columns =["longitude","latitude","depth",
                                "time","magnitude"])
    
    if truncate:

        bottom_mask = events.depth < truncate[0]
        events.loc[bottom_mask,'depth'] = truncate[0]

        top_mask = events.depth > truncate[1]
        events.loc[top_mask,'depth'] = truncate[1]

    pygmt.makecpt(**cpt_kwargs)


    xfig = pygmt.Figure()
    xfig.coast(
        region=region,
        projection='M6i',
        borders='1/1p,black',
        frame="afg",
        water='lightblue',
        land='grey',
        shorelines=True,
    )

    with pygmt.config(FONT_TITLE=5):
        # fig.basemap(rose="jTL+w1.3c+lO,E,S,N+o-0.1c/3c", map_scale="jBL+w50k+o0.5c/0.5c+f")
        xfig.basemap(map_scale="jBR+w500k+o0.5c/0.5c+f+lkm+at")

    xfig.plot(
            x=events.longitude,
            y=events.latitude,
            # sizes=0.1*events.magnitude,
            # sizes=0.2,
            color=events.depth,
            cmap=True,
            style="c0.1c",
            pen="black",
        )
    xfig.colorbar(frame='af+l"Depth (km)"')

    abc = list(string.ascii_lowercase)
    for k,point in enumerate(points):
        C,D = point
        xfig.plot(x=[C[0], D[0]], y=[C[1], D[1]], projection="M", pen="1p,white")
        xfig.text(x=C[0], y=C[1], text=f"{abc[k]}1", font="6p,Helvetica,black",
                  fill="white",offset="-0.15c/0.05c")
        xfig.text(x=D[0], y=D[1], text=f"{abc[k]}2", font="6p,Helvetica,black",
                 fill="white",offset="0.15c/0.05c")

    ##Mundo
    xfig.shift_origin(xshift='0i',yshift='0i')  # Shift for next call
    proj = 'G-70/0/1.5i'
    xfig.grdimage(
            '@earth_relief_10m',
            region='g',
            projection=proj,
            cmap='globe',
            shading=True,
            )
    xfig.coast(
            region='g',
            projection=proj,
            shorelines=True,
            water='white',
            borders='1/1p,black',
            land='grey',
            frame=True,
        )

    x_reg = [region[0], region[1], region[1], region[0], region[0]]
    y_reg = [region[2], region[2], region[3], region[3], region[2]]
    xfig.plot(x=x_reg,y=y_reg,
        pen="2p,red")

    xfig.show()


def map(catalog,points,
                    region,
                    between=None,
                    truncate=[0,200],
                    cpt_kwargs = {"cmap":'rainbow', "reverse":True,
                                "series":[0, 200] },
                    relief=False,
                    show = True,
                    save = None
                    ):
    shp = "/home/emmanuel/Tesis/data/map/qgis_map/Fallas.shp"
    all_data = gpd.read_file(shp)
    algeciras = all_data[all_data["algeciras"]==1]
    all_data = all_data[all_data["algeciras"]!=1]
    # print(all_data)
    # exit()

    events = ut.get_df(catalog,
                    between=between,
                    columns =["longitude","latitude","depth",
                                "time","magnitude"])
    events.loc[(events["depth"]<0),"depth"] = 0
    events.loc[(events["depth"]>200),"depth"] = 200
    events_sup = events[events["depth"]<=40]
    events_dep= events[events["depth"]>40]
    
    if truncate:

        bottom_mask = events.depth < truncate[0]
        events.loc[bottom_mask,'depth'] = truncate[0]

        top_mask = events.depth > truncate[1]
        events.loc[top_mask,'depth'] = truncate[1]

    pygmt.makecpt(**cpt_kwargs)


    xfig = pygmt.Figure()
    # xfig.coast(
    #     region=region,
    #     projection='M6i',
    #     borders='1/1p,black',
    #     frame=["afg","ES"],
    #     water='lightblue',
    #     land='grey',
    #     shorelines=True,
    # )

    # with pygmt.config(FONT_TITLE=5):
    #     # fig.basemap(rose="jTL+w1.3c+lO,E,S,N+o-0.1c/3c", map_scale="jBL+w50k+o0.5c/0.5c+f")
    #     xfig.basemap(map_scale="jBL+w500k+o0.5c/0.5c+f+lkm+at")

    # xfig.plot(
    #         x=events.longitude,
    #         y=events.latitude,
    #         # sizes=0.1*events.magnitude,
    #         # sizes=0.2,
    #         color=events.depth,
    #         cmap=True,
    #         style="c0.1c",
    #         pen="black",
    #     )
    

    # abc = list(string.ascii_lowercase)
    # for k,point in enumerate(points):
    #     C,D = point
    #     xfig.plot(x=[C[0], D[0]], y=[C[1], D[1]], projection="M", pen="1p,white")
    #     xfig.text(x=C[0], y=C[1], text=f"{abc[k]}1", font="6p,Helvetica,black",
    #               fill="white",offset="-0.15c/0.05c")
    #     xfig.text(x=D[0], y=D[1], text=f"{abc[k]}2", font="6p,Helvetica,black",
    #              fill="white",offset="0.15c/0.05c")

    region2 = np.array(points)
    region2 = region2.reshape(region2.shape[0]*region2.shape[1],-1)
    latm = min(region2[:,1])
    latM = max(region2[:,1])
    lonm = min(region2[:,0])
    lonM = max(region2[:,0])
    region2 = [lonm-0.5, lonM+0.5 , latm-0.5, latM+0.5]

    x_reg = [region2[0], region2[1], region2[1], region2[0], region2[0]]
    y_reg = [region2[2], region2[2], region2[3], region2[3], region2[2]]
    # xfig.plot(x=x_reg,y=y_reg,
    #         pen="2p,red,-")

    # print(region2)
    # xfig.shift_origin(xshift='0c',yshift='7.2c')  # Shift for next call

    if relief:
        # pygmt.makecpt(cmap="grayC", series=[200, 4000, 10])
        grid = pygmt.datasets.load_earth_relief(resolution="03s", 
                                                    region=region2)
        xfig.grdimage(
            '@earth_relief_03s',
            region=region2,
            projection='M4i',
            cmap="gray",
            shading=True,
            frame=["afg","WNse"]
        )
        pygmt.makecpt(**cpt_kwargs)
    else:
        xfig.coast(
            region=region2,
            projection='M4i',
            shorelines=True,
            water='lightblue',
            land='grey',
            frame=["afg","WNse"],
        )
    # xfig.plot(data=algeciras,color="white", pen=["0.02c,green"],connection="rr",label="Algeciras")
    # xfig.plot(data=all_data,color="white", pen=["0.02c,black,-"],connection="rr")
    
    
    xfig.plot(
        x=events_sup.longitude,
        y=events_sup.latitude,
        # sizes=0.1*events.magnitude,
        # sizes=0.2,
        color=events_sup.depth,
        cmap=True,
        style="c0.07c",
    )
    xfig.plot(
        x=events_dep.longitude,
        y=events_dep.latitude,
        # sizes=0.1*events.magnitude,
        # sizes=0.2,
        color=events_dep.depth,
        cmap=True,
        style="c0.1c",
        pen="black",
    )
    abc = list(string.ascii_lowercase)
    for k,point in enumerate(points):
        C,D = point
        xfig.plot(x=[C[0], D[0]], y=[C[1], D[1]], projection="M", pen="2p,white")
        xfig.text(x=C[0], y=C[1], text=f"{abc[k]}'1", font="10p,Helvetica,black",
                  fill="white",offset="-0.15c/0.05c")
        xfig.text(x=D[0], y=D[1], text=f"{abc[k]}'2", font="10p,Helvetica,black",
                 fill="white",offset="0.15c/0.05c")
    # xfig.legend()
    xfig.colorbar(frame='af+l"Depth (km)"')
    with pygmt.config(FONT_TITLE=5):
        # fig.basemap(rose="jTL+w1.3c+lO,E,S,N+o-0.1c/3c", map_scale="jBL+w50k+o0.5c/0.5c+f")
        xfig.basemap(map_scale="jBL+w100k+o0.5c/0.5c+f+lkm+at")
    if show:
        xfig.show()
    if save != None:
        xfig.savefig(save)

    return xfig,events


def profile(catalog,points,width = [-15,15],
            depth=[-3,200],color=None,
            replace_nan_magnitude=1.5,
            figsize = ("30c", "36c"),
            show = True,
            save = None):

    if replace_nan_magnitude != None:
        catalog["magnitude"] = catalog["magnitude"].fillna(replace_nan_magnitude)

    fig = pygmt.Figure()

    n = len(points)
    n_square = np.sqrt(n)
    c = int(n_square )
    r = int(n_square )
    if n%n_square != 0:
        c += 1
        r += 1


    abc = list(string.ascii_lowercase)
    with pygmt.clib.Session() as session:
        session.call_module('gmtset', 'FONT 10p')
        subplot =   fig.subplot(
                        nrows=r, ncols=c, 
                        figsize=figsize, 
                        # autolabel=True,
                        # autolabel="+jLT+gwhite",
                        autolabel="(a)+jBL+o-0c/.25c",
                        # autolabel=True,
                        frame=[f'xafg+l"Distance (km)"', 
                                f'yafg{depth[1]}+l"Depth (km)"',
                                "WSen"],
                        sharex="b",  # shared x-axis on the bottom side
                        sharey="l",  # shared y-axis on the left side
                        margins=["0.1c", "0.1c"]
                    ) 
    with subplot:
        for f,point in enumerate(points):
            C,D = point

            r,a,ba = gps2dist_azimuth(C[1],C[0],D[1],D[0])
            corte_region= [0, r/1e3] + depth

            fig.basemap(region=corte_region, 
                        projection="X?/-?", 
                        panel=f)
            with fig.set_panel(panel=f,fixedlabel=" "+"'"):

                df = pygmt.project(
                            data=catalog,
                            unit=True,
                            center=C,
                            endpoint=D,
                            convention="pz",
                            width=width,
                            verbose=True,
                        )
                df = df.rename(columns={0:"distance",1:"depth",2:"time",
                                        3:"magnitude"})
                df["magnitude"] = df["magnitude"].astype(float)
                df["time"] = pd.to_datetime(df["time"])
                df = df[["distance","depth","time","magnitude"]]
                df = df.sort_values("magnitude",ascending=False,ignore_index=True)

                # print(df)
                print(f"profile_{f}")
                if color == None:
                    pygmt.makecpt(cmap="gray", series=[df.time.min(), df.time.max()])
                    fig.plot(x=df.distance,y=df.depth,
                        projection="X?/-?", 
                        style="cc",
                        size=0.1 * np.sqrt(1.5 ** (df.magnitude*1.5)),
                        pen=0.11,
                        cmap=True, 
                        color=df.time,
                        )
                else: 
                    fig.plot(x=df.distance,y=df.depth,
                        projection="X?/-?", 
                        style="cc",
                        size=0.1 * np.sqrt(1.5 ** (df.magnitude*1.5)),
                        #  label= str(starttime.hour)+"-"+str(endtime.hour),
                        pen=0.11, 
                        color="gray")
        # cb = fig.colorbar(position="JMR+o1c/0c+w7c/0.5c+mc",frame="afg")
        if color == None:
            cb = fig.colorbar(frame=["afg","y+lHour"],
                            # position="JBC+h",
                            position="JMR+o1c/0c+w7c/0.5c",
                            scale=1,box=True)

    if show:
        fig.show()
    if save != None:
        fig.savefig(save)

    return fig

def map_with_profile(catalog,points,
                    region,
                    width = [-15,15],
                    depth=[-3,200],
                    between=None,
                    color=None,
                    truncate=[0,200],
                    cpt_kwargs = {"cmap":'rainbow', "reverse":True,
                                "series":[0, 200] },
                    profile_figsize=("30c", "36c"),
                    show = True,
                    save = None
                    ):
    if save != None:
        if not os.path.isdir(save[1]):
            os.makedirs(save[1])

        map_path = os.path.join(save[1],f"map_{save[0]}.png")
        profile_path = os.path.join(save[1],f"profile_{save[0]}.png")

    map_fig,events = map(catalog,points,
                    region,
                    between=between,
                    truncate=truncate,
                    cpt_kwargs = cpt_kwargs,
                    show = show,
                    save = map_path
                    )
    profile_fig = profile(events,points,width = width,
            depth=depth,
            figsize=profile_figsize,
            color=color,
            show = show,
            save = profile_path)


if __name__ == "__main__":

    events = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/events_20160101T20220901.csv"
    n = 25
    # p1 = [-73.5698,9.3594]
    p1 = [-77.380,4.301]
    p2 = [-71.317,11.105]
    theta = 90
    l = (500,500)
    region = [-91, -70, 0, 20]
    # width = [-1.85,1.85]
    width = [-23,23]
    # width = [-40,40]
    depth = [-3,200]
    profile_figsize = ("120c", "30c")
    # profile_figsize = None
    name = "depth_cm"
    # name = "mesetas_nodetail_pnet"
    # name = "nido_plate_sgc"
    folder = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/cm"
    folder = os.path.join(folder,name)
    # p1,p2 = ut.get_t_points(p1,p2,l[0]) ## lineas paralelas
    centers,azi = ut.get_centers(n,p1,p2) ##lineas perpendiculares

    # r,a,ba = gps2dist_azimuth(centers[0][1],centers[0][0],
    #                 centers[1][1],centers[1][0])
    # print(r/1e3)
    # exit()

    points = []
    for center in centers:
        xr,yr,xl,yl,alpha = ut.get_line_in_map(center,l,azi+theta,
                                            d2k=114,save=None)
        point_1 = (xl,yl)
        point_2 = (xr,yr)
        points.append((point_1,point_2))

    # p2 = [-73.3138,6.7503]
    # p1 = [-73.2623,6.9633]
    # n=2
    # theta = 90
    # centers,azi = ut.get_centers(n,p1,p2) ##lineas perpendiculares
    # r,a,ba = gps2dist_azimuth(centers[0][1],centers[0][0],
    #                 centers[1][1],centers[1][0])
    # for center in centers:
    #     xr,yr,xl,yl,alpha = ut.get_line_in_map(center,l,azi+theta,
    #                                         d2k=114,save=None)
    #     point_1 = (xl,yl)
    #     point_2 = (xr,yr)
    #     points.append((point_1,point_2))

    # p1 = [-73.2651,7.3209]
    # p2 = [-73.1207,8.8896]
    # n=7
    # theta = 90
    # centers,azi = ut.get_centers(n,p1,p2) ##lineas perpendiculares
    # r,a,ba = gps2dist_azimuth(centers[0][1],centers[0][0],
    #                 centers[1][1],centers[1][0])
    # for center in centers:
    #     xr,yr,xl,yl,alpha = ut.get_line_in_map(center,l,azi+theta,
    #                                         d2k=114,save=None)
    #     point_1 = (xl,yl)
    #     point_2 = (xr,yr)
    #     points.append((point_1,point_2))

    # r,a,ba = gps2dist_azimuth(centers[0][1],centers[0][0],
    #                 centers[1][1],centers[1][0])
    # print(r/1e3)
    # exit()


    # points = [((-74.26008,3.49798),(-74.15192,3.37682))]
    # points = [((-74.7956,4.09803),(-73.9116,2.9957))]

    # points = [(-74.26008,3.49798),(-74.15192,3.37682)]
    # x_coords, y_coords = zip(*points)
    # A = np.vstack([x_coords,np.ones(len(x_coords))]).T
    # m, c = np.linalg.lstsq(A, y_coords)[0]
    # # x1,x2 = -74.7956,-73.8016
    # x1,x2 = -74.5,-73.95
    # y1,y2 = (m*x1+c,m*x2+c)
    # points = [[(x1,y1),(x2,y2)]]

    # print(points)
    # exit()

    starttime = dt.datetime.strptime('20191224 190000', '%Y%m%d %H%M%S')
    endtime = dt.datetime.strptime('20191227 050000', '%Y%m%d %H%M%S')

    save = (name,folder)
    map_with_profile(events,points,
                    region,
                    width = width ,
                    depth=depth,
                    between=None,
                    color="gray",
                    # between=[starttime,endtime],
                    # color=None,
                    truncate=[0,200],
                    cpt_kwargs = {"cmap":'rainbow', "reverse":True,
                                "series":[0, 200] },
                    profile_figsize=profile_figsize,
                    show = True,
                    save = save
                    )

    