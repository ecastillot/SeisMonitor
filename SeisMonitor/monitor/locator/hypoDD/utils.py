import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from datetime import datetime
import numpy as np
import pandas as pd
from tqdm import tqdm
from glob import glob
from multiprocessing import Pool
import time
import SeisMonitor.utils as ut

from obspy import read_inventory
from math import cos
from numpy import deg2rad

"""
Some utils functions are taken from: https://github.com/wayneweiqiang/QuakeFlow/blob/master/HypoDD/gamma2hypodd.py
such us:
    - download_hypodd()
"""

CORE_HYPODD = os.path.join(os.path.dirname(__file__),"core")

DEG2KM = 111.2

PHA1 = ('# {year:4d} {month:2d} {day:2d} {hour:2d}'
        ' {min:2d} {sec:5.2f}  '
        '{lat:7.4f} {lon:9.4f}   '
        '{dep:5.2f} {mag:5.2f} {he:5.2f} {ve:5.2f} {rms:5.2f} {evid:>9}\n')
PHA2 = '{p.waveform_id.station_code:<5s}  {relt:.4f}  {weight:5.4f}  {p.phase_hint}\n'


def _map_eventid(evid, eventid_map, used_ids, counter):
    """
    Taken from https://docs.obspy.org/_modules/obspy/io/hypodd/pha.html#_write_pha
    """
    idpha = evid
    if evid in eventid_map:
        idpha = eventid_map[evid]
        if not idpha.isdigit() or len(idpha) > 9:
            msg = ('Invalid value in eventid_map, pha event id has to be '
                   'digit with max 9 digits')
            raise ValueError(msg)
        return idpha
    if not idpha.isdigit():
        idpha = ''.join(char for char in idpha if char.isdigit())
    if len(idpha) > 9:
        idpha = idpha[-9:]
    while idpha == '' or idpha in used_ids:
        idpha = str(counter[0])
        counter[0] += 1
    if idpha != evid:
        eventid_map[evid] = idpha
    used_ids.add(idpha)
    return idpha

def write_pha(catalog, filename, eventid_map=None,
               **kwargs):  # @UnusedVariable
    """
    Taken from https://docs.obspy.org/_modules/obspy/io/hypodd/pha.html#_write_pha

    Write a HypoDD PHA file.

    .. warning::
        This function should NOT be called directly, it registers via the
        the :meth:`~obspy.core.event.Catalog.write` method of an
        ObsPy :class:`~obspy.core.event.Catalog` object, call this instead.

    :type catalog: :class:`~obspy.core.event.catalog.Catalog`
    :param catalog: The ObsPy Catalog object to write.
    :type filename: str or file-like object
    :param filename: Filename to write or open file-like object.
    :param dict eventid_map: Desired mapping of event resource ids (dict keys)
        to hypodd event ids (dict values).
        HYPODD expects integer event ids with maximal 9 digits. If the event
        resource id is not present in the mapping,
        the event resource id is stripped of all non-digit characters and
        truncated to a length of 9 chars. If this method does not generate a
        valid hypodd event id, a counter starting at 1000 is used.

    :returns: Dictionary eventid_map with mapping of event resource id to
        hypodd event id. Items are only present if both ids are different.
    """
    if len(catalog) >= 10**10:
        warn('Writing a very large catalog will use event ids that might not '
             'be readable by HypoDD.')
    lines = []
    if eventid_map is None:
        eventid_map = {}
    args_map_eventid = (eventid_map, set(eventid_map.values()), [1])
    for event in catalog:
        try:
            ori = event.preferred_origin() or event.origins[0]
        except IndexError:
            warn(f'Skipping writing event with missing origin: {event}')
            continue
        try:
            mag = event.preferred_magnitude() or event.magnitudes[0]
        except IndexError:
            warn('Missing magnitude will be set to 0.0')
            mag = 0.
        else:
            mag = mag.mag
        evid = event.resource_id.id
        evid = evid.split('/')[-1] if '/' in evid else evid
        evid = _map_eventid(evid, *args_map_eventid)
        rms = (ori.quality.standard_error if 'quality' in ori and ori.quality
               else None)
        rms = rms if rms is not None else 0.0
        he1 = ori.latitude_errors.uncertainty if ori.latitude_errors else None
        he2 = (ori.longitude_errors.uncertainty if ori.longitude_errors
               else None)
        shortening = cos(deg2rad(ori.latitude))
        he = max(0. if he1 is None else he1 * DEG2KM,
                 0. if he2 is None else he2 * DEG2KM * shortening)
        ve = ori.depth_errors.uncertainty if ori.depth_errors else None
        ve = 0. if ve is None else ve / 1000

        year, month, day, hour, min, sec = (
            ori.time.year,
            ori.time.month,
            ori.time.day,
            ori.time.hour,
            ori.time.minute,
            float(ori.time.strftime("%S.%f")),
        )
        lat = ori.latitude
        lon = ori.longitude
        line = PHA1.format(year=year,month=month,day=day,hour=hour,
                            min=min,sec=sec,lat=lat,lon=lon,
                            dep=ori.depth / 1000, 
                            mag=mag,
                           he=he, ve=ve, rms=rms,
                           evid=evid)
        lines.append(line)
        weights = {str(arrival.pick_id): arrival.time_weight
                   for arrival in ori.arrivals if arrival.time_weight}
        for pick in event.picks:
            weight = weights.get(str(pick.resource_id), 1.)
            line = PHA2.format(p=pick, relt=pick.time - ori.time,
                               weight=weight)
            lines.append(line)
    data = ''.join(lines)
    try:
        with open(filename, 'w') as fh:
            fh.write(data)
    except TypeError:
        filename.write(data)
    return None if len(eventid_map) == 0 else eventid_map

def resp2df(resp):
    """
    Parameters:
    -----------
    resp: str
        RESP filepath

    Returns: DataFrame
        Dataframe with the next columns
        network,station,latitude,longitude,elevation
    """
    networks = []
    stations = []
    longitudes = []
    latitudes = []
    elevations = []
    inv = read_inventory(resp)
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
            elevations.append(sta.elevation)
            stations.append(sta.code)
            networks.append(net.code)

    df = {"network":networks,"station":stations,
        "latitude":latitudes,"longitude":longitudes,
        "elevation":elevations}
    df = pd.DataFrame(df)
    # print(df)
    return df

def write_hypoDDstation(df,out_folder):
    """
    Parameters:
    ----------
    df: DataFrame
        Dataframe with the next columns
        network,station,latitude,longitude,elevation.
        Review resp2df function.

    Returns: 
    --------
    msgs: list
        List of each message by station.
        The station message contains the next information:
        
        station lat lon elevation
        ---example---
        CA01A 353.16N 7341.07W 436

    """

    stations_hypodd_inp_path = os.path.join(out_folder, "stations_hypoDD.dat")
    if os.path.exists(stations_hypodd_inp_path):
        print("stations_hypoDD.inp input file already exists.")
        return

    df = df.drop_duplicates("station")
    df = df.sort_values(by="station")
    print(df)
    msgs = []
    for i,row in df.iterrows(): 

        line_hypoDD = f"{row.station} {row.latitude:.3f} {row.longitude:.3f}"
        msgs.append(line_hypoDD)
        # print(msg)
    msgs = "\n".join(msgs)
    with open(stations_hypodd_inp_path, "w") as open_file:
        open_file.write(msgs)
        print("Created stations_hypoDD.dat input file.")
    return msgs

class HypoDDException(Exception):
    """
    Taken from https://github.com/krischer/hypoDDpy/blob/master/hypoddpy/hypodd_relocator.py
    """
    pass

def write_hypoDD_inp_file(vel_model,out_folder,**kwargs):
        """
        Modified from from https://github.com/krischer/hypoDDpy/blob/master/hypoddpy/hypodd_relocator.py
        Writes the hypoDD.inp file.
        """
        hypodd_inp_path = os.path.join(out_folder, "hypoDD.inp")
        if os.path.exists(hypodd_inp_path):
            print("hypoDD.inp input file already exists.")
            return
        # Use this way of defining the string to avoid leading whitespaces.
        hypodd_inp = "\n".join(
            [   "*--- INPUT file selection",
                "* cross correlation diff times:",
                "",
                "*catalog P diff times",
                "dt.ct",
                "* event file:",
                "event.sel",
                "* station file:",
                "stations_hypoDD.dat",
                "*--- OUTPUT file selection",
                "* original locations:",
                "hypoDD.loc",
                "* relocations:",
                "hypoDD.reloc",
                "* station information:",
                "hypoDD.sta",
                "* residual information:",
                "hypoDD.res",
                "* source paramater information:",
                "hypoDD.src",
                "*--- DATA type selection",
                "* IDAT:  0 = synthetics; 1= cross corr; 2= catalog; 3= cross & cat ",
                "* IPHA: 1= P; 2= S; 3= P&S ",
                "* DIST:max dist [km] between cluster centroid and station ",
                "* IDAT   IPHA   DIST",
                "    {IDAT}     {IPHA}     {DIST}",
                "*--- event clustering:",
                "* OBSCC:    min # of obs/pair for crosstime data (0= no clustering)",
                "* OBSCT:    min # of obs/pair for network data (0= no clustering)",
                "* OBSCC  OBSCT    ",
                "     {OBSCC}     {OBSCT}",
                "*--- solution control:",
                "* ISTART:  	1 = from single source; 2 = from network sources",
                "* ISOLV:	1 = SVD, 2=lsqr",
                "* NSET:      	number of sets of iteration with specifications following",
                "*  ISTART  ISOLV  NSET",
                "    {ISTART}        {ISOLV}      {NSET}",
                "*--- data weighting and re-weighting: ",
                "* NITER: 		last iteration to used the following weights",
                "* WTCCP, WTCCS:		weight cross P, S ",
                "* WTCTP, WTCTS:		weight catalog P, S ",
                "* WRCC, WRCT:		residual threshold in sec for cross, catalog data ",
                "* WDCC, WDCT:  		max dist [km] between cross, catalog linked pairs",
                "* DAMP:    		damping (for lsqr only) ",
                "*       ---  CROSS DATA ----- ----CATALOG DATA ----",
                "* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP",
                "{DATA_WEIGHTING_AND_REWEIGHTING}",
                "*--- 1D model:",
                "* NLAY:		number of model layers  ",
                "* RATIO:	vp/vs ratio ",
                "* TOP:		depths of top of layer (km) ",
                "* VEL: 		layer velocities (km/s)",
                "{FORWARD_MODEL}",
                "*--- event selection:",
                "* CID: 	cluster to be relocated (0 = all)",
                "* ID:	cuspids of event to be relocated (8 per line)",
                "* CID    ",
                "    {CID}",
                "* ID    ",
                "{ID}",
            ]
        )
        # Determine all the values.
        values = {}
        # Default 2 
        values["IDAT"] = kwargs.get("IDAT",2)
        # Default 3
        values["IPHA"] = kwargs.get("IDAT",3)
        # Max distance between centroid of event cluster and stations.
        values["DIST"] = kwargs.get("DIST",300)
        # Always set it to 8.
        values["OBSCC"] = kwargs.get("OBSCC",0)
        # If IDAT=3, the sum of OBSCC and OBSCT is taken for both.
        values["OBSCT"] = kwargs.get("OBSCT",4)
        # Start from catalog locations
        values["ISTART"] = kwargs.get("ISTART",2)
        # Least squares solution via conjugate gradients
        values["ISOLV"] = kwargs.get("ISOLV",2)
        # Create the data_weighting and reweightig scheme. Currently static.
        iterations = [
            "   4     -9     -9   -9    -9   1     1      8   -9  70",
            "   4     -9     -9   -9    -9   1     1      6    4  70",
            "   4     -9     -9   -9    -9   1    0.8     4    2  70",
            "   4     -9     -9   -9    -9   1    0.8     3    2  70",
        ]

        values["NSET"] = kwargs.get("NSET",len(iterations))
        values["DATA_WEIGHTING_AND_REWEIGHTING"] = kwargs.get("DATA_WEIGHTING_AND_REWEIGHTING","\n".join(iterations)) 
        # values["DATA_WEIGHTING_AND_REWEIGHTING"] = "\n".join(iterations)
        values["FORWARD_MODEL"] = vel_model
        # Allow relocating of all clusters.
        values["CID"] = kwargs.get("CID",0)
        # Also of all events.
        values["ID"] = kwargs.get("ID","")
        hypodd_inp = hypodd_inp.format(**values)
        with open(hypodd_inp_path, "w") as open_file:
            open_file.write(hypodd_inp)
        print("Created hypoDD.inp input file.")

def setup_velocity_model(model_type, **kwargs):
        """
        Modified from https://github.com/krischer/hypoDDpy/blob/master/hypoddpy/hypodd_relocator.py

        Defines the used velocity model for the forward simulation. The chosen
        model_type determines the available kwargs.

        Possible model_types:

        * layered_p_velocity_with_constant_vp_vs_ratio

        :param vp_vs_ratio: The vp/vs ratio for all layers
        :param layer_tops: List of (depth_of_top_of_layer, layer_velocity_km_s)
           e.g. to define five layers:
            [(0.0, 3.77), (1.0, 4.64), (3.0, 5.34), (6.0, 5.75), (14.0, 6.0)]
            
        * layered_variable_vp_vs_ratio (IMOD 1 in hypodd2.1)
        
        :param layer_tops: List of (depth_of_top_of_layer, layer_velocity_km_s,
        layer_ratios)
            e.g. to define five layers:
            [(-3.0, 3.42, 2.38), 
             (0.5, 4.6, 1.75), 
             (1.0, 5.42, 1.74), 
             (1.5, 5.52, 1.75),
             (2.13, 5.67, 1.77)]
        """
        if model_type == "layered_p_velocity_with_constant_vp_vs_ratio":
            # Check the kwargs.
            if not "layer_tops" in kwargs:
                msg = "layer_tops need to be defined"
                raise HypoDDException(msg)
            if not "vp_vs_ratio" in kwargs:
                msg = "vp_vs_ratio need to be defined"
                raise HypoDDException(msg)
            ratio = float(kwargs.get("vp_vs_ratio"))
            layers = kwargs.get("layer_tops")
            if len(layers) > 30:
                msg = "Model must have <= 30 layers"
                raise HypoDDException(msg)
            depths = [str(_i[0]) for _i in layers]
            velocities = [str(_i[1]) for _i in layers]
            # Use imod 5 which allows for negative station elevations by using
            # straight rays.
            forward_model = [
                # "0",  # IMOD
                "  %i     %.2f" % (len(layers), ratio),
                " ".join(depths),
                " ".join(velocities),
            ]
            # forward_model = [ \
            # "0",  # IMOD
            ## If IMOD=0, number of layers and v_p/v_s ration.
            # "{layer_count} {ratio}".format(layer_count=len(layers),
            # ratio=ratio),
            ## Depth of the layer tops.
            # " ".join(depths),
            ## P wave velocity of layers.
            # " ".join(velocities)]
            forward_model_string = "\n".join(forward_model)
            
        elif model_type == "layered_variable_vp_vs_ratio":
            """
            *--- 1D model, variable  vp/vs ratio:
            * TOP:          depths of top of layer (km)
            * VEL:          layer velocities (km/s) end w/ -9
            * RATIO:        layer ratios  end w/ -9
            * IMOD: 1
            """
            if not "layer_tops" in kwargs:
                msg = "layer_tops need to be defined"
                raise HypoDDException(msg)
            layers = kwargs.get("layer_tops")
            if len(layers) > 30:
                msg = "Model must have <= 30 layers"
                raise HypoDDException(msg)
            depths = [str(_i[0]) for _i in layers]
            velocities = [str(_i[1]) for _i in layers]
            ratios = [str(_i[2]) for _i in layers]
            depths.append('-9')
            velocities.append('-9')
            ratios.append('-9')
            # Use imod 5 which allows for negative station elevations by using
            # straight rays.
            forward_model = [
                # "1",  # IMOD
                " ".join(depths),
                " ".join(velocities),
                " ".join(velocities)
            ]
            forward_model_string = "\n".join(forward_model)
        else:
            msg = "Model type {model_type} unknown."
            msg.format(model_type=model_type)
            raise HypoDDException(msg)

        return forward_model_string

def write_ph2dt_inp_file(out_folder,**kwargs):
        """
        Modified from from https://github.com/krischer/hypoDDpy/blob/master/hypoddpy/hypodd_relocator.py
        Create the ph2dt.inp file.
        """
        # MAXDIST is reused in the hypoDD.inp file. It always needs to be
        # calculated. Fake a forced configuration
        # value.
        ph2dt_inp_file = os.path.join(out_folder, "ph2dt.inp")
        if os.path.exists(ph2dt_inp_file):
            print("ph2dt.inp input file already exists.")
            return
        # Determine the necessary variables. See the documentation of the
        # set_forced_configuration_value method for the reasoning.
        values = {}
        values["MINWGHT"] = kwargs.get("MINWGHT",0.0)
        values["MAXDIST"] = kwargs.get("MAXDIST",300)
        values["MAXSEP"] = kwargs.get("MAXSEP",70)
        values["MAXNGH"] = kwargs.get("MAXNGH",50)
        values["MINLNK"] = kwargs.get("MINLNK",4)
        values["MINOBS"] = kwargs.get("MINOBS",4)
        values["MAXOBS"] = kwargs.get("MAXOBS",100)

        # Use this construction to get rid of leading whitespaces.
        ph2dt_string = [
            "* ph2dt.inp - input control file for program ph2dt",
            "* Input station file:",
            "stations_hypoDD.dat",
            "* Input phase file:",
            "hypoDD.pha",
            "*MINWGHT: min. pick weight allowed [0]",
            "*MAXDIST: max. distance in km between event pair and stations [200]",
            "*MAXSEP: max. hypocentral separation in km [10]",
            "*MAXNGH: max. number of neighbors per event [10]",
            "*MINLNK: min. number of links required to define a neighbor [8]",
            "*MINOBS: min. number of links per pair saved [8]",
            "*MAXOBS: max. number of links per pair saved [20]",
            "*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS",
            "   {MINWGHT}      {MAXDIST}     {MAXSEP}     {MAXNGH}     {MINOBS}      {MINLNK}     {MAXOBS}"
        ]
        ph2dt_string = "\n".join(ph2dt_string)
        ph2dt_string = ph2dt_string.format(**values)
        with open(ph2dt_inp_file, "w") as open_file:
            open_file.write(ph2dt_string)
        print("Writing ph2dt.inp successful")

def download_hypodd():
    '''
    HypoDD can be downloaded from https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Helpful compiling flags: FFLAGS = -O -I${INCLDIR} -mcmodel=large
    '''
    os.system("wget -O HYPODD_1.3.tar.gz http://www.ldeo.columbia.edu/~felixw/HYPODD/HYPODD_1.3.tar.gz")
    os.system("tar -xf HYPODD_1.3.tar.gz")
    os.system("ln -s $(which gfortran) f77")
    os.system("ln -s $(which gfortran) g77")
    os.environ['PATH'] += os.pathsep + os.getcwd()
    os.system("make -C HYPODD/src")

def _download_hypodd():
    '''
    HypoDD can be downloaded from https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Helpful compiling flags: FFLAGS = -O -I${INCLDIR} -mcmodel=large
    '''
    name = "HYPODD_1.3.tar.gz"
    gz_path = os.path.join(CORE_HYPODD,name)
    hypoDD_path = os.path.join(CORE_HYPODD,"HYPODD")
    make_path = os.path.join(hypoDD_path,"src")
    f77_path = os.path.join(CORE_HYPODD,"f77")
    g77_path = os.path.join(CORE_HYPODD,"g77")

    isfile = ut.isfile(gz_path)
    if not isfile:
        os.system(f"wget -O {gz_path} http://www.ldeo.columbia.edu/~felixw/HYPODD/HYPODD_1.3.tar.gz")
    
    if not os.path.isdir(hypoDD_path):
        os.system(f"tar -xf {gz_path} -C {CORE_HYPODD}")

    isfile = ut.isfile(f77_path)
    if not isfile:
        os.system(f"ln -s $(which gfortran) {f77_path}")

    isfile = ut.isfile(g77_path)
    if not isfile:
        os.system(f"ln -s $(which gfortran) {g77_path}")

    print(os.pathsep + CORE_HYPODD)
    os.environ['PATH'] += os.pathsep + CORE_HYPODD


    # print(os.environ['PATH'])
    # print(f"make -C {make_path}")
    os.system(f"make -C HYPODD/src")
    # os.system(f"make -C {make_path}")

def get_vel_layers(df):
    depths = df["depth"].to_list()
    vp = df["vp"].to_list()

    vel_layers = list(zip(depths,vp))
    return vel_layers



if __name__ == "__main__":
    download_hypodd()
    # events = "/home/emmanuel/EDCT/test/magnitude/Ml_magnitude.xml"
    # hypodd = "/home/emmanuel/EDCT/test/hypoDD/hypodd.pha"

    # xml_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
    # df = resp2df(xml_path)
    # sta2hypoDDstation(df)
    # catalog = read_events(events)
    # write_pha(catalog,hypodd)


    # catalog.write(hypodd,format="HYPODDPHA" )
    # print(catalog)

    ####

    out_folder = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/hypoDD"
    vel_model = setup_velocity_model("layered_p_velocity_with_constant_vp_vs_ratio",vp_vs_ratio=1.82,
                        layer_tops=[(0.0,4.80),(4.0,6.60),(25.0,7.00),
                        (32.0,8.00),(40.0,8.10),(100.0,8.20),
                            (200.0,8.30)] 
                        )
    write_hypoDD_inp_file(vel_model,out_folder)
    write_ph2dt_inp_file(out_folder)
    # print(vel_model)