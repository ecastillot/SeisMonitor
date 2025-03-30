#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 20:00:00 2020
@author: Emmanuel_Castillo
last update: 14-11-2020 
"""
import os
import json
import time
import logging
import datetime as dt
from SeisMonitor.utils import printlog
from obspy.core.stream import Stream
from obspy.clients.fdsn.mass_downloader import Restrictions
from obspy.core.inventory.inventory import Inventory, read_inventory
from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain


class DownloadRestrictions:
    """Class defining restrictions for downloading seismic data."""
    
    def __init__(
        self,
        mseed_storage,
        chunklength_in_sec=None,
        threshold=60,
        overlap_in_sec=0,
        picker_args={},
        groupby='{network}.{station}.{channel}',
        n_processor=None
    ):
        """Initialize DownloadRestrictions with download parameters.
        
        Args:
            mseed_storage (str): Path template for waveform storage
            chunklength_in_sec (int, optional): Length of each time chunk in seconds
            threshold (int): Minimum length in seconds for download
            overlap_in_sec (int): Overlap between chunks in seconds
            picker_args (dict): Picker parameters (batch_size, overlap, length)
            groupby (str): Grouping pattern for traces
            n_processor (int, optional): Number of parallel processors
        """
        self.mseed_storage = mseed_storage
        self.chunklength_in_sec = chunklength_in_sec
        self.threshold = threshold
        self.overlap_in_sec = overlap_in_sec
        self.picker_args = picker_args
        self.groupby = groupby
        self.n_processor = n_processor


def sanitize_provider_times(providers):
    """Ensure all providers have the same time interval.
    
    Args:
        providers (list): List of provider objects
        
    Returns:
        list: Sanitized providers list
        
    Raises:
        Exception: If providers have different time intervals
    """
    provider_times = [(p.waveform_restrictions.starttime,
                      p.waveform_restrictions.endtime) for p in providers]
    if provider_times.count(provider_times[0]) == len(provider_times):
        return providers
    raise Exception("Providers must have the same interval time")


def get_max_allowed_batch_size(data_length, segment_length, overlap):
    """Calculate maximum allowed batch size for given parameters.
    
    Args:
        data_length (float): Length of data in seconds
        segment_length (float): Length of each batch segment in seconds
        overlap (float): Overlap fraction (0-1)
        
    Returns:
        int: Maximum batch size
    """
    max_batch_size = (data_length - (overlap * segment_length)) / (segment_length * (1 - overlap))
    max_batch_size = int(max_batch_size)
    return max(1, max_batch_size)


def write_stream(st, mseed_storage, threshold=None, picker_args={}, ppc_and_comment=[False, ""]):
    """Write seismic stream to file with restrictions.
    
    Args:
        st (Stream): Obspy Stream object to write
        mseed_storage (str): Path template for storage
        threshold (int, optional): Minimum length threshold in seconds
        picker_args (dict): Picker parameters (batch_size, overlap, length)
        ppc_and_comment (list): [preprocessed_flag, comment]
    """
    ppc, comment = ppc_and_comment
    tr = st[0]
    mseed_filename = get_mseed_filename(mseed_storage, tr, ppc)

    download = True
    if threshold is not None:
        length = abs(tr.stats.endtime - tr.stats.starttime)
        if length < threshold:
            comment = f"length:{length} < threshold:{threshold}"
            printlog("info", "Downloader: False", f"{mseed_filename}  {comment}")
            download = False

    if picker_args:
        overlap = picker_args["overlap"]
        batch_size = picker_args["batch_size"]
        segment_length = picker_args["length"]
        data_length = abs(tr.stats.endtime - tr.stats.starttime)
        max_batch_size = get_max_allowed_batch_size(data_length, segment_length, overlap)
        
        if max_batch_size < batch_size:
            comment = (
                f"This mseed only can be used with {max_batch_size} batches. "
                f"In order to download the data, the batch size must be >= {batch_size}. "
                f"Modify this condition changing 'picker_args':{picker_args}' parameter."
            )
            printlog("info", "Downloader: False", f"{mseed_filename}  {comment}")
            download = False

    if not os.path.isfile(mseed_filename) and download:
        mseed_dir = os.path.dirname(mseed_filename)
        if not os.path.isdir(mseed_dir):
            os.makedirs(mseed_dir)
        st.write(mseed_filename, format="MSEED")
        printlog("info", "Downloader: True", f"{mseed_filename}  {comment}")
    elif os.path.isfile(mseed_filename):
        printlog("info", "Downloader: Exist", f"{mseed_filename}  {comment}")


def get_mseed_filename(_str, tr, ppc=False):
    """Generate MSEED filename from template and trace info.
    
    Args:
        _str (str): Path template with wildcards
        tr (Trace): Obspy Trace object
        ppc (bool): Flag for preprocessed data
        
    Returns:
        str: Generated filename
        
    Raises:
        TypeError: If resulting path is not a string
    """
    strftime = "%Y%m%dT%H%M%SZ"
    params = {
        'network': tr.stats.network,
        'station': tr.stats.station,
        'location': tr.stats.location,
        'channel': tr.stats.channel,
        'starttime': tr.stats.starttime.strftime(strftime),
        'endtime': tr.stats.endtime.strftime(strftime),
        'year': tr.stats.starttime.year,
        'month': tr.stats.starttime.month,
        'day': tr.stats.endtime.day,
        'julday': tr.stats.endtime.julday
    }
    path = _str.format(**params)
    if ppc:
        path += ".ppc"
    
    if not isinstance(path, (str, bytes)):
        raise TypeError(f"'{path}' is not a filepath.")
    return path


def get_chunktimes(starttime, endtime, chunklength_in_sec, overlap_in_sec=0):
    """Generate list of chunk time intervals.
    
    Args:
        starttime (UTCDateTime): Start time
        endtime (UTCDateTime): End time
        chunklength_in_sec (int): Chunk length in seconds
        overlap_in_sec (int): Overlap in seconds
        
    Returns:
        list: List of (start, end) time tuples
        
    Raises:
        Exception: If chunklength_in_sec is 0
    """
    if chunklength_in_sec == 0:
        raise Exception("chunklength_in_sec must be different than 0")
    if chunklength_in_sec is None:
        return [(starttime, endtime)]
    if overlap_in_sec is None:
        overlap_in_sec = 0

    deltat = starttime
    dtt = dt.timedelta(seconds=chunklength_in_sec)
    overlap_dt = dt.timedelta(seconds=overlap_in_sec)
    times = []
    
    while deltat < endtime:
        if deltat + dtt > endtime:
            break
        times.append((deltat, deltat + dtt))
        deltat += dtt - overlap_dt
    
    if deltat < endtime:
        times.append((deltat, endtime))
    return times


def get_st_according2preference(st, location_list, channel_list):
    """Filter stream based on location and channel preferences.
    
    Args:
        st (Stream): Input stream
        location_list (list): Preferred locations in order
        channel_list (list): Preferred channel prefixes in order
        
    Returns:
        Stream: Filtered stream
    """
    logger = logging.getLogger('Downloader: preference')
    if len(st) >= 1000:
        stats = st[0].stats
        printlog(
            "error", "Downloader: preference",
            f"{stats.network}-{stats.station} No filter by preference. gaps: {len(st)}"
        )
        return st

    if not location_list:
        new_st = st
    else:
        locations = [x[2] for x in st._get_common_channels_info().keys()]
        loc_pref = next((loc for loc in location_list if loc in locations), None)
        stats = st[0].stats
        new_st = st if loc_pref is None else st.select(
            network=stats.network, station=stats.station, location=loc_pref
        )

    if channel_list:
        common_channels = [x[3][:2] for x in new_st._get_common_channels_info().keys()]
        cha_pref = next((cha for cha in channel_list if cha in common_channels), None)
        if cha_pref:
            new_st = new_st.select(
                network=stats.network, station=stats.station,
                location=loc_pref, channel=f'{cha_pref}?'
            )

    preference = [(x[2], x[3]) for x in st._get_common_channels_info().keys()]
    selected = [(x[2], x[3]) for x in new_st._get_common_channels_info().keys()]
    printlog(
        "debug", f"Downloader: {stats.network}-{stats.station}: ",
        f"available:{preference}, selected {selected} according to preference"
    )
    return new_st


def get_filenames(mseed, filter_net=[], filter_sta=[], filter_cha=[]):
    """Get MSEED filenames with filtering.
    
    Args:
        mseed (str): Directory path
        filter_net (list): Networks to exclude
        filter_sta (list): Stations to exclude
        filter_cha (list): Channels to include (if specified)
        
    Returns:
        list: Filtered filenames
        
    Raises:
        Exception: If filter parameters are not lists
    """
    for filter_param in (filter_net, filter_sta, filter_cha):
        if not isinstance(filter_param, list):
            raise Exception(f"{filter_param} must be a list")

    filenames = []
    file_list = [f for f in os.listdir(mseed) if f.lower().endswith('.mseed')]
    for filename in file_list:
        net, sta, loc, cha = filename.split('__')[0].split('.')
        if net in filter_net or sta in filter_sta:
            continue
        if filter_cha and cha not in filter_cha:
            continue
        filenames.append(filename)
    return filenames


def get_all_sdswaveforms(client, **kwargs):
    """Get waveforms from client with multiple parameters.
    
    Args:
        client (Client): Obspy client
        ``**kwargs``: network, station, location, channel, starttime, endtime
        
    Returns:
        Stream: Combined waveforms
    """
    args = {k: v.split(",") if k in ("network", "station", "location", "channel") else v
            for k, v in kwargs.items()}
    st = Stream()
    
    for net in args["network"]:
        for sta in args["station"]:
            for loc in args["location"]:
                for cha in args["channel"]:
                    try:
                        myst = client.get_waveforms(
                            net, sta, loc, cha,
                            args['starttime'], args['endtime']
                        )
                        if myst:
                            st += myst
                    except Exception as e:
                        seedname = f"{net}.{sta}.{loc}.{cha}.{args['starttime']}.{args['endtime']}"
                        printlog("error", seedname, str(e))
    return st


def select_inventory(inv, network, station, location, channel, starttime, endtime):
    """Filter inventory based on specified criteria.
    
    Args:
        inv (Inventory): Input inventory
        network (str): Comma-separated network codes
        station (str): Comma-separated station codes
        location (str): Comma-separated location codes
        channel (str): Comma-separated channel codes
        starttime (UTCDateTime): Start time filter
        endtime (UTCDateTime): End time filter
        
    Returns:
        Inventory: Filtered inventory
    """
    networks, stations, locations, channels = (
        network.split(','), station.split(','), 
        location.split(','), channel.split(',')
    )
    
    inventory = inv.select(
        network=networks[0], station=stations[0],
        location=locations[0], channel=channels[0],
        starttime=starttime, endtime=endtime
    )
    
    for net in networks:
        for sta in stations:
            for loc in locations:
                for cha in channels:
                    one_inv = inv.select(
                        network=net, station=sta, location=loc, channel=cha,
                        starttime=starttime, endtime=endtime
                    )
                    inventory = inventory.__add__(one_inv)
    return inventory


def get_client_waveforms(client, bulk, waveform_restrictions, processing):
    """Get and process waveforms from client.
    
    Args:
        client (Client): Obspy client
        bulk (tuple): (net, sta, loc, cha, start, end)
        waveform_restrictions: Waveform restrictions object
        processing: Processing object or None
        
    Returns:
        tuple: (Stream, preprocessed_flag, comment)
    """
    strftime = "%Y%m%dT%H%M%SZ"
    net, sta, loc, cha, starttime, endtime = bulk
    why = "-".join((net, sta, loc, cha, starttime.strftime(strftime), endtime.strftime(strftime)))
    
    try:
        st = client.get_waveforms(net, sta, loc, cha, starttime, endtime)
    except Exception as e:
        printlog("info", "Downloader: False", f"{why}->{e}")
        return Stream(), False, ""
    
    if not st:
        printlog("warning", "Downloader: False", f"0 Trace(s) in Stream: {why}")
        return Stream(), False, ""
    
    st = get_st_according2preference(
        st,
        waveform_restrictions.location_preferences,
        waveform_restrictions.channel_preferences
    )
    
    return processing.run(st) if processing else (st, False, "")


def write_client_waveforms(client, bulk, waveform_restrictions, download_restrictions, processing):
    """Get and write waveforms to file.
    
    Args:
        client (Client): Obspy client
        bulk (tuple): (net, sta, loc, cha, start, end)
        waveform_restrictions: Waveform restrictions object
        download_restrictions (DownloadRestrictions): Download parameters
        processing: Processing object or None
    """
    st, ppc, comment = get_client_waveforms(client, bulk, waveform_restrictions, processing)
    st_dict = st._groupby(download_restrictions.groupby)
    for st in st_dict.values():
        write_stream(
            st,
            download_restrictions.mseed_storage,
            download_restrictions.threshold,
            download_restrictions.picker_args,
            [ppc, comment]
        )


def get_merged_inv_and_json(providers):
    """Merge inventory and JSON info from providers.
    
    Args:
        providers (list): List of provider objects
        
    Returns:
        tuple: (Inventory, dict, list, list) - inventory, JSON info, updated providers, stations outside domains
    """
    json_info = {}
    stations_outside_domains = []
    inventory = Inventory()
    updated_providers = []
    
    for provider in providers:
        client = provider.client
        restrictions = provider.waveform_restrictions
        try:
            inv = (read_inventory(provider.xml) if provider.xml else
                   client.get_stations(
                       network=restrictions.network, station=restrictions.station,
                       location=restrictions.location, channel=restrictions.channel,
                       starttime=restrictions.starttime, endtime=restrictions.endtime,
                       level='channel'
                   ))
            if provider.xml:
                inv = select_inventory(
                    inv, restrictions.network, restrictions.station,
                    restrictions.location, restrictions.channel,
                    restrictions.starttime, restrictions.endtime
                )
        except Exception:
            printlog("error", "Inventory", f"No get_stations with {restrictions.__dict__}")
            inv = Inventory()
        
        inv, jinf, sod = get_inv_and_json(
            inv,
            restrictions.filter_networks,
            restrictions.filter_stations,
            restrictions.filter_domain
        )
        
        bulk_info = list(set([(net.code, sta.code) for net in inv.networks for sta in net.stations]))
        provider.waveform_restrictions.bulk_info = bulk_info
        
        inventory = inventory.__add__(inv)
        json_info.update(jinf)
        stations_outside_domains.append(sod)
        updated_providers.append(provider)
    
    stations_outside_domains = [x for xs in stations_outside_domains for x in xs]
    return inventory, json_info, updated_providers, stations_outside_domains


def inside_the_polygon(p, pol_points):
    """Check if a point is inside a polygon.
    
    Args:
        p (tuple): Point coordinates (lon, lat)
        pol_points (list): List of polygon points (lon, lat)
        
    Returns:
        bool: True if point is inside polygon
    """
    V = tuple(pol_points) + (pol_points[0],)
    cn = 0
    for i in range(len(V) - 1):
        if (V[i][1] <= p[1] < V[i+1][1] or V[i][1] > p[1] >= V[i+1][1]):
            vt = (p[1] - V[i][1]) / float(V[i+1][1] - V[i][1])
            if p[0] < V[i][0] + vt * (V[i+1][0] - V[i][0]):
                cn += 1
    return cn % 2 == 1


def get_inv_and_json(inventory, filter_networks=[], filter_stations=[], filter_domain=[-180, 180, -90, 90]):
    """Generate inventory and JSON info with filtering.
    
    Args:
        inventory (Inventory): Input inventory
        filter_networks (list): Networks to exclude
        filter_stations (list): Stations to exclude
        filter_domain (list): [minlon, maxlon, minlat, maxlat]
        
    Returns:
        tuple: (Inventory, dict, list) - filtered inventory, station info, stations outside domain
    """
    if not filter_domain:
        filter_domain = [-180, 180, -90, 90]
    
    polygon = [
        (filter_domain[0], filter_domain[2]), (filter_domain[0], filter_domain[3]),
        (filter_domain[1], filter_domain[3]), (filter_domain[1], filter_domain[2]),
        (filter_domain[0], filter_domain[2])
    ]
    
    station_list = {}
    logger = logging.getLogger('json')
    toprint = []
    stations_outside_domain = []
    
    for ev in inventory:
        net = ev.code
        if net not in filter_networks:
            for st in ev:
                sta = st.code
                msg = f'{net}-{sta}'
                logger.debug(msg)
                
                if sta not in filter_stations:
                    lon, lat, elv = st.longitude, st.latitude, st.elevation
                    if not inside_the_polygon((lon, lat), polygon):
                        inventory = inventory.remove(network=net, station=sta)
                        stations_outside_domain.append(sta)
                        continue
                    
                    channels = list(set(ch.code for ch in st.channels))
                    if channels and sta not in station_list:
                        sample_rates = [st.select(channel=ch).channels[0].sample_rate for ch in channels]
                        station_list[sta] = {
                            "network": net,
                            "channels": channels,
                            "coords": [lat, lon, elv],
                            "sampling_rate": sample_rates
                        }
                    toprint.append(msg)
                else:
                    inventory = inventory.remove(network=net, station=sta)
        else:
            inventory = inventory.remove(network=net)
    
    logger.info(str(toprint) + ' ok')
    return inventory, station_list, stations_outside_domain

if __name__ == "__main__":
	from obspy.clients.fdsn import Client as FDSN_Client
	from obspy.core.utcdatetime import UTCDateTime
	# from .restrictions import PreprocRestrictions

	# IRIS_client = FDSN_Client(base_url="IRIS", user='gaprietogo@unal.edu.co',password="DaCgmn3hNjg")
	# st = IRIS_client.get_waveforms(network="YU",station="FC04,FC01",location="",channel="*",starttime=UTCDateTime("2016-04-22T00:00:00.0"),endtime=UTCDateTime("2016-04-23T00:00:00.0"))
	# st = get_st_according2preference(st,["","00","20","10"],["HH","BH"])
	# prove = get_all_sdswaveforms(client=IRIS_client,network="YU",station="FC04,FC01",location="",channel="HH?",starttime=UTCDateTime("2016-04-22T00:00:00.0"),endtime=UTCDateTime("2016-04-23T00:00:00.0"))
	# print(prove)


	# client = FDSN_Client('http://sismo.sgc.gov.co:8080')
	# inv = client.get_stations(network="CM",
	# 					  station="BAR2",
	# 					  location="*",
	# 					  channel="*",
	# 					  starttime=UTCDateTime("2019-04-23T00:00:00.0"),
	# 					  endtime=UTCDateTime("2019-04-23T00:02:00.0"))
	# st = client.get_waveforms(network="CM",
	# 					  station="BAR2",
	# 					  location="*",
	# 					  channel="*",
	# 					  starttime=UTCDateTime("2019-04-23T00:00:00.0"),
	# 					  endtime=UTCDateTime("2019-04-23T00:02:00.0"))
	# ppc_restrictions = PreprocRestrictions(["CM.BAR2"],detrend={'type':'simple'})
	# preproc_stream(st,ppc_restrictions)
	######### inventory
	# json_path = "/home/ecastillo/repositories/AIpicker_modules/onejson.json"
	# client_baseurl = "http://sismo.sgc.gov.co:8080"

	# restrictions = DownloadRestrictions(network="CM",
	#					   station="BAR2",
	#					   location="*",
	#					   channel="*",
	#					   starttime=UTCDateTime("2019-04-23T00:22:34.5"),
	#					   endtime=UTCDateTime("2019-04-25T00:23:39.5"),
	#					   chunklength_in_sec=5000,
	#					   overlap_in_sec=None,
	#					   groupby='{network}.{station}.{channel}')
	# xml = "/home/ecastillo/repositories/AIpicker_modules/CM.xml"

	# makeStationList(json_path,client_baseurl,restrictions,from_xml=xml)

	######## get stations
	# json_path = "/home/ecastillo/repositories/AIpicker_modules/onejson.json"
	# client_baseurl = "IRIS"
	# restrictions = DownloadRestrictions(network="CI",
	#				   station="BAK,ARV",
	#				   location="*",
	#				   channel="BH*",
	#				   starttime=UTCDateTime("2020-09-01 00:00:00.00"),
	#				   endtime=UTCDateTime("2020-09-02 00:00:00.00"),
	#				   chunklength_in_sec=3600,
	#				   overlap_in_sec=None,
	#				   groupby='{network}.{station}.{channel}')
	# makeStationList(json_path=json_path,client_baseurl="IRIS",restrictions=restrictions)