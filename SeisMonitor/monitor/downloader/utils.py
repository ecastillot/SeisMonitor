"""
Utility functions for seismic waveform downloading and inventory handling.

This module provides helper classes and functions used for:

* Waveform downloading
* MiniSEED file management
* Inventory filtering
* Station metadata extraction
* Stream preference selection
* Chunked time-window generation

The implementation is designed for ObsPy-based workflows.
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
    """
    Container class defining waveform download restrictions.

    Parameters
    ----------
    mseed_storage : str
        MiniSEED output path template.

    chunklength_in_sec : int, optional
        Length of each download chunk in seconds.

    threshold : int, default=60
        Minimum waveform duration threshold in seconds.

    overlap_in_sec : int, default=0
        Chunk overlap in seconds.

    picker_args : dict, optional
        Picker configuration dictionary.

    groupby : str, default="{network}.{station}.{channel}"
        Grouping rule for traces.

    n_processor : int, optional
        Number of parallel workers.

    Attributes
    ----------
    mseed_storage : str
        MiniSEED output path template.

    chunklength_in_sec : int or None
        Chunk duration in seconds.

    threshold : int
        Minimum waveform duration threshold.

    overlap_in_sec : int
        Overlap between chunks.

    picker_args : dict
        Picker configuration.

    groupby : str
        Trace grouping rule.

    n_processor : int or None
        Number of parallel workers.
    """
    
    def __init__(
        self,
        mseed_storage,
        chunklength_in_sec=None,
        threshold=60,
        overlap_in_sec=0,
        picker_args={},
        groupby='{network}.{station}.{channel}',
        n_processor=None):
        """
        Initialize download restrictions.
        """
        self.mseed_storage = mseed_storage
        self.chunklength_in_sec = chunklength_in_sec
        self.threshold = threshold
        self.overlap_in_sec = overlap_in_sec
        self.picker_args = picker_args
        self.groupby = groupby
        self.n_processor = n_processor


def sanitize_provider_times(providers):
    """
    Ensure all providers share the same time interval.

    Parameters
    ----------
    providers : list
        Provider list.

    Returns
    -------
    list
        Validated provider list.

    Raises
    ------
    ValueError
        If providers use different time intervals.
    """
    provider_times = [(p.waveform_restrictions.starttime,
                      p.waveform_restrictions.endtime) for p in providers]
    if provider_times.count(provider_times[0]) == len(provider_times):
        return providers
    raise Exception("Providers must have the same interval time")


def get_max_allowed_batch_size(data_length, segment_length, overlap):
    """
    Compute the maximum valid batch size.

    Parameters
    ----------
    data_length : float
        Total waveform duration in seconds.
    segment_length : float
        Segment duration in seconds.
    overlap : float
        Fractional overlap between segments.

    Returns
    -------
    int
        Maximum allowed batch size.
    """
    max_batch_size = (data_length - (overlap * segment_length)) / (segment_length * (1 - overlap))
    max_batch_size = int(max_batch_size)
    return max(1, max_batch_size)


def write_stream(st, mseed_storage, threshold=None, picker_args={}, ppc_and_comment=[False, ""]):
    """
    Write a waveform stream to MiniSEED format.

    Parameters
    ----------
    st : Stream
        ObsPy stream object.
    mseed_storage : str
        Output storage template.
    threshold : int, optional
        Minimum waveform duration threshold.
    picker_args : dict, optional
        Picker validation parameters.
    ppc_and_comment : sequence, optional
        Tuple/list containing:

        * Preprocessed flag
        * Log comment
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
    """
    Generate a MiniSEED filename from a template.

    Parameters
    ----------
    template : str
        Output path template.
    trace : Trace
        ObsPy trace object.
    ppc : bool, default=False
        Append preprocessing suffix.

    Returns
    -------
    str
        Generated file path.

    Raises
    ------
    TypeError
        If the generated path is invalid.
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
    """
    Generate chunked time intervals.

    Parameters
    ----------
    starttime : UTCDateTime
        Start time.
    endtime : UTCDateTime
        End time.
    chunklength_in_sec : int, optional
        Chunk duration in seconds.
    overlap_in_sec : int, optional
        Chunk overlap in seconds.

    Returns
    -------
    list
        List of ``(starttime, endtime)`` tuples.

    Raises
    ------
    ValueError
        If ``chunklength_in_sec`` equals zero.
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


def get_st_according2preference(st,location_list,channel_list):
    """
    Suppose that your preference_type is "location" and your 
    location_list is ["00","20","10"], then this function first
    filter the stream according to location and returns
    a  new stream only with location "00", if no exist "00" will 
    continue with the next preference "20", and otherwise, "10". 
    After that, it is going to take new stream and it will go to 
    filter according to channel_list preference if the new stream 
    has more than one channel type ("HH or BH). 
 
    Parameters
    ----------
    st : Stream
        ObsPy stream object.
    location_list : list
        Preferred locations ordered by priority.
    channel_list : list
        Preferred channel prefixes ordered by priority.

    Returns
    -------
    Stream
        Filtered stream.
    """
    logger =logging.getLogger(f'Downloader: preference')
    if len(st) >= 1000:
        stats = st[0].stats
        # logger.error(f"{stats.network}-{stats.station}: "+
        # f"No filter by preference. gaps: {len(st)}  ")
 
        printlog("error",'Downloader: preference',
                f"{stats.network}-{stats.station}"+
                f"No filter by preference. gaps: {len(st)}")
        return st
 
    preference = list(map(lambda x: (x[2],x[3]),st._get_common_channels_info().keys() ))
 
    if not location_list:
        stats = st[0].stats
        new_st = st
    else:
        locations = list(map(lambda x: x[2],st._get_common_channels_info().keys() ))
    
        index = 0
        loc_pref = None
        # print(len(location_list))
        while index < len(location_list):
            loc_pref = location_list[index]
            if loc_pref in locations:
                index = len(location_list)
            else:
                loc_pref = None
                index += 1
 
        if loc_pref == None:
            return st
 
        stats = st[0].stats
        new_st = st.select(network=stats.network, station = stats.station,
                            location=loc_pref)
 
    if not channel_list:
        pass
    else:
        ## If the same location has two differents sensor
        ## example : 00.HH? or 00.BH?, then choose 
        common_channels = new_st._get_common_channels_info().keys()
 
        common_channels = list(common_channels)
        common_channels = list(map(lambda x: x[3][:2],common_channels))
        index = 0
        cha_pref = None
        while index < len(channel_list):
            cha_pref = channel_list[index]
            if cha_pref in common_channels:
                index = len(channel_list)
            else:
                cha_pref = None
                index += 1
 
        if cha_pref == None:
            return st
        else:
            cha_pref = f'{cha_pref}?'
        
 
        new_st = new_st.select(network=stats.network, station = stats.station,
                                location=loc_pref, channel=cha_pref)
 
    common_channels = list(map(lambda x: (x[2],x[3]),new_st._get_common_channels_info().keys()))
 
    printlog("debug",f"Downloader: {stats.network}-{stats.station}: ",
                f"available:{preference},"+
                f" selected {common_channels} according to preference")
    # logger.info(f"{stats.network}-{stats.station}: "+
    #    f"available:{preference},"+
    #    f" selected {common_channels} according to preference")
    return new_st


def get_filenames(mseed, filter_net=[], filter_sta=[], filter_cha=[]):
    """
    Retrieve MiniSEED filenames with optional filtering.

    Parameters
    ----------
    mseed : str
        Directory containing MiniSEED files.

    filter_net : list of str, optional
        Network codes to exclude.

    filter_sta : list of str, optional
        Station codes to exclude.

    filter_cha : list of str, optional
        Channel codes to include.

    Returns
    -------
    list of str
        Filtered MiniSEED filenames.

    Raises
    ------
    TypeError
        If any filter parameter is not a list.

    Examples
    --------
    >>> files = get_filenames("/data/mseed")
    >>> isinstance(files, list)
    True
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
    """
    Download waveforms for multiple parameter combinations.

    Parameters
    ----------
    client : Client
        ObsPy-compatible client instance.

    **kwargs : dict
        Waveform request parameters. Expected keys include:

        * ``network``
        * ``station``
        * ``location``
        * ``channel``
        * ``starttime``
        * ``endtime``

        Comma-separated values are automatically expanded.

    Returns
    -------
    Stream
        Combined waveform stream.

    Notes
    -----
    Each combination of network, station, location, and channel
    is requested independently.
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
    """
    Filter an inventory using waveform selection criteria.

    Parameters
    ----------
    inventory : Inventory
        Input ObsPy inventory.

    network : str
        Comma-separated network codes.

    station : str
        Comma-separated station codes.

    location : str
        Comma-separated location codes.

    channel : str
        Comma-separated channel codes.

    starttime : UTCDateTime
        Selection start time.

    endtime : UTCDateTime
        Selection end time.

    Returns
    -------
    Inventory
        Filtered inventory object.
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
    """
    Retrieve and optionally process waveform data.

    Parameters
    ----------
    client : Client
        ObsPy-compatible waveform client.

    bulk : tuple
        Waveform request tuple containing:

        * Network
        * Station
        * Location
        * Channel
        * Starttime
        * Endtime

    waveform_restrictions : object
        Waveform restriction configuration.

    processing : object, optional
        Processing pipeline object.

    Returns
    -------
    tuple
        Tuple containing:

        * Stream
        * Preprocessed flag
        * Comment string
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
    """
    Download and write waveform data to MiniSEED files.

    Parameters
    ----------
    client : Client
        ObsPy-compatible client.

    bulk : tuple
        Waveform request tuple.

    waveform_restrictions : object
        Waveform restriction configuration.

    download_restrictions : DownloadRestrictions
        Download restriction configuration.

    processing : object, optional
        Processing pipeline.
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
    """
    Merge inventories and station metadata from providers.

    Parameters
    ----------
    providers : list
        List of provider objects.

    Returns
    -------
    tuple
        Tuple containing:

        * Merged inventory
        * Station metadata dictionary
        * Updated providers
        * Stations outside requested domains
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
    """
    Determine whether a point lies inside a polygon.

    Parameters
    ----------
    point : tuple
        Point coordinates ``(longitude, latitude)``.

    polygon_points : list of tuple
        Polygon vertices.

    Returns
    -------
    bool
        ``True`` if the point is inside the polygon.
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
    """
    Filter an inventory and generate station metadata.

    Parameters
    ----------
    inventory : Inventory
        ObsPy inventory object.

    filter_networks : list of str, optional
        Networks to exclude.

    filter_stations : list of str, optional
        Stations to exclude.

    filter_domain : list of float, optional
        Geographic domain defined as:

        ``[minlon, maxlon, minlat, maxlat]``

    Returns
    -------
    tuple
        Tuple containing:

        * Filtered inventory
        * Station metadata dictionary
        * Stations outside domain
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

# if __name__ == "__main__":
#     from obspy.clients.fdsn import Client as FDSN_Client
#     from obspy.core.utcdatetime import UTCDateTime
    # from .restrictions import PreprocRestrictions

    # IRIS_client = FDSN_Client(base_url="IRIS", user='gaprietogo@unal.edu.co',password="DaCgmn3hNjg")
    # st = IRIS_client.get_waveforms(network="YU",station="FC04,FC01",location="",channel="*",starttime=UTCDateTime("2016-04-22T00:00:00.0"),endtime=UTCDateTime("2016-04-23T00:00:00.0"))
    # st = get_st_according2preference(st,["","00","20","10"],["HH","BH"])
    # prove = get_all_sdswaveforms(client=IRIS_client,network="YU",station="FC04,FC01",location="",channel="HH?",starttime=UTCDateTime("2016-04-22T00:00:00.0"),endtime=UTCDateTime("2016-04-23T00:00:00.0"))
    # print(prove)


    # client = FDSN_Client('http://sismo.sgc.gov.co:8080')
    # inv = client.get_stations(network="CM",
    #                     station="BAR2",
    #                     location="*",
    #                     channel="*",
    #                     starttime=UTCDateTime("2019-04-23T00:00:00.0"),
    #                     endtime=UTCDateTime("2019-04-23T00:02:00.0"))
    # st = client.get_waveforms(network="CM",
    #                     station="BAR2",
    #                     location="*",
    #                     channel="*",
    #                     starttime=UTCDateTime("2019-04-23T00:00:00.0"),
    #                     endtime=UTCDateTime("2019-04-23T00:02:00.0"))
    # ppc_restrictions = PreprocRestrictions(["CM.BAR2"],detrend={'type':'simple'})
    # preproc_stream(st,ppc_restrictions)
    ######### inventory
    # json_path = "/home/ecastillo/repositories/AIpicker_modules/onejson.json"
    # client_baseurl = "http://sismo.sgc.gov.co:8080"

    # restrictions = DownloadRestrictions(network="CM",
    #                       station="BAR2",
    #                       location="*",
    #                       channel="*",
    #                       starttime=UTCDateTime("2019-04-23T00:22:34.5"),
    #                       endtime=UTCDateTime("2019-04-25T00:23:39.5"),
    #                       chunklength_in_sec=5000,
    #                       overlap_in_sec=None,
    #                       groupby='{network}.{station}.{channel}')
    # xml = "/home/ecastillo/repositories/AIpicker_modules/CM.xml"

    # makeStationList(json_path,client_baseurl,restrictions,from_xml=xml)

    ######## get stations
    # json_path = "/home/ecastillo/repositories/AIpicker_modules/onejson.json"
    # client_baseurl = "IRIS"
    # restrictions = DownloadRestrictions(network="CI",
    #                   station="BAK,ARV",
    #                   location="*",
    #                   channel="BH*",
    #                   starttime=UTCDateTime("2020-09-01 00:00:00.00"),
    #                   endtime=UTCDateTime("2020-09-02 00:00:00.00"),
    #                   chunklength_in_sec=3600,
    #                   overlap_in_sec=None,
    #                   groupby='{network}.{station}.{channel}')
    # makeStationList(json_path=json_path,client_baseurl="IRIS",restrictions=restrictions)