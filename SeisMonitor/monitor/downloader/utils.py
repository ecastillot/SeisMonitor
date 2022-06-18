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
from obspy.core.inventory.inventory import (Inventory,read_inventory)
from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain

class Provider():
	def __init__(self,client,download_restrictions,
				preproc_restrictions=None,xml=None) -> None:
		self.client = client
		self.download_restrictions = download_restrictions
		self.preproc_restrictions = preproc_restrictions
		self.xml = xml

class DownloadRestrictions(Restrictions):
	def __init__(self,network,station,location,channel,
			  starttime,endtime,
			  chunklength_in_sec=None,
			  overlap_in_sec=0,
			  groupby='{network}.{station}.{channel}',
			  threshold= None,
			  location_preferences=[],
			  channel_preferences=[],
			  filter_networks=[], 
			  filter_stations=[],
			  filter_domain=[-180,180,-90,90], #lonw,lone,lats,latn
			#   to_pick=None
			  ):
		"""
		Restrictions to download mseed 
		
		Parameters:
		-----------
		network: str
			Select one or more network codes. 
			Multiple codes are comma-separated (e.g. "IU,TA"). 
			Wildcards are allowed.
		station: str
			Select one or more SEED station codes. 
			Multiple codes are comma-separated (e.g. "ANMO,PFO"). 
			Wildcards are allowed.
		location: str
			Select one or more SEED location identifiers. 
			Multiple identifiers are comma-separated (e.g. "00,01"). 
			Wildcards are allowed.
		channel: str
			Select one or more SEED channel codes. 
			Multiple codes are comma-separated (e.g. "BHZ,HHZ").
		starttime: obspy.UTCDateTime
			Limit results to time series samples on or 
			after the specified start time.
		endtime: obspy.UTCDateTime
			Limit results to time series samples on or 
			before the specified end time.
		chunklength_in_sec: None or int
			The length of one chunk in seconds. 
			If set, the time between starttime and endtime will be divided 
			into segments of chunklength_in_sec seconds.
		overlap_in_sec: None or int
			For more than one chunk, each segment will have overlapping seconds
		groupby: str
			Download group traces together which have the same metadata given by this parameter. 
			The parameter should name the corresponding keys of the stats object, e.g. '{network}.{station}'. 
			This parameter can take the value 'id' which groups the traces by SEED id.
		threshold: int
			limit of length in seconds, length less than threshold will not be downloaded.
		location_preferences: list
			list of location in the order of the preference. If select the
			location of the first element, then the rest of elements will not be 
			downloaded.
		
		"""
		Restrictions.__init__(self,network=network,station=station,
							location=location,channel=channel,
							starttime=starttime,endtime=endtime,
							chunklength_in_sec=chunklength_in_sec)
		
		self._name = "DownloadRestrictions"
		self.overlap_in_sec = overlap_in_sec
		self.groupby = groupby
		self.threshold = threshold
		self.location_preferences = location_preferences
		self.channel_preferences = channel_preferences
		self.filter_networks = filter_networks
		self.filter_stations = filter_stations
		self.filter_domain = filter_domain

class PreprocRestrictions(object):
	def __init__(self,seed_ids,order=['merge','detrend','taper','normalized'],
				decimate=None,detrend=None,applyfilter=None,
				merge=None,normalize=None,remove_response=None,
				resample=None,taper=None):
		"""
		Restrictions to preprocess a stream selected by seed_ids
		
		Parameters:
		-----------
		seed_ids: list
			Contains each seed_id in the next way: network.station"
			ex: ["IU.ANMO","CM.BAR2"]
		order: list of str
			Order to preprocess the stream.
			ex: ['merge','detrend','taper','normalized']
		decimate: dict
			Contains the parameters of the decimate Stream method 
		detrend: dict 
			Contains the parameters of the detrend Stream method 
		filter: dict
			Contains the parameters of the filter Stream method 
		merge: dict
			Contains the parameters of the merge Stream method 
		normalize: dict
			Contains the parameters of the normalize Stream method 
		remove_response: dict
			Contains the parameters of the remove_response Stream method 
		resample: dict
			Contains the parameters of the resample Stream method 
		taper: dict
			Contains the parameters of the taper Stream method 

		--------
		"""
		self.seed_ids = seed_ids
		self.order = order
		self.decimate = decimate
		self.detrend = detrend
		self.applyfilter = applyfilter
		self.merge = merge
		self.normalize = normalize
		self.remove_response = remove_response
		self.resample = resample
		self.taper = taper

def write_stream(one_st,ppc_restrictions,mseed_storage,
				threshold=None, to_pick=None):
	"""
	Write a stream in a specific storage given by mseed_storage

	Parameters:
	-----------
	one_st: obspy.Stream object
		Stream that will be written.
	ppc_restrictions: PreprocRestrictions object
		Restrictions to preprocess a stream.
	mseed_storage:
		Where to store the waveform files.
		The parameter should name the corresponding keys 
		of the stats object,
		e.g. '{network}.{station}.{location}.{channel}__{starttime}__{endtime}
	threshold: int
		data length in seconds
	to_pick: tuple
		(batch_size,overlap)
		It's used to know if the stream can be downloaded according to the 
		picker batch size. If the the segments given by the length of the stream
		and overlap parammeter are less than batcg_size, then no download the stream.

	Returns:
		write one stream
	"""
	one_st,ppc,comment = preproc_stream(one_st,ppc_restrictions)
	one_st = one_st.merge(method=0,fill_value='latest')
	tr = one_st[0]

	mseed_filename = get_mseed_filename(_str=mseed_storage, 
									network=tr.stats.network, 
									station=tr.stats.station,
									location=tr.stats.location, 
									channel=tr.stats.channel,
									starttime=tr.stats.starttime, 
									endtime=tr.stats.endtime,
									ppc=ppc)
	if threshold != None:
		length = abs(tr.stats.endtime - tr.stats.starttime)
		if length < threshold:
			comment = f"length:{length} < threshold:{threshold}"
			# logger =logging.getLogger('Downloaded: False') 
			# logger.info(f'{mseed_filename}  {comment}')
			printlog("info",'Downloader: False',
					f'{mseed_filename}  {comment}')
			download = False
		else:
			download = True

	if to_pick != None:
		batch_size, overlap = to_pick
		length = abs(tr.stats.endtime - tr.stats.starttime)
		mybatch = (length + overlap*length)/60
		if mybatch < batch_size:
			comment = (f"This mseed only can be used with {int(mybatch)} batchs. "+\
					f"In order to download the data, the batch size must be >= {batch_size}."+\
					f" Modify this condition changing 'to_pick':{to_pick}' parameter.")
			# logger =logging.getLogger('Downloaded: False') 
			# logger.info(f'{mseed_filename}  {comment}')
			printlog("info",'Downloader: False',
					f'{mseed_filename}  {comment}')
			download = False
		else:
			download = True
	if (threshold == None) or (to_pick == None):
		download = True


	filename = os.path.basename(mseed_filename)
	if os.path.isfile(mseed_filename) == False:

		if download == False:
			return
		else:
			mseed_dir = os.path.dirname(mseed_filename)
			if os.path.isdir(mseed_dir) == False:
				os.makedirs(mseed_dir)
			else:
				pass
			one_st.write(mseed_filename,format="MSEED")
			# logger =logging.getLogger('Downloaded: True')
			# logger.info(f'{mseed_filename}  {comment}') 

			printlog("info",'Downloader: True',
					f'{mseed_filename}  {comment}')

	else:
		# logger =logging.getLogger('Downloaded: Exist')
		# logger.info(f'{mseed_filename}  {comment}') 
		printlog("info",'Downloader: Exist',
					f'{mseed_filename}  {comment}')

def get_mseed_filename(_str, network, station, location, channel,
					   starttime, endtime,ppc=False):
	"""
	Function taken from obspy

	Helper function getting the filename of a MiniSEED file.

	If it is a string, and it contains ``"{network}"``,  ``"{station}"``,
	``"{location}"``, ``"{channel}"``, ``"{starttime}"``, and ``"{endtime}"``
	formatting specifiers, ``str.format()`` is called.

	Otherwise it is considered to be a folder name and the resulting
	filename will be
	``"FOLDER_NAME/NET.STA.LOC.CHAN__STARTTIME__ENDTIME.mseed"``

	In the last two cases, the times will be formatted with
	``"%Y%m%dT%H%M%SZ"``.
	"""
	strftime = "%Y%m%dT%H%M%SZ"

	if ppc == True:
		ppc_str = ".ppc"
	else:
		ppc_str = ""

	if ("{network}" in _str) and ("{station}" in _str) and \
			("{location}" in _str) and ("{channel}" in _str) and \
			("{starttime}" in _str) and ("{endtime}" in _str) and \
			("{ppc}" in _str):
		
		path = _str.format(
			network=network, station=station, location=location,
			channel=channel, starttime=starttime.strftime(strftime),
			endtime=endtime.strftime(strftime), ppc=ppc_str)
	elif ("{network}" in _str) and ("{station}" in _str) and \
			("{location}" in _str) and ("{channel}" in _str) and \
			("{starttime}" in _str) and ("{endtime}" in _str):
		path = _str.format(
			network=network, station=station, location=location,
			channel=channel, starttime=starttime.strftime(strftime),
			endtime=endtime.strftime(strftime))
	else:
		path = os.path.join(
			_str,
			"{network}.{station}.{location}.{channel}__{s}__{e}.{ppc}.mseed".format(
				network=network, station=station, location=location,
				channel=channel, s=starttime.strftime(strftime),
				e=endtime.strftime(strftime), ppc=ppc_str) )
	if path is True:
		return True
	elif not isinstance(path, (str, bytes)):
		raise TypeError("'%s' is not a filepath." % str(path))
	return path

def preproc_stream(st,ppc_restrictions):
	"""
	Parameters:
	-----------
	st: obspy.Stream object
		Stream object to preprocessing
	ppc_restrictions: PreprocRestrictions object
		Restrictions to preprocess a stream

	Returns:
	--------
	st: obspy.Stream object
		Preprocessed stream according to the order of the parameters. 
	processed: True or False
		True if was processed, False if not.
	"""
	if ppc_restrictions == None:
		processed = False
		comment = ""
	else:
		tr = st[0]
		seed_id = f"{tr.stats.network}.{tr.stats.station}"
		comment = ""
		if seed_id in ppc_restrictions.seed_ids:
			for i,process in enumerate(ppc_restrictions.order):
				try:
					if process == 'decimate':
						st.decimate(**ppc_restrictions.decimate)
					elif process == 'detrend':
						st.detrend(**ppc_restrictions.detrend)
					elif process == 'applyfilter':
						st.filter(**ppc_restrictions.applyfilter)
					elif process == 'merge':
						st.merge(**ppc_restrictions.merge)
					elif process == 'normalize':
						st.normalize(**ppc_restrictions.normalize)
					elif process == 'remove_response':
						st.remove_response(**ppc_restrictions.remove_response)
					elif process == 'resample':
						st.resample(**ppc_restrictions.resample)
					elif process == 'taper':
						st.taper(**ppc_restrictions.taper)
					else:
						str_failed = (f"Failed ppc: {seed_id}-> NO {process}")
						raise Exception( str_failed)

					## only for print comments	
					if i == len(ppc_restrictions.order)-1:
						comment += f"({process}:ok)"
					else:
						comment += f"({process}:ok)->"
				except:
					if i == len(ppc_restrictions.order)-1:
						comment += f"({process}:Failed)"
					else:
						comment += f"({process}:Failed)->"
				processed = True
			comment = f"[{comment}]"
		else:
			processed = False

	return st, processed, comment

def get_chunktimes(starttime,endtime,chunklength_in_sec, overlap_in_sec=0):
	"""
	Make a list that contains the chunktimes according to 
	chunklength_in_sec and overlap_in_sec parameters.

	Parameters:
	-----------
	starttime: obspy.UTCDateTime object
		Start time
	endtime: obspy.UTCDateTime object
		End time
	chunklength_in_sec: None or int
		The length of one chunk in seconds. 
		The time between starttime and endtime will be divided 
		into segments of chunklength_in_sec seconds.
	overlap_in_sec: None or int
		For more than one chunk, each segment will have overlapping seconds

	Returns:
	--------
	times: list
		List of tuples, each tuple has startime and endtime of one chunk.
	"""

	if chunklength_in_sec == 0:
		raise Exception("chunklength_in_sec must be different than 0")
	elif chunklength_in_sec == None:
		return [(starttime,endtime)]

	if overlap_in_sec == None:
		overlap_in_sec = 0

	deltat = starttime
	dtt = dt.timedelta(seconds=chunklength_in_sec)
	overlap_dt = dt.timedelta(seconds=overlap_in_sec)

	times = []
	while deltat < endtime:
		# chunklength can't be greater than (endtime-startime)
		if deltat + dtt > endtime:
			break
		else:
			times.append((deltat,deltat+dtt))
			deltat += dtt - overlap_dt

	if deltat < endtime:	
		times.append((deltat,endtime))
	# print(times)
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

	Parameters:
	-----------
	st: stream object
		stream 
	location_list: list
		locations in order of the preference ["00","20","10"]
	channel_list: list
		channels in order of the preference ["HH","BH"]

	results:
	--------
	new_st : stream object
		stream according to the preference
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
	# 	f"available:{preference},"+
	# 	f" selected {common_channels} according to preference")
	return new_st

def get_filenames(mseed, filter_net=[],
				filter_sta=[],filter_cha=[]):
	"""
	Parameters:
	-----------
	mseed: str
		path of the folder where you need to get all mseed filenames.
	filter_net: list
		avoid mseed if its filename contain one network specified in
		filter net
	filter_sta: list
		avoid mseed if its filename contain one station specified in
		filter net
	filter_cha: list
		avoid mseed if its filename contain one channel specified in
		filter net

	Returns:
	--------
	list
	Get all filenames contained in mseed folder.
	"""

	if not isinstance(filter_net, list):
		raise Exception("filter_net is a list")
	elif not isinstance(filter_sta, list):
		raise Exception("filter_sta is a list")
	elif not isinstance(filter_cha, list):
		raise Exception("filter_cha is a list")

	filenames = []
	file_list = [ev for ev in os.listdir(mseed) \
				if ev.split('/')[-1].split('.')[-1].lower() == 'mseed']
	for onefile in file_list:
		filename = onefile.split('/')[-1]
		info = filename.split('__')[0]
		net,sta,loc,cha = info.split('.')
		if net in filter_net:
			pass
		if sta in filter_sta:
			pass
		if len(filter_cha) != 0:
			if cha not in filter_cha:
				pass
		else:
			filenames.append(filename)
	return filenames

def get_all_sdswaveforms(client,**kwargs):
	"""
	client: Client from obspy
		Client to get_waveforms
	**kwargs: 
		network,station,location,channel

	"""
	myargs = {}
	for key,value in kwargs.items():
		if key in ("network","station","location","channel"):
			myvalue = value.split(",")
		else: 
			myvalue = value
		myargs[str(key)] = myvalue

	st = Stream()
	for net in myargs["network"]:
		for sta in myargs["station"]:
			for loc in myargs["location"]:
				for cha in myargs["channel"]:
					try:
						myst = client.get_waveforms(net,sta,loc,cha,
								myargs['starttime'],myargs['endtime'])

						if len(myst) != 0:
							st += myst

					except Exception as e:
						# logger = logging.getLogger('get waveform')
						# logger.warning(f'not obtained: {net}.{sta}.{loc}.{cha}.
						# 		{myargs['starttime']}.{myargs['endtime']}:{e}') 
						# print(f"not obtained: {net}.{sta}.{loc}.{cha}."+\
						#	   f"{myargs['starttime']}.{myargs['endtime']}:{e}")
						seedname = f"{net}.{sta}.{loc}.{cha}.{myargs['starttime']}.{myargs['endtime']}"
						printlog("warning",seedname,f"{e}")

	return st

def select_inventory(inv,network,station,location,
                    channel,starttime,endtime):
    """
    Filter inventory according to parameters
    Parameters
    ----------
    inv:Inventory
        Inventory that will be filtered
    network: str
        Select one or more network codes. 
        Multiple codes are comma-separated (e.g. "IU,TA"). 
        Wildcards are allowed.
    station: str
        Select one or more SEED station codes. 
        Multiple codes are comma-separated (e.g. "ANMO,PFO"). 
        Wildcards are allowed.
    location: str
        Select one or more SEED location identifiers. 
        Multiple identifiers are comma-separated (e.g. "00,01"). 
        Wildcards are allowed.
    channel: str
        Select one or more SEED channel codes. 
        Multiple codes are comma-separated (e.g. "BHZ,HHZ").
    starttime: obspy.UTCDateTime
        Limit results to time series samples on or 
        after the specified start time.
    endtime: obspy.UTCDateTime
        Limit results to time series samples on or 
        before the specified end time.
    Returns:
    --------
    inventory: Inventory
        Filtered Inventory
    """
    networks = network.split(',')
    stations = station.split(',')
    locations = location.split(',')
    channels = channel.split(',')

    cha0 = channels[0]
    inventory = inv.select(network=networks[0],
                            station=stations[0],
                            location=locations[0],
                            channel=cha0,
                            starttime=starttime,
                            endtime=endtime)
    for net in networks:
        for sta in stations:
            for loc in locations:
                for cha in channels:
                    one_inv = inv.select(network=net,station=sta,
                                            location=loc,channel=cha,
                                            starttime=starttime,
                                            endtime=endtime)
                    inventory = inventory.__add__(one_inv) 
    return inventory

def get_fdsn_waveforms(client,bulk,
						dld_restrictions,
						ppc_restrictions,
						mseed_storage):
	
	net,sta,loc,cha,starttime,endtime = bulk
	loc_preference = dld_restrictions.location_preferences
	cha_preference = dld_restrictions.channel_preferences
	groupby = dld_restrictions.groupby

	try:
		st = client.get_waveforms(net,sta,loc,cha,starttime,endtime)
	except:
		st = Stream()
		return st
	
	if len(st)==0:
		return st

	st= get_st_according2preference(st,
								loc_preference,
								cha_preference)
	st_dict = st._groupby(groupby)
	
	for one_st in st_dict.values():
		write_stream(one_st,ppc_restrictions,
						mseed_storage,
						dld_restrictions.threshold,
						dld_restrictions.to_pick)
	
	return st_dict

def get_merged_inv_and_json(providers):
	json_info = {}
	stations_outside_domains = []
	inventory = Inventory()
	updated_providers = []
	for provider in providers:
		client = provider.client
		restrictions = provider.download_restrictions
		if provider.xml != None:
			inv = read_inventory(provider.xml)
			inv = select_inventory(inv=inv,network=restrictions.network,
										station=restrictions.station,
										location=restrictions.location,
										channel=restrictions.channel,
										starttime=restrictions.starttime, 
										endtime=restrictions.endtime)

		else:
			inv=  client.get_stations(network=restrictions.network,
														station=restrictions.station,
														location=restrictions.location,
														channel=restrictions.channel,
														starttime=restrictions.starttime, 
														endtime=restrictions.endtime, 
														level='channel')

		inv,jinf,sod = get_inv_and_json(inv, 
									restrictions.filter_networks, 
									restrictions.filter_stations,
									restrictions.filter_domain)

		bulk_info = list(set([(net.code,sta.code) for net in inv.networks for sta in net.stations ]))

		provider.download_restrictions.bulk_info = bulk_info

		inventory = inventory.__add__(inv) 
		json_info.update(jinf)
		stations_outside_domains.append(sod)
		updated_providers.append(provider)

	stations_outside_domains = [x for xs in stations_outside_domains for x in xs]
	return inventory,json_info,updated_providers,stations_outside_domains

def inside_the_polygon(p,pol_points):
    """
    Parameters:
    -----------
    p: tuple
        Point of the event. (lon,lat)
    pol_points: list of tuples
        Each tuple indicates one polygon point (lon,lat).
    Returns: 
    --------
    True inside 
    """
    V = pol_points

    cn = 0  
    V = tuple(V[:])+(V[0],)
    for i in range(len(V)-1): 
        if ((V[i][1] <= p[1] and V[i+1][1] > p[1])   
            or (V[i][1] > p[1] and V[i+1][1] <= p[1])): 
            vt = (p[1] - V[i][1]) / float(V[i+1][1] - V[i][1])
            if p[0] < V[i][0] + vt * (V[i+1][0] - V[i][0]): 
                cn += 1  
    condition= cn % 2  
    
    if condition== 1:   
        return True
    else:
        return False

def get_inv_and_json(inventory, 
			filter_networks=[], 
			filter_stations=[],
			filter_domain=[-180,180,-90,90]):


	"""
	
	Uses fdsn to find available stations in a specific geographical location and time period.  
	Parameters
	----------
	json_path: str
		Path of the json file that will be returned
	providers:list
		list of Client object
	restrictions: DownloadRestrictions or MassDownloaderRestrictions
		Restrictions to download mseed
	from_xml : str
		Path of xml file
	**kwargs: str
		get_stations kwargs


	Returns
	----------
	stations_list.json: A dictionary containing information for the available stations.	  
		
	"""  
	if not filter_domain:
		filter_domain=[-180,180,-90,90]
		
	polygon = [(filter_domain[0],filter_domain[2]),
				(filter_domain[0],filter_domain[3]),
				(filter_domain[1],filter_domain[3]),
				(filter_domain[1],filter_domain[2]),
				(filter_domain[0],filter_domain[2])
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

					elv = st.elevation
					lat = st.latitude
					lon = st.longitude

					is_in_domain = inside_the_polygon((lon,lat),polygon)

					if not is_in_domain:
						inventory = inventory.remove(network=net,station=sta)
						stations_outside_domain.append(sta)
						continue

					new_chan = [ch.code for ch in st.channels]

					if len(new_chan) > 0 and (sta not in station_list):
						channels = list(set(new_chan))
						sample_rates_gen = lambda x: st.select(channel=x).channels[0].sample_rate
						sample_rates = list(map(sample_rates_gen,channels))
						station_list[str(sta)] ={"network": net,
												"channels": list(set(new_chan)),
												"coords": [lat, lon, elv],
												"sampling_rate": sample_rates
												}

					toprint.append(msg)
				else:
					inventory = inventory.remove(network=net,station=sta)
		else:
			inventory = inventory.remove(network=net)
			
	logger.info(str(toprint) + ' ok')

	return inventory,station_list,stations_outside_domain

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