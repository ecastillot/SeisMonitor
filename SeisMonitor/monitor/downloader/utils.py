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

class DownloadRestrictions():
	def __init__(self,mseed_storage, 
                chunklength_in_sec=None,
                threshold= 60,
                overlap_in_sec=0,
                picker_args={}, # {"batch_size":10,"overlap":0.3,"lenght":60}
                groupby='{network}.{station}.{channel}',
                n_processor=None) -> None:
		self.mseed_storage = mseed_storage
		self.chunklength_in_sec = chunklength_in_sec
		self.threshold = threshold
		self.overlap_in_sec = overlap_in_sec
		self.picker_args = picker_args
		self.groupby = groupby
		self.n_processor=n_processor
		

def sanitize_provider_times(providers):
    provider_times = [(provider.waveform_restrictions.starttime,
                        provider.waveform_restrictions.endtime)\
                        for provider in providers]
    if provider_times.count(provider_times[0]) == len(provider_times):
        return  providers
    else:
        raise Exception("Providers must have the same interval time")

def get_max_allowed_batch_size(data_length,segment_length,overlap):
    """
    data_length : length of the data in seconds
    segment_length: length of each batch segment in seconds 
    overlap: overlap (0-1)
    """
    max_batch_size = (data_length-(overlap*segment_length)) / (segment_length*(1-overlap))
    max_batch_size = int(max_batch_size)

    if max_batch_size == 0:
        max_batch_size = 1

    return max_batch_size

def write_stream(st,mseed_storage, 
				threshold=None, 
				picker_args={}, # {"batch_size":10,"overlap":0.3,"lenght":60}
				ppc_and_comment = [False,""]):
	"""
	Write a stream in a specific storage given by mseed_storage

	Parameters:
	-----------
	st: obspy.Stream object
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
	picker_args : dict
		keys: batch_size,overlap,length
		It's used to know if the stream can be downloaded according to the 
		picker keys. If the the segments given by the length of the stream
		and overlap parammeter are less than batch_size, then no download the stream.

	Returns:
		write one stream
	"""
	

	ppc,comment = ppc_and_comment
	tr = st[0]

	mseed_filename = get_mseed_filename(_str=mseed_storage, 
									tr=tr,
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

	if picker_args:
		overlap,batch_size,segment_length = picker_args["overlap"],picker_args["batch_size"],picker_args["length"]
		data_length = abs(tr.stats.endtime - tr.stats.starttime)
		max_batch_size = get_max_allowed_batch_size(data_length,segment_length,overlap)
		
		if max_batch_size < batch_size:
			comment = (f"This mseed only can be used with {max_batch_size} batchs. "+\
					f"In order to download the data, the batch size must be >= {batch_size}."+\
					f" Modify this condition changing 'picker_args':{picker_args}' parameter.")
			# logger =logging.getLogger('Downloaded: False') 
			# logger.info(f'{mseed_filename}  {comment}')
			printlog("info",'Downloader: False',
					f'{mseed_filename}  {comment}')
			download = False
		else:
			download = True
	if (threshold == None) or (bool(picker_args) == False):
		download = True

	if os.path.isfile(mseed_filename) == False:

		if download == False:
			return
		else:
			mseed_dir = os.path.dirname(mseed_filename)
			if os.path.isdir(mseed_dir) == False:
				os.makedirs(mseed_dir)
			else:
				pass
			st.write(mseed_filename,format="MSEED")
			# logger =logging.getLogger('Downloaded: True')
			# logger.info(f'{mseed_filename}  {comment}') 

			printlog("info",'Downloader: True',
					f'{mseed_filename}  {comment}')

	else:
		# logger =logging.getLogger('Downloaded: Exist')
		# logger.info(f'{mseed_filename}  {comment}') 
		printlog("info",'Downloader: Exist',
					f'{mseed_filename}  {comment}')
	del tr; del st


def get_mseed_filename(_str, tr,ppc=False):
	"""
	It generates the path of the file.

	_str : str
		string with the wildcards
		possible wildcards: network,station,location,channel,starttime,endtime,year,month,day,julday
	tr: Trace
		trace to take the information
	ppc: bool
		True if it was preprocesed. If true .pcc is added in the path.
	"""
	strftime = "%Y%m%dT%H%M%SZ"
	network=tr.stats.network 
	station=tr.stats.station
	location=tr.stats.location 
	channel=tr.stats.channel
	starttime=tr.stats.starttime 
	_starttime=starttime.strftime(strftime) 
	endtime=tr.stats.endtime
	_endtime=endtime.strftime(strftime)
	year = starttime.year
	month = starttime.month
	day = endtime.day
	julday = endtime.julday
	path = _str.format(
                    network=network, station=station, location=location,
                    channel=channel, year=year, month=month, 
                    day=day, julday=julday,starttime=_starttime,endtime=_endtime)
	if ppc:
		path = path + ".ppc"

	if not isinstance(path, (str, bytes)):
		raise TypeError("'%s' is not a filepath." % str(path))
	return path

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
						printlog("error",seedname,f"{e}")

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

def get_client_waveforms(client,bulk,
						waveform_restrictions,
						processing):
	strftime = "%Y%m%dT%H%M%SZ"
	net,sta,loc,cha,starttime,endtime = bulk
	loc_preference = waveform_restrictions.location_preferences
	cha_preference = waveform_restrictions.channel_preferences

	why = "-".join((net,sta,loc,cha,starttime.strftime(strftime),endtime.strftime(strftime)))
	try:
		st = client.get_waveforms(net,sta,loc,cha,starttime,endtime)
	except Exception as e:
		printlog("info","Downloader: False",why+"->"+str(e))
		st = Stream()
		ppc = False
		comment = ""
		return st,ppc,comment
	# print(st)
	if len(st)==0:
		printlog("warning","Downloader: False","0 Trace(s) in Stream: "+why)
		ppc = False
		comment = ""
		return st,ppc,comment

	st= get_st_according2preference(st,
								loc_preference,
								cha_preference)	
	if processing == None:
		ppc = False
		comment = ""
	else:
		st,ppc,comment = processing.run(st)
	return st,ppc,comment

def write_client_waveforms(client,bulk,
						waveform_restrictions,
						download_restrictions,
						processing):
	st,ppc,comment = get_client_waveforms(client,bulk,
						waveform_restrictions,
						processing)
	groupby = download_restrictions.groupby
	st_dict = st._groupby(groupby)
	for st in st_dict.values():
		write_stream(st,
					download_restrictions.mseed_storage,
					download_restrictions.threshold,
					download_restrictions.picker_args,
					[ppc,comment])
	# del st;st_dict;del client;del bulk;del waveform_restrictions;del processing
	# return st_dict

def get_merged_inv_and_json(providers):
	json_info = {}
	stations_outside_domains = []
	inventory = Inventory()
	updated_providers = []
	for provider in providers:
		client = provider.client
		restrictions = provider.waveform_restrictions
		if provider.xml != None:
			# print("with xml")
			# print(restrictions.__dict__)
			try:
				inv = read_inventory(provider.xml)
				inv = select_inventory(inv=inv,network=restrictions.network,
											station=restrictions.station,
											location=restrictions.location,
											channel=restrictions.channel,
											starttime=restrictions.starttime, 
											endtime=restrictions.endtime)
			except:
				printlog("error","Inventory",f"No get_stations with {restrictions.__dict__}")
				inv = Inventory()
		else:
			# print("no xml")
			# print(restrictions.__dict__)
			try:
				inv=  client.get_stations(network=restrictions.network,
															station=restrictions.station,
															location=restrictions.location,
															channel=restrictions.channel,
															starttime=restrictions.starttime, 
															endtime=restrictions.endtime, 
															level='channel')
			except:
				printlog("error","Inventory",f"No get_stations with {restrictions.__dict__}")
				inv = Inventory()
				

		inv,jinf,sod = get_inv_and_json(inv, 
									restrictions.filter_networks, 
									restrictions.filter_stations,
									restrictions.filter_domain)

		bulk_info = list(set([(net.code,sta.code) for net in inv.networks for sta in net.stations ]))

		provider.waveform_restrictions.bulk_info = bulk_info

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