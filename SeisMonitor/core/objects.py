from . import utils as ut
from obspy import read_inventory
import copy

class Provider():
	def __init__(self,client,waveform_restrictions, 
				processing=None,xml=None) -> None:
		self.client = client
		self.waveform_restrictions = waveform_restrictions
		self.processing = processing
		self.xml = xml

		if xml != None:
			self.inventory = read_inventory(xml)
		else:
			self.inventory = client.get_stations(network=waveform_restrictions.network,
												station=waveform_restrictions.station,
												location="*",
												channel="*",level='response')


	def copy(self):
		return copy.deepcopy(self)

class WaveformRestrictions():
	def __init__(self,network,station,
			  location,channel,
			  starttime,endtime,
			  location_preferences=[],
			  channel_preferences=[],
			  filter_networks=[], 
			  filter_stations=[],
			  filter_domain=[-180,180,-90,90] #lonw,lone,lats,latn
			  ):
		"""
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
		location_preferences: list
			list of location in the order of the preference. If select the
			location of the first element, then the rest of elements will not be 
			downloaded.
		"""
		self.network = network
		self.station = station
		self.location = location
		self.channel = channel
		self.starttime = starttime
		self.endtime = endtime
		self.location_preferences = location_preferences
		self.channel_preferences = channel_preferences
		self.filter_networks = filter_networks
		self.filter_stations = filter_stations
		self.filter_domain = filter_domain

class Processing(object):
	def __init__(self,
				order=['normalize','merge','detrend',
						'taper',"filter"],
				decimate={"factor":2},
				detrend={"type":"demean"},
				filter={"type":'bandpass', 
						"freqmin" : 1, 
						"freqmax" : 45, 
						"corners":2, 
						"zerophase":True},
				merge={"method":0,"fill_value":'latest'},
				normalize={"global_max":False},
				resample={"sampling_rate":200},
				taper={"max_percentage":0.001, 
						"type":"cosine", 
						"max_length":2} ,
				select_networks=[], 
				select_stations=[],
				filter_networks=[], 
			  	filter_stations=[]):
		"""
		Restrictions to preprocess a stream selected by seed_ids
		
		Parameters:
		-----------
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
		self.order = order
		self.decimate = decimate
		self.detrend = detrend
		self.filter = filter
		self.merge = merge
		self.normalize = normalize
		self.resample = resample
		self.taper = taper
		self.select_networks = select_networks
		self.select_stations = select_stations
		self.filter_networks = filter_networks
		self.filter_stations = filter_stations
	
	def run(self,st):
		return ut.preproc_stream(st,self.order ,
							self.decimate ,
							self.detrend ,
							self.filter ,
							self.merge ,
							self.normalize ,
							self.resample ,
							self.taper ,
							self.select_networks ,
							self.select_stations ,
							self.filter_networks ,
							self.filter_stations)
