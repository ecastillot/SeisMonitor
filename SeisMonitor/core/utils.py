
"""
 * @author Emmanuel David Castillo Taborda
 * @email ecastillot@unal.edu.co
 * @create date 2021-03-03 23:51:21
 * @modify date 2021-06-02 04:31:09
 * @desc [description]
"""

# import sys
# AIpicker_path = "/home/ecastillo/AIpicker"
# sys.path.insert(0,AIpicker_path)

import datetime as dt
import time
import os
import ast
import json
import pandas as pd
import concurrent.futures
from obspy import UTCDateTime
from obspy.core.event.base import CreationInfo
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.catalog import read_events
from obspy.core.event import Catalog


def add_aditional_origin_info(origin,
				agency = None,
				region = None,
				method_id = None,
				earth_model_id = None,
				evaluation_mode = "automatic",
				evaluation_status = "preliminary"):

	if region != None:
		origin.region = region
	if earth_model_id != None:
		origin.earth_model_id = ResourceIdentifier(id=earth_model_id)
	if method_id != None:
		origin.method_id = ResourceIdentifier(id=method_id)
		
	origin.evaluation_mode = evaluation_mode
	origin.evaluation_status = evaluation_status
	origin.creation_info = CreationInfo(agency_id=agency,
                                                    agency_uri=ResourceIdentifier(id=agency),
                                                    author="SeisMonitor",
                                                    author_uri=ResourceIdentifier(id="SeisMonitor"),
                                                    creation_time=UTCDateTime.now())
	return origin

def add_aditional_event_info(event,
				agency = None,
				event_type = "earthquake",
				event_type_certainty = "suspected"):
	
	ori_pref  = event.preferred_origin()
	if ori_pref == None:
		event.preferred_origin_id = event.origins[0].resource_id.id

	event.event_type = event_type
	event.event_type_certainty = event_type_certainty
	event.creation_info = CreationInfo(agency_id=agency,
										agency_uri=ResourceIdentifier(id=agency),
										author="SeisMonitor",
										author_uri=ResourceIdentifier(id="SeisMonitor"),
										creation_time=UTCDateTime.now())
	return event

def add_aditional_catalog_info(catalog,
				agency = None):
	catalog.creation_info = CreationInfo(agency_id=agency,
										agency_uri=ResourceIdentifier(id=agency),
										author="SeisMonitor",
										author_uri=ResourceIdentifier(id="SeisMonitor"),
										creation_time=UTCDateTime.now())
	return catalog


def preproc_stream(st,
				order=['normalize','merge','detrend',
						'taper',"filter"],
				decimate=None,detrend=None,filter=None,
				merge=None,normalize=None,
				resample=None,taper=None,
				select_networks=[], 
				select_stations=[],
				filter_networks=[], 
			  	filter_stations=[]):
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
	tr = st[0]
	network = tr.stats.network
	station = tr.stats.station
	comment = ""

	if len(order) == 0:
		processed = False
		return st, processed, comment

	if (network in filter_networks) or\
		(station in filter_stations):
		processed = False

	else:
		if not select_networks:
			select_networks.append(network)
		if not select_stations:
			select_stations.append(station)
		
		if (network in select_networks) or\
		(station in select_stations):
			for i,process in enumerate(order):
				try:
					if process == 'decimate':
						st.decimate(**decimate)
					elif process == 'detrend':
						st.detrend(**detrend)
					elif process == 'filter':
						st.filter(**filter)
						# print("filter applied")
					elif process == 'merge':
						st.merge(**merge)
					elif process == 'normalize':
						st.normalize(**normalize)
					elif process == 'resample':
						st.resample(**resample)
					elif process == 'taper':
						st.taper(**taper)

					## only for print comments	
					if i == len(order)-1:
						comment += f"({process}:ok)"
					else:
						comment += f"({process}:ok)->"
				except:
					if i == len(order)-1:
						comment += f"({process}:Failed)"
					else:
						comment += f"({process}:Failed)->"
				processed = True
			comment = f"[{comment}]"
		else:
			processed = False

	return st, processed, comment

def get_csv_events(seiscomp_file, version="0.9", with_magnitude=True, 
					export = None, 
					from_format="SC3ML", inside_polygon=False ):
	"""
	parameters
	----------
	seiscom_file : str
		path of the sc3ml seiscomp file that contains only one event
	version: str
		String of the xml version number.
	with_magnitude: Bolean
		If True return the magnitude and magnitude information
	picker: None
		In 'eqt' or 'phasenet' doesn't matter the channel where is located the pick.
		While in others matters the channel. None is select to  have importance in the channel
	export: str (deault : None)
		Path to export in a csv file. None don't create any file, then only returns.
	returns
	-------
	appended_events : list
		It's a list of Pandas DataFrames where each dataframe contains all information about the picks from the respective event
	"""
	def change_xml_version(ev_file,new_version="0.9"):
		lines = open(ev_file).readlines()
		lines[1] = f'<seiscomp xmlns="http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/{new_version}" version="{new_version}">\n'
		with open(ev_file, 'w') as f:
			f.write(''.join(lines))

	if from_format == "SC3ML":
		if version not in ["0.5", "0.6", "0.7", "0.8", "0.9", "0.10"]:
			change_xml_version(seiscomp_file,new_version="0.10")

	datefmt = "%Y-%m-%d %H:%M:%S.%f"
	seiscomp_file_path = os.path.splitext(seiscomp_file)[0]
	catalog = read_events(seiscomp_file,format = from_format)
	event_list = catalog.events 


	event_colname= ["n_event","event_id","event_time","latitude","latitude_uncertainty",
					"longitude","longitude_uncertainty","depth","depth_uncertainty",
					"rms","region","method","earth_model","event_type","magnitude",
					"magnitude_type","n_P_phases","n_S_phases"]
	pick_colname = ["n_event","event_id","pick_id","phasehint","arrival_time",
					"probability","snr","detection_probability",
					"network","station","location","channel","picker"
					]
	events_df = []
	picks_df = []
	for n_ev,event in enumerate(event_list):
		loc_id = os.path.basename(str(event.resource_id))
		ev_type = event.event_type
		# if event.creation_info != None:
		# 	agency = event.creation_info.agency_id
		# else:
		# 	agency = None

		if event.event_descriptions:
			region = event.event_descriptions[0].text
		else:
			region = None

		## Preferred Origin
		pref_origin = event.preferred_origin() 
		time = pref_origin.time
		latitude = pref_origin.latitude
		latitude_error = pref_origin.latitude_errors.uncertainty
		longitude = pref_origin.longitude
		longitude_error = pref_origin.longitude_errors.uncertainty

		depth = pref_origin.depth
		depth_error = pref_origin.depth_errors
		rms = pref_origin.quality.standard_error
		method = os.path.basename(str(pref_origin.method_id))
		earth_model = os.path.basename(str(pref_origin.earth_model_id))
		evaluation_mode = pref_origin.evaluation_mode

		if depth != None:
			depth = float(depth) #in km
		if depth_error != None:
			if depth_error.uncertainty != None:
				depth_error = depth_error.uncertainty #in km
			else:
				depth_error = None
		else:
			depth_error = None
		## Preferred Magnitude
		if with_magnitude:
			pref_magnitude = event.preferred_magnitude()
			if pref_magnitude != None:
				magnitude = pref_magnitude.mag
				magnitude_type = pref_magnitude.magnitude_type
				if magnitude != None:
					magnitude = round(magnitude,2)
			else:
				magnitude = None
				magnitude_type = None

		else:
			magnitude = None
			magnitude_type = None
		
		
		## Dictionary with the picks information
		picks = {}
		for pick in event.picks:
			if pick.resource_id.id not in picks.keys():
				
				if pick.creation_info == None:
					_author = None
				else:
					_author = pick.creation_info.author
				
				if pick.comments:
					comment = ast.literal_eval(pick.comments[0].text)
					prob = comment["probability"]
					if _author == "EQTransformer":
						snr = comment["snr"]
						ev_prob = comment["detection_probability"]
					else:
						snr = None
						ev_prob = None
				else:
					prob = None
					snr = None
					ev_prob = None


				

				picks[pick.resource_id.id] = {
											"id":os.path.basename(str(pick.resource_id)),
											"network_code":pick.waveform_id.network_code,
											"station_code":pick.waveform_id.station_code,
											"location_code":pick.waveform_id.location_code,
											"channel_code":pick.waveform_id.channel_code,
											"phase_hint":pick.phase_hint,
											"time":pick.time,
											"author":_author,
											"probability": prob,
											"snr":snr,
											"detection_probability":ev_prob,
											"time_errors":pick.time_errors,
											"filter_id":pick.filter_id,
											"method_id":pick.method_id,
											"polarity":pick.polarity,
											"evaluation_mode":pick.evaluation_mode,
											"evaluation_status":pick.evaluation_status } 
		
		p_count = 0
		s_count = 0
		for i,arrival in enumerate(pref_origin.arrivals):

			pick = picks[arrival.pick_id.id]
			pick_row = [n_ev,loc_id,pick["id"],pick["phase_hint"],
						pick["time"].datetime,pick["probability"],
						pick["snr"],pick["detection_probability"],
						pick["network_code"],pick["station_code"],
						pick["location_code"],pick["channel_code"],
						pick["author"]
						]
			
			picks_df.append(pick_row)

			if pick["phase_hint"].upper() == "P":
				p_count += 1
			elif pick["phase_hint"].upper() == "S":
				s_count += 1

		
		event_row = [n_ev,loc_id,time.datetime,latitude,latitude_error,\
					longitude,longitude_error,depth,depth_error,rms,region,\
					method, earth_model,ev_type,  magnitude, magnitude_type,
					p_count,s_count]

		events_df.append(event_row)

	events_df = pd.DataFrame(events_df,columns=event_colname)
	picks_df = pd.DataFrame(picks_df,columns=pick_colname)


	events_df = events_df.sort_values(by="event_time",
										ascending=True,
										ignore_index=True)
	picks_df = picks_df.sort_values(by="arrival_time",
										ascending=True,
										ignore_index=True)

	if export != None:
		if os.path.isdir(os.path.dirname(export)) == False:
			os.makedirs(os.path.dirname(export))
		events_df.to_csv(export,index=False)
		print(f"Events_csv_file: {export}")
	return events_df,picks_df
