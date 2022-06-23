
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
import json
import pandas as pd
import concurrent.futures
from obspy import UTCDateTime
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.catalog import read_events
from obspy.core.event import Catalog

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
					picker=None,sort='time_event',export = None, 
					from_format="SC3ML",pick_counts=False, inside_polygon=False ):
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


	seiscomp_file_path = os.path.splitext(seiscomp_file)[0]
	catalog = read_events(seiscomp_file,format = from_format)
	event_list = catalog.events 

	events = []
	for event in event_list:
		loc_id = os.path.basename(str(event.resource_id))
		ev_type = event.event_type
		agency = event.creation_info.agency_id

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
			depth = float(depth)/1000 #in km
		if depth_error != None:
			if depth_error.uncertainty != None:
				depth_error = depth_error.uncertainty/1000 #in km
			else:
				depth_error = None
		else:
			depth_error = None
		## Preferred Magnitude
		if with_magnitude:
			pref_magnitude = event.preferred_magnitude()
			# print("aca",pref_magnitude)
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
				
				## conditions to get the filter_id
				if pick.filter_id == None:
					prob = None
					snr = None
					ev_prob = None
					pick.filter_id = ResourceIdentifier(id=None)
				else:
					if "Probability" in str(pick.filter_id):
						prob,snr,ev_prob = str(pick.filter_id).split("+")
						prob = float(str(prob).split("_")[-1])
						snr = str(snr).split("_")[-1]
						ev_prob = str(ev_prob).split("_")[-1]
						if snr != 'None':
							snr = float(snr)
						if ev_prob != 'None':
							ev_prob = float(ev_prob)
						pick.filter_id= ResourceIdentifier(id=None)
					else:
						prob = None
						snr = None
						ev_prob = None

				if pick.creation_info == None:
					_author = None
				else:
					_author = pick.creation_info.author

				picks[pick.resource_id.id] = {"network_code":pick.waveform_id.network_code,
											"station_code":pick.waveform_id.station_code,
											"location_code":pick.waveform_id.location_code,
											"channel_code":pick.waveform_id.channel_code,
											"phase_hint":pick.phase_hint,
											"time":pick.time,
											"author":_author,
											"probability": prob,
											"snr":snr,
											"ev_prob":ev_prob,
											"time_errors":pick.time_errors,
											"filter_id":pick.filter_id,
											"method_id":pick.method_id,
											"polarity":pick.polarity,
											"evaluation_mode":pick.evaluation_mode,
											"evaluation_status":pick.evaluation_status } 
		
		## Preferred picks
		columns = ["agency","id","time_event","latitude","latitude_uncertainty",
					"longitude","longitude_uncertainty","depth","depth_uncertainty",
					"rms","region","method","earth_model","event_type","magnitude",
					"magnitude_type","picker","network","station","location","channel",
					"pick","time_pick","probability","snr","detection_probability"]
		data = []
		fmt = "%Y-%m-%d %H:%M:%S.%f" #2020-01-01 00:02:38
		for i,arrival in enumerate(pref_origin.arrivals):

			pick = picks[arrival.pick_id.id]

			line = [agency,loc_id,time.strftime(fmt),latitude,latitude_error,\
					longitude,longitude_error,depth,depth_error,rms,region,\
					method, earth_model,ev_type,  magnitude, magnitude_type,\
					pick["author"],pick["network_code"],pick["station_code"],\
					pick["location_code"],pick["channel_code"],pick["phase_hint"],\
					pick["time"].strftime(fmt),pick["probability"],pick["snr"],pick["ev_prob"]]
			data.append(line)

		picks_df = pd.DataFrame(data,columns=columns)
		picks_p = picks_df[ picks_df["pick"] == 'P']
		picks_s = picks_df[ picks_df["pick"] == 'S']
		if picker in ("eqt","EQTransformer","eqtransformer",
					"phasenet","PhaseNet"):
			merge_on = columns[:-5]
		else:
			merge_on = columns[:-6]
		picks_df = pd.merge(picks_p,picks_s,how='left',suffixes=("_p", "_s"),on=merge_on)
		events.append(picks_df)

	events_df = pd.concat(events)

	events_df.dropna(subset=['latitude','longitude'],inplace=True)

	events_df.index+=1
	events_df.index.name= 'No'

	if sort != None:   
		events_df = events_df.sort_values(by=sort,ascending=True,ignore_index=True)
	
	if pick_counts:
		events_df = get_pick_counts(events_df)

	if export != None:
		if os.path.isdir(os.path.dirname(export)) == False:
			os.makedirs(os.path.dirname(export))
		events_df.to_csv(export)
		print(f"Events_csv_file: {export}")
	return events_df
