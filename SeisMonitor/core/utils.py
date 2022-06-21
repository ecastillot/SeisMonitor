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