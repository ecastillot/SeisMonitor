import pandas as pd
import datetime as dt
from obspy.core.inventory.inventory import (Inventory,read_inventory)

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
	# logger = logging.getLogger('json') 
	toprint = []
	stations_outside_domain = []
	for ev in inventory:
		net = ev.code
		if net not in filter_networks:
			for st in ev:
				sta = st.code
				msg = f'{net}-{sta}'
				# logger.debug(msg)

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
			
	# logger.info(str(toprint) + ' ok')

	return inventory,station_list,stations_outside_domain

# inv = read_inventory("/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml")
# inv,sta,sta_out = get_inv_and_json(inv,filter_domain=[-77.282,-71.467,0.785,6.046])
# print(sta_out)

# -77.282,0.785
# -71.467,6.046

##picks
eqt = "/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/test/seismonitor_picks.csv"
df = pd.read_csv(eqt)

# df["p_arrival_time"] = pd.to_datetime(df["p_arrival_time"])
# df = df[(df["p_arrival_time"] >= dt.datetime(2019,12,1)) & \
#         (df["p_arrival_time"] <= dt.datetime(2020,1,1))]

filter_domain=[-75.3105,-73.6362,3.2033,4.4337]
# filter_domain=[-77.955,-72.442,2.666,7.874]
polygon = [(filter_domain[0],filter_domain[2]),
				(filter_domain[0],filter_domain[3]),
				(filter_domain[1],filter_domain[3]),
				(filter_domain[1],filter_domain[2]),
				(filter_domain[0],filter_domain[2])
				]
is_in_domain = lambda x: inside_the_polygon((x.station_lon,x.station_lat),polygon)
df = df[df.apply(is_in_domain,axis=1)]
df = df.reset_index(drop=True)
sta = df["station"].to_list()
sta= list(set(sta))
print(len(sta))
print(sta)
df.to_csv("/home/emmanuel/EDCT/SeisMonitor/asso_mesetas_test/picks_in_mesetas.csv",index=False)
print(df)

# ##eventos
# eqt = "/home/emmanuel/Tesis/auto/aipicker/events/events_1d/csv/CM/CM_eqt_events__20191201__20210101.csv"
# df = pd.read_csv(eqt)
# df["time_event"] = pd.to_datetime(df["time_event"])
# df = df[(df["time_event"] >= dt.datetime(2019,12,1)) & \
#         (df["time_event"] <= dt.datetime(2020,1,1))]
# # filter_domain=[-75.3105,-73.6362,3.2033,4.4337]
# filter_domain=[-77.955,-72.442,2.666,7.874]
# polygon = [(filter_domain[0],filter_domain[2]),
# 				(filter_domain[0],filter_domain[3]),
# 				(filter_domain[1],filter_domain[3]),
# 				(filter_domain[1],filter_domain[2]),
# 				(filter_domain[0],filter_domain[2])
# 				]
# is_in_domain = lambda x: inside_the_polygon((x.longitude,x.latitude),polygon)
# df = df[df.apply(is_in_domain,axis=1)]
# df = df.reset_index(drop=True)
# df = df.drop_duplicates(subset="id")
# print(len(df))
# df.to_csv("/home/emmanuel/Tesis/events_in_mesetas_region_600_600.csv",index=False)
# print(df)