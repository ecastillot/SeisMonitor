# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-10-12 10:58:25
#  * @modify date 2021-10-12 10:58:25
#  * @desc [description]
#  */
######## 
# Moment magnitude adapted from 
# https://github.com/krischer/moment_magnitude_calculator/tree/master/scripts
# :copyright:
# Lion Krischer (krischer@geophysik.uni-muenchen.de), 2012
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.base import TimeWindow,WaveformStreamID
from obspy.core.inventory.inventory import Inventory
from obspy.core.stream import Stream
from obspy.io.xseed.parser import Parser
from obspy.core.event import Magnitude as Mag
from obspy.core.event import Catalog
from obspy.core.event import read_events, Comment
from SeisMonitor.utils import isfile
from SeisMonitor.core import utils as scut
from obspy.core.event.base import CreationInfo
import concurrent.futures as cf
from obspy import UTCDateTime
# import mtspec
import numpy as np
import scipy
import os
from SeisMonitor.utils import printlog
from . import utils as ut

class MwPhysicalMagParams():
    def __init__(self,
                vp=4800,
                vsp_factor=1.73,
                density=2650.0,
                waterlevel=10,
                p_radiation_pattern=0.52,
                s_radiation_pattern=0.63):
        self.vp = vp
        self.vsp_factor = vsp_factor
        self.vs = vp/vsp_factor
        self.density = density
        self.waterlevel = waterlevel
        self.p_radiation_pattern = p_radiation_pattern
        self.s_radiation_pattern = s_radiation_pattern

class MwProcessingMagParams():
    def __init__(self,
                time_before_pick = 5,
                time_after_pick = 30,
                padding = 40,
                only_proc_p_pick=False,
                only_proc_s_pick=False):
        self.time_before_pick = time_before_pick
        self.time_after_pick = time_after_pick
        self.padding = padding
        self.only_proc_p_pick = only_proc_p_pick
        self.only_proc_s_pick = only_proc_s_pick

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
    if (st == None) or (len(st)==0) :
        return None

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

    printlog("debug",f"Magnitude: {stats.network}-{stats.station}: ",
                f"available:{preference},"+
                f" selected {common_channels} according to preference")
    # logger.info(f"{stats.network}-{stats.station}: "+
    #     f"available:{preference},"+
    #     f" selected {common_channels} according to preference")
    return new_st

class Magnitude():
    def __init__(self,providers,catalog,out_dir) -> None:
        # self.client = client
        # super().__init__(sds_archive,sds_type,format)
        # self.response = response
        self.providers = providers
        print("Reading catalog... ")
        if isinstance(catalog,Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)

        self.agency = self.catalog.creation_info.agency_id
        self.out_dir = out_dir
        self.xml_ml_out_file = os.path.join(out_dir,"Ml_magnitude.xml")    
        self.xml_mw_out_file = os.path.join(out_dir,"Mw_magnitude.xml")    

    def get_Ml(self,
                mag_type="RSNC",
                trimmedtime=50,
                padding=20,waterlevel=10,
                zone=None,
                out_format="QUAKEML"):
        events_mag = []
        for n_ev,event in enumerate(self.catalog,1):
            print(f"Event No. {n_ev}:",event.resource_id,f"from {len(self.catalog)} events")

            if not event.origins:
                print ("No origin for event %s" % event.resource_id)
                continue

            ori_pref  = event.preferred_origin()
            if ori_pref == None:
                event.preferred_origin_id = event.origins[0].resource_id.id
                ori_pref  = event.preferred_origin()

            origin_time = ori_pref.time
            latitude = ori_pref.latitude
            longitude = ori_pref.longitude
            depth = ori_pref.depth
            # depth = ori_pref.depth *1e3

            if latitude ==None or longitude==None:
                continue
            
            amplitudes = []
            station_magnitudes = []
            Mls = []

            def _get_maginfo_by_station(pick):
            # for pick in event.picks:
                
                if pick.phase_hint.upper() != "S":
                    # continue
                    return None
                
                st = self._get_corresponding_st(pick.waveform_id, pick.time,
                                                padding)
                if st == None:
                    # continue
                    return None

                stats = st[0].stats

                ev_params = {"picktime":pick.time, "latitude":latitude,
                            "longitude":longitude}

                resp = Inventory()
                for provider in self.providers:
                    if len(resp) > 0:
                        continue
                    response = provider.inventory
                    # selected_response = response.select(network = pick.waveform_id.network_code,
                    #                                     station = pick.waveform_id.station_code,
                    #                                     location = "*",
                    #                                     channel = pick.waveform_id.channel_code[0:2]+"*")
                    selected_response = response.select(network = stats.network,
                                                        station = stats.station,
                                                        location = stats.location,
                                                        channel = stats.channel[:2]+"*" )
                                    
                    resp = selected_response.__add__(selected_response)

                ampl,epi_dist,tr_id = ut.get_Ml_magparams_by_station(st,resp,
                                ev_params,
                                trimmedtime,
                                waterlevel)
                net,sta,loc,cha = tr_id.split(".")
                amp = ut.write_amplitude_values(ampl,
                                                time_window=TimeWindow(begin=0,
                                                                        end=trimmedtime,
                                                                        reference=ev_params["picktime"]),
                                                pick_id=pick.resource_id,
                                                waveform_id=WaveformStreamID(net,sta,loc,cha),
                                                magnitude_hint="Ml",
                                                agency=self.agency)
                
                amplitudes.append(amp)

                if (ampl == None) or (ampl==0):
                    # continue
                    return None

                
                Ml = ut.get_Ml(ampl,epi_dist,mag_type,zone)
                Mls.append(Ml)

                stamag = ut.write_magsta_values(value=Ml,uncertainty=None,
                                                mag_type="Ml",
                                                origin_id =ori_pref.resource_id,
                                                amplitude_id=amp.resource_id,
                                                method_id=ResourceIdentifier(id="https://docs.obspy.org/tutorial/advanced_exercise/advanced_exercise_solution_3b.html"),
                                                waveform_id=WaveformStreamID(net,sta,loc,cha),
                                                agency=self.agency
                                                 )
                station_magnitudes.append(stamag)

                staname = ".".join((stats.network, stats.station,
                              stats.location ,stats.channel[:2]+"*"))
                print(f"\t-> Ml | {staname}-{pick.phase_hint.upper()} | {Ml}")

            with cf.ThreadPoolExecutor() as executor:
                executor.map(_get_maginfo_by_station,event.picks)
                
            if not Mls:
                Ml = 0
                Ml_std = 0
                print(f"No magnitude for {event.resource_id}")
            else:
                Ml = scipy.stats.trim_mean(Mls,0.25)
                Ml_std = np.array(Mls).std()

            mag = ut.write_magnitude_values(Ml,Ml_std,len(Mls),"Ml",
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           method_id=ResourceIdentifier(id="https://docs.obspy.org/tutorial/advanced_exercise/advanced_exercise_solution_3b.html"),
                           origin_id=ori_pref.resource_id,
                           agency=self.agency,
                           comments=None)

            event.amplitudes = amplitudes
            event.station_magnitudes = station_magnitudes
            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

            events_mag.append(event)
            print(f"{event.resource_id} | Ml: {round(Ml,2)} | Ml_std {round(Ml_std,2)} | stations:{len(Mls)}\n\n")

        # _get_maginfo_by_station(pick)
        catalog = Catalog(events = events_mag )
        catalog = scut.add_aditional_catalog_info(catalog,agency=self.agency)

        if self.xml_ml_out_file != None:
            print (f"Writing output file in {self.xml_ml_out_file}")
            isfile(self.xml_ml_out_file)
            catalog.write(self.xml_ml_out_file, out_format)
        return self.catalog

    def get_Mw(self,physparams,procparams,
            out_format="QUAKEML"):

        Mws = []
        Mws_std = []

        for n_ev,event in enumerate(self.catalog,1):
            print(f"Event No. {n_ev}:",event.resource_id,f"from {len(self.catalog)} events.")

            if not event.origins:
                print ("No origin for event %s" % event.resource_id)
                continue

            origin_time = event.origins[0].time
            latitude = event.origins[0].latitude
            longitude = event.origins[0].longitude
            depth = event.origins[0].depth

            moments = []
            corner_frequencies = []
            for pick in event.picks:

                if procparams.only_proc_s_pick:
                    if pick.phase_hint.upper() == "P":
                        continue

                elif procparams.only_proc_p_pick:
                    if pick.phase_hint.upper() == "S":
                        continue
                    

                st = self._get_corresponding_st(pick.waveform_id, pick.time,
                                                procparams.padding )
                
                if st == None:
                    continue

                resp = Inventory()
                for provider in self.providers:
                    if len(resp) > 0:
                        continue

                    if pick.waveform_id.channel_code == None:
                        channel = "*"
                    else:
                        channel = pick.waveform_id.channel_code[0:2] +"*"

                    response = provider.inventory
                    selected_response = response.select(network = pick.waveform_id.network_code,
                                                        station = pick.waveform_id.station_code,
                                                        location = "*",
                                                        channel = channel )
                                    
                    resp = selected_response.__add__(selected_response)
                st = ut.Mw_st_processing(st, resp,
                                            physparams.waterlevel,
                                             pick.time-procparams.padding)

                if st == None:
                    continue

                traveltime = pick.time - origin_time
                M_0,corner_freqs = ut.get_M0_magnitude_by_pick(st,pick.time,traveltime,
                                                    pick.phase_hint,
                                                    physparams,
                                                    procparams)

                if M_0 != None:
                    staname = ".".join((pick.waveform_id.network_code, pick.waveform_id.station_code,
                              pick.waveform_id.location_code ,channel))
                    
                    Mw_sta = 2.0 / 3.0 * (np.log10(M_0) - 9.1)
                    print(f"\t-> Mw | {staname}-{pick.phase_hint.upper()} | {Mw_sta}")

                    moments.append(M_0)
                    corner_frequencies.extend(corner_freqs)
            
            if not len(moments):
                print (f"No moments could be calculated for event: {event.resource_id}") 
                continue

            # Calculate the seismic moment via basic statistics.
            moments = np.array(moments)
            # moment = moments.mean()
            moment = scipy.stats.trim_mean(moments,0.25)
            moment_std = moments.std()

            corner_frequencies = np.array(corner_frequencies)

            Mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
            Mw_std = 2.0 / 3.0 * moment_std / (moment * np.log(10))
            Mws_std.append(Mw_std)
            Mws.append(Mw)

            print(f"{event.resource_id} | Mw: {round(Mw,2)} | Mw_std {round(Mw_std,2)} | stations:{len(moments)}\n\n")

            mag = ut.write_magnitude_values(Mw,Mw_std,len(moments),"Mw",
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           method_id=ResourceIdentifier("smi:com.github/krischer/moment_magnitude_calculator/automatic/1"),
                           origin_id=event.origins[0].resource_id,
                           agency=self.agency,
                           comments=None)
            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

        if self.xml_mw_out_file != None:
            print (f"Writing output file in {self.xml_mw_out_file}")
            isfile(self.xml_mw_out_file)
            self.catalog.write(self.xml_mw_out_file, out_format)
        return self.catalog
                
    def _get_corresponding_st(self,waveform_id, 
                                pick_time,
                                padding=20.0,
                                ):
        """
        Helper function to find a requested waveform in the previously created
        waveform_index file.
        Also performs the instrument correction.
        Returns None if the file could not be found.
        """
        
        start = pick_time - padding
        end = pick_time + padding

        if (waveform_id.network_code == None) or\
            (waveform_id.network_code == ""):
            # net = "*"
            net = [ x.waveform_restrictions.network for x in self.providers]
            net = ",".join(net)
        else:
            net = waveform_id.network_code
        if (waveform_id.channel_code == None):
            cha = ""
        else:
            cha = waveform_id.channel_code[0:2]

        stream = Stream()
        for provider in self.providers:

            if len(stream) > 0:
                continue

            client = provider.client
            try:
                st= client.get_waveforms(net,
                                        waveform_id.station_code,
                                        "*",
                                        cha+"*",start,end
                                        )
                st = get_st_according2preference(st,
                                    provider.waveform_restrictions.location_preferences,
                                    provider.waveform_restrictions.channel_preferences)
                stream += st
            except:
                pass

        if len(stream) > 0:
            return st
            
        else:
            print(f"\t-> Not found: {net}-{waveform_id.station_code}"+\
                    f"-*-{cha}* |"+\
                        f"{start} - {end}  ")
            return None
        


if __name__ =="__main__":
    import math
    from obspy.clients.fdsn import Client
    from obspy.core.inventory.inventory import read_inventory
    client = 'http://sismo.sgc.gov.co:8080'
    client = Client(client)

    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022knomqj.xml" #3.7
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022cvykgs.xml"
    catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022kszqlt.xml" #2.2 143km
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022ktsrvi.xml" #2 33km
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022krnhiu.xml" #4.3
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022kszqlt.xml" #2.9 103km
    # resp = "/home/emmanuel/Ecopetrol/SeisMonitor/data/metadata/CM.dataless"
    resp = "/home/emmanuel/EDCT/SeisMonitor/data/events/public_CM.xml"
    resp = read_inventory(resp)
    
    out="./test_magnitude,xml"
    mag = Magnitude(client,catalog,resp) 

## prof
    physparams = MwPhysicalMagParams(vp=8200,
                            p_radiation_pattern=0.05,
                            )
    procparams = MwProcessingMagParams(
                            time_before_pick = 2,
                            time_after_pick = 15,
                            only_proc_p_pick=True)

## sup
    # physparams = MwPhysicalMagParams(vp=4800,
    #                         )
    # procparams = MwProcessingMagParams(
    #                         time_before_pick = 0.2,
    #                         time_after_pick = 0.8,
    #                         only_proc_p_pick=True)
    # mag.get_Mw(physparams,procparams,out)

    # ml_params = {"a":1.019,"b":0.0016,"r_ref":140}
    # k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
    #         ml_params["b"]* (ml_params["r_ref"]-100) +3
    # mt = lambda ampl,epi_dist: math.log10(ampl * 1000) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
    #                     ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k
    mag.get_Ml(mag_type="RSNC")