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


from obspy.core.inventory.inventory import Inventory
from obspy.core.stream import Stream
from obspy.io.xseed.parser import Parser
from obspy.core.event import Magnitude as Mag
from obspy.core.event import Catalog
from obspy.core.event import read_events, Comment
from SeisMonitor.utils import isfile
from obspy.core.event.base import CreationInfo
from obspy import UTCDateTime
import mtspec
import numpy as np
import scipy
import os
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

class Magnitude():
    def __init__(self,providers,catalog,out_dir) -> None:
        # self.client = client
        # super().__init__(sds_archive,sds_type,format)
        # self.response = response
        self.providers = providers

        if isinstance(catalog,Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)
        self.out_dir = out_dir
        self.xml_ml_out_file = os.path.join(out_dir,"magnitude","Ml_magnitude.xml")    
        self.xml_mw_out_file = os.path.join(out_dir,"magnitude","Mw_magnitude.xml")    

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

            origin_time = event.origins[0].time
            latitude = event.origins[0].latitude
            longitude = event.origins[0].longitude
            depth = event.origins[0].depth

            if latitude ==None or longitude==None:
                continue

            Mls = []
            for pick in event.picks:
                
                if pick.phase_hint.upper() != "S":
                    continue

                st = self._get_corresponding_st(pick.waveform_id, pick.time,
                                                padding)
                if st == None:
                    continue

                station = st[0].stats.station

                ev_params = {"picktime":pick.time, "latitude":latitude,
                            "longitude":longitude}


                resp = Inventory()
                for provider in self.providers:
                    if len(resp) > 0:
                        continue

                    response = provider.inventory
                    selected_response = response.select(network = pick.waveform_id.network_code,
                                                        station = pick.waveform_id.station_code,
                                                        location = "*",
                                                        channel = pick.waveform_id.channel_code[0:2] +"*" )
                                    
                    resp = selected_response.__add__(selected_response)

                ampl,epi_dist = ut.get_Ml_magparams_by_station(st,resp,
                                ev_params,
                                trimmedtime,
                                waterlevel)
                if ampl == None:
                    continue


                Ml = ut.get_Ml(ampl,epi_dist,mag_type,zone)
                Mls.append(Ml)

                staname = ".".join((pick.waveform_id.network_code, pick.waveform_id.station_code,
                              pick.waveform_id.location_code ,pick.waveform_id.channel_code))
                print(f"\t-> Ml | {staname}-{pick.phase_hint.upper()} | {Ml}")

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
                           method="Rengifo&Ojeda(2004)",
                           origin_id=event.origins[0].resource_id,
                           comments=None)

            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id
            events_mag.append(event)
            print(f"{event.resource_id} | Ml: {round(Ml,2)} | Ml_std {round(Ml_std,2)} | stations:{len(Mls)}\n\n")

        catalog = Catalog(events = events_mag,
                          creation_info= CreationInfo(author="SeisMonitor",
                                            creation_time=UTCDateTime.now())  )
        if self.xml_ml_out_file != None:
            print ("Writing output file...")
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

                    response = provider.inventory
                    selected_response = response.select(network = pick.waveform_id.network_code,
                                                        station = pick.waveform_id.station_code,
                                                        location = "*",
                                                        channel = pick.waveform_id.channel_code[0:2] +"*" )
                                    
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
                              pick.waveform_id.location_code ,pick.waveform_id.channel_code))
                    
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
                           method="smi:com.github/krischer/moment_magnitude_calculator/automatic/1",
                           origin_id=event.origins[0].resource_id,
                           comments=None)
            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

        if self.xml_mw_out_file != None:
            print ("Writing output file...")
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
            net = "CM"
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