# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 15:13:49
#  * @modify date 2021-12-22 15:13:49
#  * @desc [description]
#  */

"""
Concurrent downloader based on obspy's Mass Downloader class. 

Author:
    Emmanuel Castillo (ecastillot@unal.edu.co), 2020
"""

import time
import numpy as np
import itertools
import concurrent.futures 
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import Inventory
from obspy.clients.fdsn.mass_downloader import (Restrictions,MassDownloader)
from SeisMonitor.utils import printlog

def _run_subprocess(mdl,domain,restriction,
                    mseed_storage,stationxml_storage):
    """
    Parameters
    ----------
    domain: class:'obspy.mass_downloader.domain'
        The download domain.
    restrictions: class:'obspy.mass_downloader.restrictions.Restrictions'
        Non-spatial downloading restrictions.
    mseed_storage: str
        Where to store the waveform files. 
    stationxml_storage: str
        Where to store the StationXML files.

    Returns:
    Download the mseed files.
    """
    try:
        tic = time.time()
        mdl.download(domain=domain,
                restrictions=restriction, 
                download_chunk_size_in_mb=20, 
                threads_per_client=3, 
                mseed_storage=mseed_storage,
                stationxml_storage=stationxml_storage)
        time.sleep(np.random.randint(25, 30))
        toc = time.time()

        printlog("info","Downloader",
        f"Thread done:{restriction.network} -- {restriction.station}"+
        f"\ttime: {round(toc-tic,2)} s")

    except Exception:
        printlog("error","Downloader",
        f"Thread failure: {restriction.network} -- {restriction.station}")
        pass

def process(args):
    """
    Parameters:
    -----------
    args: dict
        Contains the massdownloader parameters to download
        by process. 
    """

    mdl,domain,restriction,\
    mseed_storage,stationxml_storage = args

    try:
        tic = time.time()
        mdl.download(domain=domain,
                restrictions=restriction, 
                download_chunk_size_in_mb=20, 
                threads_per_client=3, 
                mseed_storage=mseed_storage,
                stationxml_storage=stationxml_storage)
        toc = time.time()
        printlog("info","Downloader",
        f"Process done:{restriction.network} -- {restriction.station}"+
        f"\ttime: {round(toc-tic,2)} s")                 

    except Exception:
        printlog("error","Downloader",
        f"Process failure: {restriction.network} -- {restriction.station}")
        pass

class MseedDownloader(object):
    def __init__(self,providers):
        self.providers= providers

    """Concurrent bulk downloader based on obspy's Mass Downloader class.

    parameters
    ----------
    providers: list
        list of Client instances

    returns
    -------
        MseedDownloader object

    Warnings:
    ---------
    client instances must have available the get_stations method    
    """

    def _get_stations_info(self,bulk):
        """
        Parameters:
        -----------
        bulk: str, file or list of lists
            Information about the requested data. 
            See get_stations_bulk from obspy for details.

        returns: 
        --------
            stations_info: list
            list of tuples with stations information
            (network,station)
        """

        inv = Inventory()
        for provider in self.providers:
            one_inv = provider.get_stations_bulk(bulk,level="station")
            # one_inv = provider.get_stations_bulk(bulk,level="response")
            inv = inv.__add__(one_inv) 

        stations_info = [(net.code, sta.code) 
                        for net in inv
                        for sta in net ]
        return stations_info

    def _build_station_restrictions(self,restrictions,
                                    stations_info):
        """
        Parameters:
        -----------
        restrictions: class:'obspy.mass_downloader.restrictions.Restrictions'
            Download restrictions from obspy
        stations_info: list
            list of tuples with stations information
            (network,station)

        Returns:
        --------
        rest_list: list
            List of restrictions by each tuple given by stations_info
        """
        rest_dict = restrictions.__dict__
        new_rest_dict = {}
        for key,value in rest_dict.items():
            if key[0] == '_':
                pass
            elif key == 'chunklength':
                new_rest_dict[key+'_in_sec'] = value
            else:
                new_rest_dict[key] = value

        rest_list = []
        for info in stations_info:
            net,sta = info

            new_rest_dict['network'] = net
            new_rest_dict['station'] = sta
            station_rest = Restrictions(**new_rest_dict)
            rest_list.append(station_rest)

        return rest_list

    def _prepare_args_for_process(self,domain,restrictions_list,
                                    mseed_storage,stationxml_storage):
        
        """
        method to join all massdownloader parameters and save it in 
        the list to map parallelly.
        """
        mdl = MassDownloader(providers=self.providers)
        args = zip(itertools.repeat(mdl), 
            itertools.repeat(domain),
            restrictions_list,
            itertools.repeat(mseed_storage),  
            itertools.repeat(stationxml_storage))
        return list(args)

    def download(self,domain,restrictions,
                    mseed_storage,stationxml_storage,
                    workers=None,parallel_mode="thread"):
        """
        Parallel download by each station. 

        Parameters
        ----------
        domain: class:'obspy.mass_downloader.domain'
            The download domain.
        restrictions: class:'obspy.mass_downloader.restrictions.Restrictions'
            Non-spatial downloading restrictions.
        mseed_storage: str
            Where to store the waveform files. 
        stationxml_storage: str
            Where to store the StationXML files.
        workers: int
            Number of subprocess that will be used.
        parallel_mode: str
            It can be 'thread' or 'process'. However It's recommended
            to use 'process' because by 'thread' you only can use 6 workers
            because for greater workers the downloading is not complete 
            due to no FDSN response.

        returns
        -------
        mseed files and xml files in mseed_storage and stationxml_storage
        respectively.
        """

        ## get the stations information and prepare the restrisctions_list
        ## for the parallel downloading

        bulk = [(restrictions.network, restrictions.station,
                '*', '*',
                 restrictions.starttime, restrictions.endtime) ]
        # bulk = [(restrictions.network, restrictions.station,
        #         restrictions.location, restrictions.channel,
        #          restrictions.starttime, restrictions.endtime) ]

        stations_info = self._get_stations_info(bulk)
        restrictions_list = self._build_station_restrictions(restrictions,
                                                            stations_info)


        mdl = MassDownloader(providers=self.providers)
        ## Go to the downloading: 
        if workers == 1:
            for rest in restrictions_list:
                args = (mdl,domain,rest,mseed_storage,stationxml_storage)
                process(args)
        elif workers == 0:
            raise Exception("workers must be grater than 1")
        else:
            total_tic = time.time()

            #Thread mode
            if parallel_mode in ("thread","t","T"):

                def subprocess(restriction):
                    _run_subprocess(mdl,domain,restriction,
                                    mseed_storage, stationxml_storage)

                with concurrent.futures.ThreadPoolExecutor(
                    max_workers=workers) as executor:
                    executor.map(subprocess,restrictions_list)

            #Process mode
            elif parallel_mode in ("process","p","P"):
                args = self._prepare_args_for_process(domain,
                                    restrictions_list,
                                    mseed_storage,
                                    stationxml_storage)
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=workers) as executor:
                    executor.map(process,args)
                
            else:
                raise Exception(f"Doesn't exist {parallel_mode} mode"
                                "only 1)thread or 2)process")

            total_toc = time.time()
            printlog("info","Downloader"," Total download time: "+
                    f"\t{round(total_toc-total_tic,2)} s")


if __name__ == "__main__":
    from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain
    from obspy.clients.fdsn.client import Client


    IRIS_client = Client(base_url="IRIS", user="gaprietogo@unal.edu.co",
                        password="DaCgmn3hNjg")  

    CM_client = Client(base_url='http://10.100.100.232:8091')

    YU_restrictions = Restrictions(starttime=UTCDateTime(2016, 8, 1),
                            endtime=UTCDateTime(2016, 8, 2),
                            network="YU", 
                            station="CS01,CS02",
                            # location="*", channel="*",
                            chunklength_in_sec=86400,
                            reject_channels_with_gaps=False,
                            minimum_length=0.0,
                            minimum_interstation_distance_in_m=0.0,
                            channel_priorities=["HH[ZNE]", "BH[ZNE]","EH[ZNE]","HN[ZNE]"],
                            location_priorities=["", "00", "20", "10"])

    VMM_restrictions = Restrictions(starttime=UTCDateTime(2016, 8, 1),
                            endtime=UTCDateTime(2016, 8, 2),
                            network="CM", 
                            station="BAR2,VMM*",
                            # location="*", channel="*",
                            chunklength_in_sec=86400,
                            sanitize = False,
                            reject_channels_with_gaps=False,
                            minimum_length=0.0,
                            minimum_interstation_distance_in_m=1000,
                            channel_priorities=["HH[ZNE]", "BH[ZNE]","EH[ZNE]","HN[ZNE]"],
                            location_priorities=["", "00", "20", "10"])

    Colombian_domain = RectangularDomain(minlatitude=-2, maxlatitude=16,
                           minlongitude=-84, maxlongitude=-68)

    mseed_storage = "/home/ecastillo/download/prove/waveforms"
    stationxml_storage = "/home/ecastillo/download/prove/stations"

    mseed_dl = MseedDownloader([CM_client])
    mseed_dl.download(domain=Colombian_domain,
                        restrictions=VMM_restrictions, 
                        mseed_storage=mseed_storage,
                        stationxml_storage=stationxml_storage,
                        workers=1,parallel_mode="process",picker='phasenet')
