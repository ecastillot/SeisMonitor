import time
import numpy as np
import itertools
import concurrent.futures
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.inventory import Inventory
from obspy.clients.fdsn.mass_downloader import Restrictions, MassDownloader
from SeisMonitor.utils import printlog


def _run_subprocess(mdl, domain, restriction, mseed_storage, stationxml_storage):
    """Execute download in a subprocess with error handling.
    
    Args:
        mdl (MassDownloader): MassDownloader instance
        domain (Domain): Download domain specification
        restriction (Restrictions): Non-spatial download restrictions
        mseed_storage (str): Path for waveform storage
        stationxml_storage (str): Path for StationXML storage
    """
    try:
        tic = time.time()
        mdl.download(
            domain=domain,
            restrictions=restriction,
            download_chunk_size_in_mb=20,
            threads_per_client=3,
            mseed_storage=mseed_storage,
            stationxml_storage=stationxml_storage
        )
        time.sleep(np.random.randint(25, 30))
        toc = time.time()
        printlog(
            "info", "Downloader",
            f"Thread done: {restriction.network} -- {restriction.station}"
            f"\ttime: {round(toc-tic, 2)} s"
        )
    except Exception:
        printlog(
            "error", "Downloader",
            f"Thread failure: {restriction.network} -- {restriction.station}"
        )


def process(args):
    """Process download parameters for a single station.
    
    Args:
        args (tuple): Contains (mdl, domain, restriction, mseed_storage, stationxml_storage)
            mdl (MassDownloader): MassDownloader instance
            domain (Domain): Download domain specification
            restriction (Restrictions): Non-spatial download restrictions
            mseed_storage (str): Path for waveform storage
            stationxml_storage (str): Path for StationXML storage
    """
    mdl, domain, restriction, mseed_storage, stationxml_storage = args
    
    try:
        tic = time.time()
        mdl.download(
            domain=domain,
            restrictions=restriction,
            download_chunk_size_in_mb=20,
            threads_per_client=3,
            mseed_storage=mseed_storage,
            stationxml_storage=stationxml_storage
        )
        toc = time.time()
        printlog(
            "info", "Downloader",
            f"Process done: {restriction.network} -- {restriction.station}"
            f"\ttime: {round(toc-tic, 2)} s"
        )
    except Exception:
        printlog(
            "error", "Downloader",
            f"Process failure: {restriction.network} -- {restriction.station}"
        )


class MseedDownloader:
    """Concurrent bulk downloader based on Obspy's Mass Downloader class.
    
    Attributes:
        providers (list): List of FDSN client instances
        
    Warnings:
        Client instances must implement the get_stations method
    """
    
    def __init__(self, providers):
        """Initialize MseedDownloader with FDSN providers.
        
        Args:
            providers (list): List of Client instances
        """
        self.providers = providers

    def _get_stations_info(self, bulk):
        """Retrieve station information from providers.
        
        Args:
            bulk (str or list): Station request information
                See obspy.clients.fdsn.client.Client.get_stations_bulk for details
                
        Returns:
            list: List of tuples with (network, station) information
        """
        inv = Inventory()
        for provider in self.providers:
            one_inv = provider.get_stations_bulk(bulk, level="station")
            inv = inv.__add__(one_inv)

        stations_info = [(net.code, sta.code)
                        for net in inv
                        for sta in net]
        return stations_info

    def _build_station_restrictions(self, restrictions, stations_info):
        """Build restriction objects for each station.
        
        Args:
            restrictions (Restrictions): Base download restrictions
            stations_info (list): List of (network, station) tuples
            
        Returns:
            list: List of Restrictions objects per station
        """
        rest_dict = {k: v for k, v in restrictions.__dict__.items() if k[0] != '_'}
        if 'chunklength' in rest_dict:
            rest_dict['chunklength_in_sec'] = rest_dict.pop('chunklength')

        rest_list = []
        for net, sta in stations_info:
            new_rest_dict = rest_dict.copy()
            new_rest_dict['network'] = net
            new_rest_dict['station'] = sta
            rest_list.append(Restrictions(**new_rest_dict))

        return rest_list

    def _prepare_args_for_process(self, domain, restrictions_list, mseed_storage, stationxml_storage):
        """Prepare arguments for parallel processing.
        
        Args:
            domain (Domain): Download domain specification
            restrictions_list (list): List of Restrictions objects
            mseed_storage (str): Path for waveform storage
            stationxml_storage (str): Path for StationXML storage
            
        Returns:
            list: List of argument tuples for parallel execution
        """
        mdl = MassDownloader(providers=self.providers)
        args = zip(
            itertools.repeat(mdl),
            itertools.repeat(domain),
            restrictions_list,
            itertools.repeat(mseed_storage),
            itertools.repeat(stationxml_storage)
        )
        return list(args)

    def download(self, domain, restrictions, mseed_storage, stationxml_storage, workers=None, parallel_mode="thread"):
        """Perform parallel download of seismic data by station.
        
        Args:
            domain (Domain): Download domain specification
            restrictions (Restrictions): Non-spatial download restrictions
            mseed_storage (str): Path for waveform storage
            stationxml_storage (str): Path for StationXML storage
            workers (int, optional): Number of parallel workers
            parallel_mode (str): 'thread' or 'process' parallel execution mode
                'thread' recommended for <=6 workers due to FDSN response limitations
                'process' recommended for higher worker counts
                
        Raises:
            Exception: If workers < 1 or invalid parallel_mode specified
            
        Notes:
            Downloads mseed and StationXML files to specified storage locations
        """
        bulk = [(restrictions.network, restrictions.station, '*', '*',
                restrictions.starttime, restrictions.endtime)]
        
        stations_info = self._get_stations_info(bulk)
        restrictions_list = self._build_station_restrictions(restrictions, stations_info)
        mdl = MassDownloader(providers=self.providers)

        if workers == 1:
            for rest in restrictions_list:
                process((mdl, domain, rest, mseed_storage, stationxml_storage))
        elif workers == 0 or workers is None:
            raise Exception("workers must be greater than 0")
        else:
            total_tic = time.time()

            if parallel_mode.lower() in ("thread", "t"):
                with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
                    executor.map(
                        lambda r: _run_subprocess(mdl, domain, r, mseed_storage, stationxml_storage),
                        restrictions_list
                    )
            elif parallel_mode.lower() in ("process", "p"):
                args = self._prepare_args_for_process(
                    domain, restrictions_list, mseed_storage, stationxml_storage
                )
                with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                    executor.map(process, args)
            else:
                raise Exception(f"Invalid parallel_mode '{parallel_mode}' - use 'thread' or 'process'")

            total_toc = time.time()
            printlog(
                "info", "Downloader",
                f"Total download time: {round(total_toc-total_tic, 2)} s"
            )


if __name__ == "__main__":
    from obspy.clients.fdsn.mass_downloader.domain import RectangularDomain
    from obspy.clients.fdsn.client import Client

    IRIS_client = Client(
        base_url="IRIS",
        user="gaprietogo@unal.edu.co",
        password="DaCgmn3hNjg"
    )
    CM_client = Client(base_url='http://10.100.100.232:8091')

    YU_restrictions = Restrictions(
        starttime=UTCDateTime(2016, 8, 1),
        endtime=UTCDateTime(2016, 8, 2),
        network="YU",
        station="CS01,CS02",
        chunklength_in_sec=86400,
        reject_channels_with_gaps=False,
        minimum_length=0.0,
        minimum_interstation_distance_in_m=0.0,
        channel_priorities=["HH[ZNE]", "BH[ZNE]", "EH[ZNE]", "HN[ZNE]"],
        location_priorities=["", "00", "20", "10"]
    )

    VMM_restrictions = Restrictions(
        starttime=UTCDateTime(2016, 8, 1),
        endtime=UTCDateTime(2016, 8, 2),
        network="CM",
        station="BAR2,VMM*",
        chunklength_in_sec=86400,
        sanitize=False,
        reject_channels_with_gaps=False,
        minimum_length=0.0,
        minimum_interstation_distance_in_m=1000,
        channel_priorities=["HH[ZNE]", "BH[ZNE]", "EH[ZNE]", "HN[ZNE]"],
        location_priorities=["", "00", "20", "10"]
    )

    Colombian_domain = RectangularDomain(
        minlatitude=-2,
        maxlatitude=16,
        minlongitude=-84,
        maxlongitude=-68
    )

    mseed_storage = "/home/ecastillo/download/prove/waveforms"
    stationxml_storage = "/home/ecastillo/download/prove/stations"

    mseed_dl = MseedDownloader([CM_client])
    mseed_dl.download(
        domain=Colombian_domain,
        restrictions=VMM_restrictions,
        mseed_storage=mseed_storage,
        stationxml_storage=stationxml_storage,
        workers=1,
        parallel_mode="process"
    )