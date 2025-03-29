import os
import json
import time
import logging
import concurrent.futures as cf
from datetime import timedelta
from SeisMonitor.utils import printlog, isfile
from . import utils as ut


class MseedDownloader:
    """Concurrent simple mass downloader for seismic data.
    
    Attributes:
        providers (list): List of processed Client instances
        providers_are_processed (bool): Flag indicating if providers are processed
        _stations_outside_domains (set): Stations outside requested domains
        
    Warnings:
        Client instances must implement the get_stations method
        Stations not available will not be downloaded
    """
    
    def __init__(self, providers):
        """Initialize MseedDownloader with client providers.
        
        Args:
            providers (list): List of Client instances
        """
        self.providers = ut.sanitize_provider_times(providers)
        self.providers_are_processed = False
        self._stations_outside_domains = None

    def make_inv_and_json(self, out_folder=None):
        """Create inventory and JSON files from provider metadata.
        
        Args:
            out_folder (str, optional): Directory to save output files
            
        Returns:
            tuple: (Inventory, dict) containing station inventory and JSON info
            
        Notes:
            If out_folder is provided, saves:
            - stations.json at {out_folder}/stations.json
            - inv.xml at {out_folder}/inv.xml
        """
        tic = time.time()
        printlog("info", "metadata", "running to create inventory and json files")

        inv, json_info, self.providers, sod = ut.get_merged_inv_and_json(self.providers.copy())
        self.providers_are_processed = True
        self._stations_outside_domains = sod

        if out_folder:
            json_out = os.path.join(out_folder, "stations.json")
            inv_out = os.path.join(out_folder, "inv.xml")

            isfile(inv_out)
            inv.write(inv_out, format="STATIONXML")

            isfile(json_out)
            with open(json_out, 'w') as fp:
                json.dump(json_info, fp)

            toc = time.time()
            exetime = timedelta(seconds=toc - tic)
            printlog(
                "info", "metadata",
                f"Total time of execution: {exetime.total_seconds()} seconds"
            )
        
        return inv, json_info

    def download(
        self,
        mseed_storage,
        chunklength_in_sec=None,
        threshold=60,
        overlap_in_sec=0,
        picker_args={},
        groupby='{network}.{station}.{channel}',
        n_processor=None
    ):
        """Download seismic waveforms with specified parameters.
        
        Args:
            mseed_storage (str): Path template for waveform storage
                Supports keywords: {network}, {station}, {location}, {channel},
                {year}, {month}, {day}, {julday}, {starttime}, {endtime}
            chunklength_in_sec (int, optional): Length of each time chunk in seconds
            threshold (int): Minimum length in seconds for download
            overlap_in_sec (int): Overlap between chunks in seconds
            picker_args (dict): Picker parameters (batch_size, overlap, length)
            groupby (str): Grouping pattern for traces (e.g., '{network}.{station}')
            n_processor (int, optional): Number of parallel processors
            
        Notes:
            If providers aren't processed, triggers make_inv_and_json()
        """
        if not self.providers_are_processed:
            self.make_inv_and_json()

        for provider in self.providers:
            provider = provider.copy()
            client = provider.client
            waveform_restrictions = provider.waveform_restrictions
            processing = provider.processing
            download_restrictions = ut.DownloadRestrictions(
                mseed_storage,
                chunklength_in_sec,
                threshold,
                overlap_in_sec,
                picker_args,
                groupby,
                n_processor
            )

            self._run_download(
                client,
                waveform_restrictions,
                download_restrictions,
                processing
            )

    def _run_download(self, client, waveform_restrictions, download_restrictions, processing=None):
        """Execute the download process for a single client.
        
        Args:
            client (Client): Obspy client with get_waveforms method
            waveform_restrictions (WaveformRestrictions): Waveform constraints
            download_restrictions (DownloadRestrictions): Download parameters
            processing (list, optional): Processing steps to apply
        """
        tic = time.time()
        times = ut.get_chunktimes(
            starttime=waveform_restrictions.starttime,
            endtime=waveform_restrictions.endtime,
            chunklength_in_sec=download_restrictions.chunklength_in_sec,
            overlap_in_sec=download_restrictions.overlap_in_sec
        )
        
        logger_chunktimes = logging.getLogger('Downloader: chunktime')
        logger_chunktimes.info(f'Total chunktime list: {len(times)}')

        chunktic = time.time()
        for chunkt, (starttime, endtime) in enumerate(times):
            def get_client_waveforms_by_thread(netsta):
                """Download waveforms for a single network-station pair."""
                bulk = (
                    netsta[0], netsta[1],
                    waveform_restrictions.location,
                    waveform_restrictions.channel,
                    starttime, endtime
                )
                ut.write_client_waveforms(
                    client, bulk,
                    waveform_restrictions,
                    download_restrictions,
                    processing
                )

            if download_restrictions.n_processor == 1:
                for netsta in waveform_restrictions.bulk_info:
                    get_client_waveforms_by_thread(netsta)
            else:
                with cfitools.ThreadPoolExecutor(download_restrictions.n_processor) as executor:
                    executor.map(
                        get_client_waveforms_by_thread,
                        waveform_restrictions.bulk_info
                    )

        chunktoc = time.time()
        wav_exetime = timedelta(seconds=chunktoc - chunktic)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] [%(name)s]  %(message)s',
        datefmt='%m-%d %H:%M'
    )

    from obspy.clients.fdsn import Client as FDSN_Client
    from obspy.core.utcdatetime import UTCDateTime
    from restrictions import DownloadRestrictions

    client = FDSN_Client('http://10.100.100.232:8091')
    
    restrictions = DownloadRestrictions(
        network="CM",
        station="AGCC,EZNC,SNPBC,MORC,OCNC,SML1C,VMM*,BRR*,LL*,OCA,PAM,BAR2,PTB,ZAR,RUS,SPBC,NOR,HEL",
        starttime=UTCDateTime("2017-09-22T00:00:00.000000Z"),
        endtime=UTCDateTime("2017-09-23T00:00:00.000000Z"),
        chunklength_in_sec=86400,
        overlap_in_sec=None,
        groupby='{network}.{station}.{location}',
        threshold=60,
        location_preferences=["", "00", "20", "10", "40"],
        channel_preferences=["HH", "BH", "EH", "HN", "HL"],
        to_pick=(10, 0.3)
    )

    mseed_storage = (
        "/home/ecastillo/downloads/"
        "{network}/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
    )
    
    md = MseedDownloader([client])
    md.download(
        restrictions,
        mseed_storage,
        n_processor=16,
        concurrent_feature="thread"
    )