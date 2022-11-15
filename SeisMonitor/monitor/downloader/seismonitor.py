# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2021-12-22 15:13:49
#  * @modify date 2021-12-22 15:13:49
#  * @desc [description]
#  */

"""
Concurrent downloader based on writing mseed files. 

Created on Fri Nov 12 20:00:00 2020
@author: Emmanuel_Castillo
last update: 06-01-2021 
"""
import os
import json
import time
import logging
import concurrent.futures as cf
from datetime import timedelta
from SeisMonitor.utils import printlog,isfile
from . import utils as ut

class MseedDownloader(object):
  def __init__(self,providers):
    """
    Concurrent simple mass downloader.
    

    parameters
    ----------
    providers: list
        list of Client instances
    restrictions: DownloadRestrictions
      Restrictions to download mseed. See utils file.
    ppc_restrictions: PreprocRestrictions object
      Object to know about preprocessing parameters before downloading
      See the utils file

    returns
    -------
      MseedDownloader object

    Warnings:
    ---------
    client instances must have available the get_stations method    

    If the stations are not available. It doesn't download them.
    """
    self.providers = ut.sanitize_provider_times(providers)
    self.providers_are_processed = False
    self._stations_outside_domains = None


  def make_inv_and_json(self, out_folder=None):
    """
    Returns:
    --------
    It saves the json file in the next path:
    {self.info_dir}/json/station_list.json')
    """
    tic = time.time()
    printlog("info",'metadata',"running to create inventory and json files")

    # print(inv)
    inv,json_info,self.providers,sod = ut.get_merged_inv_and_json(self.providers.copy())
    
    self.providers_are_processed = True
    self._stations_outside_domains = sod

    if out_folder != None:
      json_out = os.path.join(out_folder,"stations.json")
      inv_out = os.path.join(out_folder,"inv.xml")

      isfile(inv_out)
      inv.write(inv_out,format="STATIONXML")

      isfile(json_out)
      with open(json_out, 'w') as fp:
        json.dump(json_info, fp)

      toc = time.time()
      exetime = timedelta(seconds=toc-tic)
      printlog("info",'metadata',
          f'Total time of execution: {exetime.total_seconds()} seconds')  
    return inv,json_info

  def download(self,mseed_storage, 
                chunklength_in_sec=None,
                threshold= 60,
                overlap_in_sec=0,
                picker_args = {},
                groupby='{network}.{station}.{channel}',
                n_processor=None):
    """
    Parameters:
    -----------
    
    mseed_storage: str
      Path where to store the waveform files.
      You could use the next key words to set the path for downloading:
      -----------------------
      network,station, location,
      channel, year, month, 
      day, julday,starttime,endtime
    chunklength_in_sec: None or int
			The length of one chunk in seconds. 
			The time between starttime and endtime will be divided 
			into segments of chunklength_in_sec seconds.
		overlap_in_sec: None or int
			For more than one chunk, each segment will have overlapping seconds
		groupby: str
			Download group traces together which have the same metadata given by this parameter. 
			The parameter should name the corresponding keys of the stats object, e.g. '{network}.{station}'. 
			This parameter can take the value 'id' which groups the traces by SEED id.
		threshold: int
			limit of length in seconds, length less than threshold will not be downloaded.
    picker_args : dict
			keys: batch_size,overlap,length
			It's used to know if the stream can be downloaded according to the 
			picker keys. If the the segments given by the length of the stream
			and overlap parammeter are less than batch_size, then no download the stream.
    n_processor: int
      Number of processor
    """
    if not self.providers_are_processed:
      self.make_inv_and_json()
    else:
      pass

                
    for provider in self.providers:
      provider = provider.copy()
      client = provider.client
      waveform_restrictions = provider.waveform_restrictions
      processing = provider.processing
      download_restrictions = ut.DownloadRestrictions(mseed_storage, 
                                                chunklength_in_sec,
                                                threshold,
                                                overlap_in_sec,
                                                picker_args,
                                                groupby,
                                                n_processor)

      self._run_download(client,waveform_restrictions,
                        download_restrictions,
                        processing)
    #   del client; del waveform_restrictions; del download_restrictions; del processing
    # del provider
  
  def _run_download(self,client,waveform_restrictions,
                        download_restrictions,
                        processing=None):
    """
    Parameters:
    -----------
    dld_restrictions: DownloadRestrictions
      Restrictions to download mseed. See utils file.
    mseed_storage: str
      Where to store the waveform files.
    client: Client
      Any obspy client that have get_waveforms function
    ppc_restrictions: PreprocRestrictions object
      Object to know about preprocessing parameters before downloading
      See the utils file
    n_processor: int
      Number of processor
    concurrent_futures: thread
      "thread" or "process"
      Recommend thread
    """
    tic = time.time()

    times = ut.get_chunktimes(starttime=waveform_restrictions.starttime,
                          endtime = waveform_restrictions.endtime,
                          chunklength_in_sec=download_restrictions.chunklength_in_sec,
                          overlap_in_sec=download_restrictions.overlap_in_sec)
    logger_chunktimes = logging.getLogger('Downloader: chunktime')
    logger_chunktimes.info(f'Total chunktime list: {len(times)}')

    logger_chunkt = logging.getLogger('Downloader: chunktime')  

    chunktic = time.time()

    for chunkt,(starttime, endtime) in enumerate(times):

      def get_client_waveforms_by_thread(netsta):
        bulk = (netsta[0],netsta[1],waveform_restrictions.location,
                waveform_restrictions.channel,starttime,endtime)
        ut.write_client_waveforms(client,bulk,
                                waveform_restrictions,
                                download_restrictions,
                                processing  )

      if download_restrictions.n_processor == 1:
        for netsta in waveform_restrictions.bulk_info:
          get_client_waveforms_by_thread(netsta)
      else:
        with cf.ThreadPoolExecutor(download_restrictions.n_processor) as executor:
          executor.map(get_client_waveforms_by_thread,waveform_restrictions.bulk_info)

    chunktoc = time.time()
    wav_exetime = timedelta(seconds=chunktoc-chunktic)
    
    # print(wav_exetime)

if __name__ == "__main__":
  logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] [%(name)s]  %(message)s',
                        datefmt='%m-%d %H:%M')  
  # # set a format which is simpler for console use
  # formatter = logging.Formatter('%(asctime)-15s %(levelname)-8s [%(name)-16s]: %(message)s',
  #                               "%Y-%m-%d %H:%M:%S")
  # # tell the handler to use this format
  # console = logging.StreamHandler()
  # console.setLevel(logging.INFO)
  # console.setFormatter(formatter)
  # #  add the handler to the root logger
  # logging.getLogger().addHandler(console)


  from obspy.clients.fdsn import Client as FDSN_Client
  from obspy.clients.filesystem.sds import Client as SDS_Client
  from obspy.core.utcdatetime import UTCDateTime
  from restrictions import DownloadRestrictions


  
  # client = FDSN_Client('http://sismo.sgc.gov.co:8080')
  client = FDSN_Client('http://10.100.100.232:8091')
  # client = SDS_Client('/mnt/sc232',
  #                    sds_type='D', format='MSEED',)
  
  restrictions = DownloadRestrictions(network="CM",
                          station="AGCC,EZNC,SNPBC,MORC,OCNC,SML1C,VMM*,BRR*,LL*,OCA,PAM,BAR2,PTB,ZAR,RUS,SPBC,NOR,HEL",
                          # location="*",
                          # channel="*",
                          starttime=UTCDateTime("2017-09-22T00:00:00.000000Z"),
                          endtime=UTCDateTime("2017-09-23T00:00:00.000000Z"),
                          chunklength_in_sec=86400,
                          overlap_in_sec=None,
                          groupby='{network}.{station}.{location}',
                          threshold=60,
                          location_preferences=["","00","20","10","40"],
                          channel_preferences=["HH","BH","EH","HN","HL"],
                          to_pick=(10,0.3))

  mseed_storage = ("/home/ecastillo/downloads/"
                  "{network}/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
  
  md = MseedDownloader([client])
  md.download(restrictions,mseed_storage,
                ppc_restrictions=None,
                n_processor=16,concurrent_feature="thread")

