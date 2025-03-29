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
import json
import os
import numpy as np
import time
import logging
import concurrent.futures
from functools import partial
from datetime import timedelta
from obspy.clients.fdsn.mass_downloader.utils import get_mseed_filename
from obspy.core.stream import Stream
from SeisMonitor.utils import printlog,isfile
from . import utils as ut

def run_process(ppc_restrictions,mseed_storage,threshold,to_pick,st):
  """
  Paramaters:
  -----------
  ppc_restrictions: PreprocRestrictions object
    Object to know about preprocessing parameters before downloading
    See the utils file

  mseed_storage: str
    Where to store the waveform files. 

  threshold: int
    limit of length in seconds, length less than threshold will not be downloaded.

  to_pick: tuple
    (batch_size,overlap)
		It's used to know if the stream can be downloaded according to the 
		picker batch size. If the the segments given by the length of the stream
		and overlap parammeter are less than batch_size, then no download the stream.

  st: obspy.Stream object
    Stream that will be written.

  """
  ut.write_stream(st,ppc_restrictions,
                  mseed_storage,threshold,to_pick)

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
    self.providers = providers
    self.new_providers = None
    self._stations_outside_domains = None

  def make_inv_and_json(self, out=None):
    """
    Returns:
    --------
    It saves the json file in the next path:
    {self.info_dir}/json/station_list.json')
    """
    tic = time.time()
    printlog("info",'json',"running to create json")

    # print(inv)
    inv,json_info,uprovider,sod = ut.get_merged_inv_and_json(self.providers)
    self.new_providers = uprovider
    self._stations_outside_domains = sod

    if out != None:
      isfile(out)
      with open(out, 'w') as fp:
        json.dump(json_info, fp)

      toc = time.time()
      exetime = timedelta(seconds=toc-tic)
      printlog("info",'json',
          f'Total time of execution: {exetime.total_seconds()} seconds')  
    return inv,json_info

  def download(self,mseed_storage, 
                batch_size_guarantee = (20,0.3),
                n_processor=1,
                concurrent_feature="thread"):
    """
    Parameters:
    -----------
    
    mseed_storage: str
      Where to store the waveform files.
    batch_size_guarantee : tuple
			(batch_size,overlap)
			It's used to know if the stream can be downloaded according to the 
			picker overlap and batch size parameters. If the the segments given by the length of the stream
			and overlap parammeter are less than batch_size, then no download the stream.
    n_processor: int
      Number of processor
    concurrent_futures: thread
      "thread" or "process"
      Recommend thread
    """
    if self.new_providers == None:
      print("demora")
      self.make_inv_and_json()
    else:
      self.providers = self.new_providers

                
    for provider in self.providers:
      restrictions = provider.download_restrictions
      restrictions.to_pick = batch_size_guarantee 
      if restrictions._name != "DownloadRestrictions":
        raise Exception("Restrictions must be 'DownloadRestrictions' from AIpicker.downloader")

      ppc_restrictions = provider.preproc_restrictions
      self._run_download(mseed_storage=mseed_storage,
                  client=provider.client,dld_restrictions=restrictions,
                  ppc_restrictions=provider.preproc_restrictions,
                  n_processor=n_processor,
                  concurrent_feature=concurrent_feature)

  def _run_download(self,mseed_storage,client,dld_restrictions,
                        ppc_restrictions=None, n_processor=1, 
                        concurrent_feature="thread"):
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

    times = ut.get_chunktimes(starttime=dld_restrictions.starttime,
                          endtime = dld_restrictions.endtime,
                          chunklength_in_sec=dld_restrictions.chunklength,
                          overlap_in_sec=dld_restrictions.overlap_in_sec)
    logger_chunktimes = logging.getLogger('Downloader: chunktime')
    logger_chunktimes.info(f'Total chunktime list: {len(times)}')

    def run_thread(st):
      """
      Download by thread
      """
      ut.write_stream(st,ppc_restrictions,mseed_storage,
                      dld_restrictions.threshold,
                      dld_restrictions.to_pick)
      time.sleep(np.random.randint(25, 30))
    logger_chunkt = logging.getLogger('Downloader: chunktime')  

    for chunkt,(starttime, endtime) in enumerate(times):

      chunktic = time.time()
      logger_chunkt.info(f'{chunkt+1} of {len(times)} -> starttime={starttime} | endtime={endtime}')
      try:
        wavtic = time.time()
        # st = client.get_waveforms_bulk([(dld_restrictions.network,\
        #                           dld_restrictions.station, \
        #                           dld_restrictions.location,\
        #                           dld_restrictions.channel,\
        #                           starttime,\
        #                           endtime)])
        if client.__module__ in ("obspy.clients.filesystem.sds"):
          st = ut.get_all_sdswaveforms(client=client,
                                    network=dld_restrictions.network,
                                    station=dld_restrictions.station, 
                                    location=dld_restrictions.location,
                                    channel=dld_restrictions.channel,
                                    starttime=starttime,
                                    endtime=endtime)
        else:
          st = client.get_waveforms(network=dld_restrictions.network,
                                    station=dld_restrictions.station, 
                                    location=dld_restrictions.location,
                                    channel=dld_restrictions.channel,
                                    starttime=starttime,
                                    endtime=endtime)

        wavtoc = time.time()
        wav_exetime = timedelta(seconds=wavtoc-wavtic)
        if (not dld_restrictions.location_preferences) or \
          (not dld_restrictions.channel_preferences):
          st_dict = st._groupby(dld_restrictions.groupby)
          st_values = list(st_dict.values())
        else:
          st_by_stations = st._groupby("{network}.{station}")
          streams_by_preference = Stream()
          for no, st_by_onestation in enumerate(list(st_by_stations.values())):
            try:
              st_pref = ut.get_st_according2preference(st_by_onestation,
                                dld_restrictions.location_preferences,
                                dld_restrictions.channel_preferences)
            except Exception as e:
              st_pref = st_by_onestation
              station_name = st_pref[0].stats.station
              logger =logging.getLogger(f'Downloader: preference')
              logger.error(f"No select preference. Download all available in {station_name} ")



            if st_pref != None:
              streams_by_preference += st_pref
          st_dict = streams_by_preference._groupby(dld_restrictions.groupby)
          st_values = list(st_dict.values())

      except Exception as e:
        st_warn = (f"{dld_restrictions.network}."
                    f"{dld_restrictions.station}."
                    f"{dld_restrictions.location}."
                    f"{dld_restrictions.channel}."
                    f"{starttime}."
                    f"{endtime}")
        logger = logging.getLogger(f'client.get_waveforms not available: {e}')
        logger.error(f"{st_warn}") 
        st_values = None


      if st_values != None:
        if n_processor == 1:
          for one_st in st_values:
            ut.write_stream(one_st,ppc_restrictions,mseed_storage,
                            dld_restrictions.threshold,
                            dld_restrictions.to_pick)
        else:
          if n_processor > len(st_values):
            n_processor = len(st_values)

          if concurrent_feature in ("thread","Thread","t","T"):
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_processor) as executor:
              executor.map(run_thread,st_values) 

          elif concurrent_feature in ("process","Process","p","P"):
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_processor) as executor:
              myfunc = partial(run_process,ppc_restrictions,mseed_storage,
                              dld_restrictions.threshold,dld_restrictions.to_pick)
              executor.map(myfunc,st_values) 

      chunktoc = time.time()
      chunk_exetime = timedelta(seconds=chunktoc-chunktic)
      exetime_logger = logging.getLogger('Downloader: chunktime')
      exetime_logger.info(f'Time of execution of chunktime {chunkt+1}: {chunk_exetime.total_seconds()} seconds')

    toc = time.time()
    exetime = timedelta(seconds=toc-tic)
    exetime_logger = logging.getLogger('Downloader')
    exetime_logger.info(f'Total time of execution: {exetime.total_seconds()} seconds')

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

