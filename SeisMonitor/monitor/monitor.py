import os
from SeisMonitor.monitor.downloader.utils import get_chunktimes
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj

def get_chunktimes_by_provider(providers,chunklength_in_sec):
    times_by_provider = {}
    for provider in providers:
        restrictions = provider.waverform_restrictions
        starttime = restrictions.starttime
        endtime = restrictions.endtime

        times = get_chunktimes(starttime,endtime,
                                chunklength_in_sec)
        times_by_provider.append(times)
    return times_by_provider

class SeisMonitor():
    def __init__(self,providers,out_folder,
                chunklength_in_sec=7200):
        self.providers = providers
        self.out_folder = out_folder
        self.chunklength_in_sec = chunklength_in_sec

        self.download_folder = os.path.join(out_folder,"downloads")
        self.pick_folder = os.path.join(out_folder,"picks")
        self.association_folder = os.path.join(out_folder,"asso")
        self.location_folder = os.path.join(out_folder,"loc")
        self.mag_folder = os.path.join(out_folder,"mag")

    def download(self,
                chunklength_in_sec=None,
                threshold= 60,
                overlap_in_sec=0,
                pick_batch_size = (20,0.3),
                groupby='{network}.{station}.{channel}',
                n_processor=None):


        md = MseedDownloader(self.providers)
        md.make_inv_and_json(self.json_path)
        md.download(self.mseed_storage,
                    chunklength_in_sec,
                    threshold,
                    overlap_in_sec,
                    pick_batch_size,
                    groupby,n_processor)

    def run(self,object):
        for provider in self.providers:
            restrictions = provider.waverform_restrictions
            starttime = restrictions.starttime
            endtime = restrictions.endtime

            chunktimes = get_chunktimes(starttime,endtime,
                                    self.chunklength_in_sec)

            for starttime,endtime in chunktimes:
                provider.waverform_restrictions.starttime = starttime
                provider.waverform_restrictions.endtime = endtime

                md = MseedDownloader([provider])

                st = starttime.strftime("%Y%m%dT%H%M%S")
                et = endtime.strftime("%Y%m%dT%H%M%S")
                chunk_name = st+"__"+et
                chunk_dir = os.path.join(self.out_folder,chunk_name)
                info_dir = os.path.join(chunk_dir,'downloads')
                picks_dir = os.path.join(chunk_dir,'detections')
                events_dir = os.path.join(chunk_dir,'events')
                md.make_inv_and_json()

                
        
        
        




        # self.json_path = os.path.join(self.download_folder,"json",
        #                                             "stations.json")
        # self.mseed_folder = os.path.join(self.download_folder,"mseed")
        # if object.name == "EQTransformer":
        #     picker = EQTransformer(self.mseed_storage,
        #                            self.json_path,out_dir)

if __name__ == "__main__":
