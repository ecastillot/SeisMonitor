import os
from SeisMonitor.monitor.downloader.utils import get_chunktimes
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader
from SeisMonitor.monitor.picker import ai as ai_picker
from SeisMonitor.monitor.associator import ai as ai_asso
from SeisMonitor.monitor.locator.nlloc import nlloc


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

def get_folders_by_chunk(out_folder,starttime,endtime):
    st = starttime.strftime("%Y%m%dT%H%M%S")
    et = endtime.strftime("%Y%m%dT%H%M%S")
    chunk_name = st+"__"+et
    chunk_dir = os.path.join(out_folder,chunk_name)

    downloads_dir = os.path.join(chunk_dir,'downloads')
    metadata_dir = os.path.join(chunk_dir,'metadata')
    detections_dir = os.path.join(chunk_dir,'detections')
    asso_dir = os.path.join(chunk_dir,'associations')
    loc_dir = os.path.join(chunk_dir,'locations')
    events_dir = os.path.join(chunk_dir,'events')

    return {"metadata":metadata_dir,
            "downloads":downloads_dir,
            "detections":detections_dir,
            "associations":asso_dir,
            "locations":loc_dir,
            "events":events_dir}

def sanitize_pick_batch_size(pickers,download_args):
    overlaps = []
    batch_sizes = []
    for picker,args in pickers.items():
        batch_sizes.append(args.batch_size)
        overlaps.append(args.overlap)

    download_args["pick_batch_size"] = (min(overlaps),min(batch_sizes))
    return download_args



class SeisMonitor():
    def __init__(self,providers,out_folder,
                chunklength_in_sec=3600):
        self.providers = providers
        self.out_folder = out_folder
        self.chunklength_in_sec = chunklength_in_sec

        self.download_folder = os.path.join(out_folder,"downloads")
        self.pick_folder = os.path.join(out_folder,"picks")
        self.association_folder = os.path.join(out_folder,"asso")
        self.location_folder = os.path.join(out_folder,"loc")
        self.mag_folder = os.path.join(out_folder,"mag")

        self.process = {}

    def add_downloader(self,
                    download_args={"chunklength_in_sec":None,
                    "threshold": 60,
                    "overlap_in_sec":0,
                    "pick_batch_size": (20,0.3),
                    "groupby":'{network}.{station}.{channel}',
                    "n_processor":None}):
        self.process["download"] = download_args
        
    def add_picker(self,
                    pickers={}):
        if pickers:
            self.process["pick"] = pickers
            if "download" in list(self.process.keys()):
                self.process["download"] = sanitize_pick_batch_size(pickers,self.process["download"])

    # def add_associator(self,
    #                     associators={}):
        
    # def add_associator(self,
    #                 associators={"GaMMA":ai_asso.GaMMAObj(
    #                                     region=[-76.729, -72.315,1.55, 5.314,0, 150],
    #                                     epsg_proj="EPSG:3116",
    #                                     use_amplitude = False,
    #                                     use_dbscan=False,
    #                                     calculate_amp=False) 
    #                              }
    #                 ):

    # def add_locator(self,
    #                 locators={ "NLLoc":nlloc.NLLoc()

    #                         }
    #                 )

    # def add_task(self,
    #             download={"chunklength_in_sec":None,
    #                         "threshold": 60,
    #                         "overlap_in_sec":0,
    #                         "pick_batch_size ": (20,0.3),
    #                         "groupby":'{network}.{station}.{channel}',
    #                         "n_processor":None},
    #             pick={"EQTransformer":EQTransformerObj(model_path = ai_picker.EQTransformer_model_path,
    #                                     n_processor = 4,
    #                                     overlap = 0.3,
    #                                     detection_threshold =0.1,
    #                                     P_threshold = 0.01,
    #                                     S_threshold = 0.01,
    #                                     batch_size = 100,
    #                                     number_of_plots = 0,
    #                                     plot_mode = 1 ) } ,
    #             association,
    #             location,
    #             magnitude):


    def run(self):
        for provider in self.providers:
            restrictions = provider.waveform_restrictions
            starttime = restrictions.starttime
            endtime = restrictions.endtime

            chunktimes = get_chunktimes(starttime,endtime,
                                    self.chunklength_in_sec)

            for starttime,endtime in chunktimes:
                provider.waveform_restrictions.starttime = starttime
                provider.waveform_restrictions.endtime = endtime

                folders = get_folders_by_chunk(self.out_folder,starttime,endtime)

                for process, process_args in self.process.items():

                    if process == "download":
                        structure = os.path.join("{station}","{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
                        download_path = os.path.join(folders["downloads"],structure)
                        md = MseedDownloader([provider])
                        md.make_inv_and_json(folders["metadata"])
                        md.download(download_path,**process_args)

                    if process == "pick":
                        for picker,picker_args in process_args.items():
                            out_path = os.path.join(folders["detections"],picker)
                            if picker == "EQTransformer":
                                _picker = ai_picker.EQTransformer(folders["downloads"],
                                                                folders["metadata"],
                                                                out_path)
                                _picker.pick(picker_args)
                            elif picker == "PhaseNet":
                                _picker = ai_picker.PhaseNet(folders["downloads"],
                                                            folders["metadata"],
                                                            out_path)
                                _picker.mv_downloads2onefolder()
                                _picker.make_datalist()
                                _picker.pick(picker_args)
                                _picker.mv_downloads2stationfolder()



                
        
        
        




        # self.json_path = os.path.join(self.download_folder,"json",
        #                                             "stations.json")
        # self.mseed_folder = os.path.join(self.download_folder,"mseed")
        # if object.name == "EQTransformer":
        #     picker = EQTransformer(self.mseed_storage,
        #                            self.json_path,out_dir)

# if __name__ == "__main__":
