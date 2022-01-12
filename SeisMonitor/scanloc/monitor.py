# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:48:08
#  * @modify date 2021-12-22 09:48:08
#  * @desc [description]
#  */

"""
EQTransformer seismological monitor
"""

import os
import shutil
import time
import pandas as pd
import logging
from obspy.core.event.catalog import read_events
import tensorflow as tf
from datetime import timedelta
from .utils import makeJSON,get_one_stream,make_PhaseNet_datalist,phasenet_from_console,get_picks,rm_phasenet_duplicate_picks
from EQTransformer.core.mseed_predictor import mseed_predictor
from EQTransformer.utils.associator import run_associator
from SeisMonitor.locator.seisan import Hypocenter,STATION0,resp2df,cp_station0,check_sfile_integrity
from SeisMonitor.downloader.seismonitor import MseedDownloader as mdl
from SeisMonitor.utils import printlog, isfile
# from SeisMonitor.scripts.tools.magnitude import Magnitude


class Monitor(object):
    def __init__(self,providers,restrictions, 
                info_dir, picks_dir, events_dir):
        """
        Parameters:
        -----------
        providers: list
            list of initialized clients. FDSN or SDS Clients.
        restrictions: DownloaderRestrictions 
            To describe the downloader restrictions.
        info_dir: str
            Path where is located the input information to use PhaseNet and
            EQTransformer pickers.
        picks_dir: str
            Path where is located the output information from PhaseNet and
            EQTransformer pickers.
        events_dir: str
            Path where is located the event information after associate and
            locate the event
        """
        self.providers = providers
        self.restrictions = restrictions
        self.info_dir = info_dir
        self.picks_dir = picks_dir
        self.events_dir = events_dir

        ## Input
        self.json_path = os.path.join(self.info_dir,'json','station_list.json')
        self.datalist_dir = os.path.join(self.info_dir,'datalist')
        self.xml_storage = os.path.join(self.info_dir,"xml_storage")
        self.mseed_storage = os.path.join(self.info_dir,'mseed')
        self.all_mseed = os.path.join(self.info_dir,'all_mseed')
        self.hypocenter = os.path.join(self.events_dir,"hypocenter")
        self.sfiles = os.path.join(self.hypocenter,"sfiles")
        self.station0 = os.path.join(self.sfiles,"STATION0.HYP")

        ## Output
        self.eqt_pick_storage = os.path.join(self.picks_dir,'eqt')
        self.phasenet_pick_storage = os.path.join(self.picks_dir,'pnet')

        if os.path.isdir(self.info_dir) == False:
            os.makedirs(self.info_dir)
        if os.path.isdir(self.picks_dir) == False:
            os.makedirs(self.picks_dir)
        if os.path.isdir(self.hypocenter) == False:
            os.makedirs(self.hypocenter)

    def make_json(self, from_xml=None):
        """
        Parameters:
        -----------
        from_xml: str
            Path of xml file that contains inventory information.

            It was thinking because if the provider clients are not able
            to get stations. You could make your json information 
            from one xml file.

        Returns:
        --------
        It saves the json file in the next path:
        {self.info_dir}/json/station_list.json')
        """
        tic = time.time()
        printlog("info",'json',"running to create json")
        directory = os.path.dirname(self.json_path)

        if not os.path.exists(directory):
            os.makedirs(directory)

        json_list = makeJSON(json_path=self.json_path, 
                            providers=self.providers,
                            restrictions=self.restrictions,
                            from_xml=from_xml)
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        printlog("info",'json',
            f'Total time of execution: {exetime.total_seconds()} seconds')

    def download(self,n_processor=1,
        concurrent_feature="t"):
        """
        Parameters:
        ------------
        n_processor: int or None
            Number of threads to download the information
            according to the restriction attribute.
        concurrent_feature: str
            "t"-> threads
            "p"-> Processes

        """
        printlog('info','Downloader','Running to download')
        
        mseed_dl = mdl(providers=self.providers)
        station_storage = os.path.join(self.mseed_storage,
                "{station}","{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
        mseed_dl.download(restrictions=self.restrictions,
                        mseed_storage=station_storage,
                        ppc_restrictions=None,
                        n_processor=n_processor,
                        concurrent_feature=concurrent_feature)

    def mv_downloads2onefolder(self):
        """
        Move all downloads to folder to all_mseed folder.
        """
        get_one_stream(self.mseed_storage,
                        self.all_mseed)

    def make_datalist(self,all_in_folder=True):
        tic = time.time()
        logger = logging.getLogger('PhaseNet: datalist')
        logger.info("Running to create datalist")
        
        if not os.path.exists(self.datalist_dir):
            os.makedirs(self.datalist_dir)

        if all_in_folder:
            datalist = os.path.join(self.datalist_dir,'fname.csv')
            datalist = make_PhaseNet_datalist(datadir=self.all_mseed, 
                                    datalist_path=datalist, 
                                    groupby='{network}.{station}')

        else:

            stations = [ sta for sta in os.listdir(self.mseed_storage)]

            for sta in stations:
                datadir = os.path.join(self.mseed_storage,sta)
                datalist_path = os.path.join(self.datalist_dir,sta,'fname.csv')
                datalist = make_PhaseNet_datalist(datadir=datadir, 
                                        datalist_path=datalist_path, 
                                        groupby='{network}.{station}.{location}')
        
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        logger.info(f'Total time of execution: {exetime.total_seconds()} seconds')

    def phasenet_picker(self,pnet_obj):
        tf.compat.v1.reset_default_graph()
        tf.keras.backend.clear_session()

        tic = time.time()
        pnet_logger = logging.getLogger('PhaseNet')
        pnet_logger.info("Running PhaseNet")

        args = pnet_obj
        args.data_dir = self.all_mseed
        args.data_list = os.path.join(self.datalist_dir,'fname.csv')
        args.output_dir  = self.phasenet_pick_storage

        # run_phasenet.main(args)
        # phasenet_from_console(args)



        picks = os.path.join(self.phasenet_pick_storage,'picks.csv')
        pick2date = get_picks(phaseNet_picks=picks,jsonfile=self.json_path,
                            dt=self.restrictions.chunklength,
                            min_prob=0.3, mode='df_obj', export='csv')

        out_dir = os.path.join(os.path.dirname(picks),'picks_df.csv')
        single_out_dir = os.path.join(os.path.dirname(picks),'picks_df_ok.csv')

        try:
            rm_phasenet_duplicate_picks(path=out_dir, output_path=single_out_dir,
                            group_by_loc=True)
            printlog('info','PhaseNet: cleaner',
                    f'Remove duplicated picks. See your results in {out_dir}')
        except Exception as e:
            printlog('error','PhaseNet: cleaner',
                    f" Can't remove duplicated picks: {e}")
        

        # else: 
        #     stations = [ sta for sta in os.listdir(self.mseed_storage)]
        #     for sta in stations:
        #         static = time.time()
        #         pnet_logger.info(f"picking {sta}")

        #         tf.compat.v1.reset_default_graph()
                        
        #         args.data_dir = os.path.join(self.mseed_storage,sta)
        #         args.data_list = os.path.join(self.datalist_dir,sta,'fname.csv')
        #         args.output_dir  = os.path.join(self.phasenet_pick_storage,sta)

        #         run_phasenet.main(args)
        #         tf.keras.backend.clear_session()

        #         picks = os.path.join(self.phasenet_pick_storage,sta,'picks.csv')
        #         # pick2date = get_picks(phaseNet_picks=picks,jsonfile=self.json_path,
        #         #                 dt=self.restrictions.chunklength,
        #         #                 min_prob=0.3, mode='df_obj', export='csv')
        #         statoc = time.time()
        #         staexetime = timedelta(seconds=statoc-static)
        #         pnet_logger.info(f'{sta} execution: {staexetime.total_seconds()} seconds')

        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        pnet_logger.info(f'Total time of execution: {exetime.total_seconds()} seconds')

    def eqt_picker(self, eqt_obj,rm_downloads=False):
        """
        Parameters:
        --------------
        eqt_obj: ./utils/eqt.EQTobj
            Object where you can set your own
            EQTransformer parameter
        rm_downloads: bool
            True for remove the input mseed folder

        Returns:
        ---------
        It saves your EQtransformer solutions in the next path:
        {self.picks_dir}/eqt
        """

        tic = time.time()
        printlog('info','EQTransformer','Running EQTransformer')
        try:
            tf.compat.v1.reset_default_graph()
            mseed_predictor(
                    input_dir=self.mseed_storage,
                    input_model=eqt_obj.model_path,
                    stations_json=self.json_path,
                    output_dir=self.eqt_pick_storage,
                    detection_threshold=eqt_obj.detection_threshold,
                    P_threshold=eqt_obj.P_threshold, 
                    S_threshold=eqt_obj.S_threshold,
                    number_of_plots=eqt_obj.number_of_plots, 
                    plot_mode=eqt_obj.plot_mode, 
                    overlap=eqt_obj.overlap, 
                    batch_size=eqt_obj.batch_size,
                    overwrite = eqt_obj.overwrite)
        except ValueError as e:
            if str(e) == ("The target structure is of type `<class 'NoneType'>`\n"+
                        "  None\n"+
                        "However the input structure is a sequence (<class 'list'>) of length 0.\n"+
                        "  []\n"+
                        "nest cannot guarantee that it is safe to map one to the other."):
                printlog('info','EQTransformer',
                f"Error in batch_size: {eqt_obj.batch_size} to the"+
                " last mseed. Check 'to_pick' parameter in DownloadRestrictions")

        tf.keras.backend.clear_session()   
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        printlog('info','EQTransformer',
            f'Total time of execution: {exetime.total_seconds()} seconds')

        if rm_downloads:
            shutil.rmtree(self.mseed_storage)
            printlog('info','downloads',
            f'removed: {self.mseed_storage}')

    def eqt_associator(self,moving_window=15, pair_n=3):
        """
        Parameters:
        -----------
        moving_window: int, default=15
            The length of time window used for association in second. 
        
        pair_n: int, default=2
            The minimum number of stations used for the association. 

        Returns:
        --------
        It saves your EQtransformer association solutions in the next path:
        {self.picks_dir}/eqt
        """
        inp_dir = self.eqt_pick_storage
        out_dir = self.eqt_pick_storage
        starttime = self.restrictions.starttime.strftime("%Y-%m-%d %H:%M:%S.%f")
        endtime = self.restrictions.endtime.strftime("%Y-%m-%d %H:%M:%S.%f")
        run_associator(input_dir=inp_dir, output_dir=out_dir ,start_time=starttime, 
                        end_time=endtime,  moving_window=moving_window, 
                        pair_n=pair_n)

    def make_station0(self,xml,vel_model):
        """
        Parameters:
        ----------- 
        xml: str
            Path of xml file that contains inventory information.
        vel_model: str
            Path of csv file that contains the velocity model
            columns -> dep,vp,vs

        returns:
        ---------
        station0: str
            path where was saved the station0.
            Normally is saved in the next path:

            {self.events_dir}/hypocenter/sfiles/STATION0.HYP

        """
        tic = time.time()
        printlog("info",'station0',"running to create station0")
        sta_df = resp2df(xml)
        vel_df = pd.read_csv(vel_model)
        s0 = STATION0(sta_df,vel_df)
        s0.write(self.station0)
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        printlog("info",'station',
                f"{self.station0} ok" +\
                f'Total time of execution: {exetime.total_seconds()} seconds')
        return self.station0

    def hypocenter_locator(self,station0,
                            format="QUAKEML"):   
        """
        Parameters:
        -----------
        station0: str
            Path of the station0 file
        format: str
            Format to export the event locations. (same formats supported by obspy)

        Returns:
        ---------
        Catalog with the event locations with the format indicated in the function.
        It's located in the next path:
            {self.events_dir}/hypocenter/events.xml



        It loads the xml association file that must be located in the next path:
            {self.picks_dir}/eqt

        After it saves it in one nordic event located in the next path:
            {self.events_dir}/hypocenter/sfiles/sfile

        This is helped by the 'split' command (SEISAN) to make all the sfiles
        by each event and it's saved in the next path:
            {self.events_dir}/hypocenter/sfiles/

        Then, it copies the station0 path in the next path:
            {self.events_dir}/hypocenter/sfiles/STATION0.HYP

        Now, it uses the 'update' command (SEISAN) to locate each event and it
        writes the location information in the first lines of each sfile.

        After it collects all the sfiles information in only one sfile, and it's 
        saved in the next path:
            {self.events_dir}/hypocenter/sfiles/collect.out

        Finally, It loads the 'collect.out' nordic file and saves it in the format 
        indicated in the function, and its saved in the next path:
            {self.events_dir}/hypocenter/events.xml

        Warnings:
        ---------
        It removes events not locatable

        """
        tic = time.time()
        printlog("info",'xml_to_nordic',"running to create nordic file")
        asso_xml = os.path.join(self.eqt_pick_storage,"associations.xml") 
        catalog = read_events(asso_xml)
        sfilename = "sfile"
        sfile = os.path.join(self.sfiles,sfilename)
        isfile(sfile)
        catalog.write(sfile,format="NORDIC")
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        printlog("info",'xml_to_nordic',
                f'Total time of execution: {exetime.total_seconds()} seconds')

        if station0 != self.station0:
            cp_station0(station0,self.station0)

        tic = time.time()
        printlog("info",'hypocenter',"running")
        hyp = Hypocenter(self.sfiles,self.station0,from_sfile=sfilename)
        hyp.split()
        hyp.station0()
        hyp.update()
        
        check_sfile_integrity(self.sfiles,rm_not_locatable=True)

        out = os.path.join(self.hypocenter,"events.xml")
        hyp.collect(out,format)
        toc = time.time()
        exetime = timedelta(seconds=toc-tic)
        printlog("info",'hypocenter',
                f'Total time of execution: {exetime.total_seconds()} seconds')

        catalog = read_events(out,format)
        return catalog

    # def compute_magnitude(self,client,catalog,resp,
    #                     method="Mw",out=None):

    #     """
    #     TESTED
    #     """

    #     if os.path.isdir(self.json_path):
    #         jsonfile = self.json_path
    #     else:
    #         jsonfile = None

    #     mag = Magnitude(client,catalog,resp,
    #             jsonfile=jsonfile)
        
    #     if method == "Mw":
    #         cat = mag.estimate_moment_magnitude(out)
    #     elif method == "Ml":
    #         cat = mag.estimate_local_magnitude(out)
    #     else:
    #         raise Exception("Method not defined")

    #     return cat

    