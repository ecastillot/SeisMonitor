from genericpath import isdir
from . import utils as ut
import shutil
import os
import time
import logging
import warnings
import pandas as pd
import datetime as dt
import tensorflow as tf
from SeisMonitor.utils import printlog, isfile
from EQTransformer.core.mseed_predictor import mseed_predictor

warnings.simplefilter(action='ignore', category=FutureWarning)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# EQTransformer_core_path = os.path.join(os.path.dirname(__file__),"core","EQTransformer")
# EQTransformer_model_path = os.path.join(EQTransformer_core_path,"ModelsAndSampleData","EqT_model.h5")
# PhaseNet_core_path = os.path.join(os.path.dirname(__file__),"core","PhaseNet")
# PhaseNet_model_path = os.path.join(PhaseNet_core_path,"model","190703-214543")

class EQTransformerObj(object):
    def __init__(self,eqt_path,
                n_processor=2,overlap=0.3,
                detection_threshold=0.1, P_threshold=0.1,
                S_threshold=0.1,number_of_plots=1,
                batch_size=1,
                plot_mode=1,
                overwrite=False,
                # chunk_size=3600,
                rm_downloads=False):

        """
        EQTransformer parameters
        """
        self.eqt_path = eqt_path
        self.model_path = os.path.join(eqt_path,"ModelsAndSampleData","EqT_model.h5")
        self.n_processor = n_processor
        self.overlap = overlap
        self.detection_threshold = detection_threshold
        self.P_threshold = P_threshold
        self.S_threshold = S_threshold
        self.number_of_plots = number_of_plots
        self.batch_size = batch_size
        self.plot_mode = plot_mode
        self.overwrite = overwrite
        self.rm_downloads = rm_downloads
        self.name = "EQTransformer"

class PhaseNetObj(object):
    def __init__(self, pnet_path,
                 mode='pred', P_threshold=0.3, S_threshold=0.3,
                batch_size=2, plot = False, save_result=False,
                epochs = 100,learning_rate= 0.01,decay_step = -1,
                decay_rate = 0.9,momentum = 0.9,filters_root = 8,
                depth = 5, kernel_size = [7, 1], pool_size = [4, 1],
                drop_rate = 0, dilation_rate = [1, 1],loss_type = "cross_entropy",
                weight_decay = 0, optimizer = "adam",   summary = True,
                class_weights = [1, 1, 1],log_dir = "log",  num_plots = 10,
                input_length = None,input_mseed = True, filename_picks = "picks",
                one_single_sampling_rate=-1,
                data_dir = "./dataset/waveform_pred/",
                data_list = "./dataset/waveform.csv",
                train_dir = "./dataset/waveform_train/",
                train_list ="./dataset/waveform.csv",
                valid_dir = None, valid_list = None,
                output_dir = None,
                rm_downloads=False ):
        """
        PhaseNet parameters
        """

        self.phasenet_path = pnet_path
        self.model_dir = os.path.join(pnet_path,"model","190703-214543")
        self.mode = mode
        self.tp_prob = P_threshold
        self.ts_prob = S_threshold
        self.batch_size = batch_size
        self.plot_figure = plot
        self.save_result = save_result
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.decay_step = decay_step
        self.decay_rate = decay_rate
        self.momentum = momentum
        self.filters_root = filters_root
        self.depth = depth
        self.kernel_size = kernel_size
        self.pool_size = pool_size
        self.drop_rate = drop_rate
        self.dilation_rate = dilation_rate
        self.loss_type = loss_type
        self.weight_decay = weight_decay
        self.optimizer = optimizer
        self.summary = summary
        self.class_weights = class_weights
        self.log_dir = log_dir
        self.num_plots = num_plots
        self.input_length = input_length
        self.input_mseed = input_mseed
        self.fpred = filename_picks
        self.data_dir = data_dir
        self.data_list = data_list
        self.train_dir = train_dir
        self.train_list = train_list
        self.valid_dir = valid_dir
        self.valid_list = valid_list
        self.output_dir = output_dir
        self.overlap = 0.5
        self.rm_downloads = rm_downloads
        self.one_single_sampling_rate = one_single_sampling_rate
        self.name = "PhaseNet"

class EQTransformer():
    def __init__(self,eqt_obj):
        self.eqt_obj = eqt_obj
    
    def pick(self,mseed_storage,metadata_dir,
                out_dir):
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
        self.mseed_storage = mseed_storage
        self.json_path = os.path.join(metadata_dir,"stations.json")
        self.pick_storage = os.path.join(out_dir,"results")
        self.msg_author = "Picker: EQTransfomer"

        tic = time.time()
        printlog('info',self.msg_author,'Running EQTransformer')
        try:
            tf.compat.v1.reset_default_graph()
            mseed_predictor(
                    input_dir=self.mseed_storage,
                    input_model=self.eqt_obj.model_path,
                    stations_json=self.json_path,
                    output_dir=self.pick_storage,
                    detection_threshold=self.eqt_obj.detection_threshold,
                    P_threshold=self.eqt_obj.P_threshold, 
                    S_threshold=self.eqt_obj.S_threshold,
                    number_of_plots=self.eqt_obj.number_of_plots, 
                    plot_mode=self.eqt_obj.plot_mode, 
                    overlap=self.eqt_obj.overlap, 
                    batch_size=self.eqt_obj.batch_size,
                    overwrite = self.eqt_obj.overwrite)
        except ValueError as e:
            if str(e) == ("The target structure is of type `<class 'NoneType'>`\n"+
                        "  None\n"+
                        "However the input structure is a sequence (<class 'list'>) of length 0.\n"+
                        "  []\n"+
                        "nest cannot guarantee that it is safe to map one to the other."):
                printlog('info',self.msg_author,
                f"Error in batch_size: {self.eqt_obj.batch_size} to the"+
                " last mseed. Check 'to_pick' parameter in DownloadRestrictions")
        # eqt_merge = os.path.join(self.pick_storage,"eqt_merge.csv")
        # ut.merge_picks(self.pick_storage,"X_prediction_results.csv",
                    #    eqt_merge)
        output_path = os.path.join(self.pick_storage,"seismonitor_picks.csv")
        result = ut.eqt_picks_2_seismonitor_fmt(self.pick_storage,self.mseed_storage,
                                        output_path)
        tf.keras.backend.clear_session()   
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog('info',self.msg_author,
            f'Total time of execution: {exetime.total_seconds()} seconds')

        if self.eqt_obj.rm_downloads:
            shutil.rmtree(self.mseed_storage, ignore_errors=True)
        return result

class PhaseNet():
    def __init__(self,pnet_obj):
        self.pnet_obj = pnet_obj
        self.msg_author = "Picker: PhaseNet"

    def __mv_downloads2onefolder(self):
        """
        Move all downloads to folder to all_mseed folder.
        """
        tic = time.time()
        # ut.get_one_stream(self.mseed_storage,
        #                 self.datadir)
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)

        ut.mv_mseed2onefolder(self.mseed_storage,self.datadir)
        # ut.mv_mseed2stationfolder(self.datadir,self.mseed_storage)
        # exit()
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Moving mseed2onefolder in {exetime.total_seconds()} seconds.")

    def __mv_downloads2stationfolder(self):
        """
        Move all downloads to folder to all_mseed folder.
        """
        tic = time.time()
        ut.mv_mseed2stationfolder(self.datadir,self.mseed_storage)
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Moving mseed2stationfolder in {exetime.total_seconds()} seconds.")

    def __make_datalist(self):
        tic = time.time()
        printlog("info",self.msg_author,"Running to create datalist")
        
        datalist = ut.make_PhaseNet_datalist(datadir=self.datadir, 
                                json_path=self.json_path,
                                datalist_path=self.datalist)

        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Datalist created in {exetime.total_seconds()} seconds.")

    def pick(self,mseed_storage,metadata_dir,
                out_dir):
        self.mseed_storage = mseed_storage
        self.json_path = os.path.join(metadata_dir,"stations.json")
        self.pick_storage = os.path.join(out_dir,"results")
        self.datadir = os.path.join(out_dir,"datadir")
        self.datalist = os.path.join(out_dir,"datalist",'fname.csv')
        
        
        self.__mv_downloads2onefolder()
        self.__make_datalist()

        tf.compat.v1.reset_default_graph()
        tf.keras.backend.clear_session()

        tic = time.time()
        printlog("info",self.msg_author,"Running PhaseNet pretrained model.")

        args = self.pnet_obj
        args.data_dir = self.datadir
        args.data_list = self.datalist
        args.output_dir  = self.pick_storage

        ut.phasenet_from_console(args,self.msg_author)
        datapicks = os.path.join(self.pick_storage,'picks.csv')
        seismonitor_datapicks = os.path.join(self.pick_storage,'seismonitor_picks.csv')
        result = ut.get_picks(datapicks=datapicks,
                            datalist=self.datalist,
                            min_p_prob=self.pnet_obj.tp_prob, 
                            min_s_prob=self.pnet_obj.ts_prob, 
                            one_single_sampling_rate = self.pnet_obj.one_single_sampling_rate,
                             mode='df_obj',
                            export=seismonitor_datapicks)

        try:
            ut.rm_phasenet_duplicate_picks(path=seismonitor_datapicks, 
                                        output_path=seismonitor_datapicks)
            printlog('info',self.msg_author,
                    f'Duplicated picks had been removed. See your results in {seismonitor_datapicks}')
        except Exception as e:
            printlog('error',self.msg_author,
                    f" Can't remove duplicated picks: {e}")
        
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f'Time of prediction: {exetime.total_seconds()} seconds')

        self.__mv_downloads2stationfolder()        
        if self.pnet_obj.rm_downloads:
            shutil.rmtree(self.mseed_storage, ignore_errors=True)

        return result