from . import utils as ut

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

class EQTransformerObj(object):
    def __init__(self,model_path,
                # chunk_size=3600,
                n_processor=2,overlap=0.3,
                detection_threshold=0.1, P_threshold=0.1,
                S_threshold=0.1,number_of_plots=1,
                batch_size=1,
                plot_mode=1,
                overwrite=False):

        """
        EQTransformer parameters
        """

        self.model_path = model_path
        # self.chunk_size = chunk_size
        self.n_processor = n_processor
        self.overlap = overlap
        self.detection_threshold = detection_threshold
        self.P_threshold = P_threshold
        self.S_threshold = S_threshold
        self.number_of_plots = number_of_plots
        self.batch_size = batch_size
        self.plot_mode = plot_mode
        self.overwrite = overwrite
        self.name = "EQTransformer"

class PhaseNetObj(object):
    def __init__(self, model_path,
                 mode='pred', P_threshold=0.3, S_threshold=0.3,
                batch_size=2, plot = False, save_result=False,
                epochs = 100,learning_rate= 0.01,decay_step = -1,
                decay_rate = 0.9,momentum = 0.9,filters_root = 8,
                depth = 5, kernel_size = [7, 1], pool_size = [4, 1],
                drop_rate = 0, dilation_rate = [1, 1],loss_type = "cross_entropy",
                weight_decay = 0, optimizer = "adam",   summary = True,
                class_weights = [1, 1, 1],log_dir = "log",  num_plots = 10,
                input_length = None,input_mseed = True, filename_picks = "picks",
                data_dir = "./dataset/waveform_pred/",
                data_list = "./dataset/waveform.csv",
                train_dir = "./dataset/waveform_train/",
                train_list ="./dataset/waveform.csv",
                valid_dir = None, valid_list = None,
                output_dir = None ):
        """
        PhaseNet parameters
        """

        self.model_dir = model_path
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
        self.name = "PhaseNet"


class EQTransformer():
    def __init__(self,mseed_storage,json_path,
                out_dir):
        self.mseed_storage = mseed_storage
        self.json_path = json_path
        self.pick_storage = os.path.join(out_dir,"results")
        self.msg_author = "Picker: EQTransfomer"
    
    def picker(self,eqt_obj):
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
        printlog('info',self.msg_author,'Running EQTransformer')
        try:
            tf.compat.v1.reset_default_graph()
            mseed_predictor(
                    input_dir=self.mseed_storage,
                    input_model=eqt_obj.model_path,
                    stations_json=self.json_path,
                    output_dir=self.pick_storage,
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
                printlog('info',self.msg_author,
                f"Error in batch_size: {eqt_obj.batch_size} to the"+
                " last mseed. Check 'to_pick' parameter in DownloadRestrictions")

        tf.keras.backend.clear_session()   
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog('info',self.msg_author,
            f'Total time of execution: {exetime.total_seconds()} seconds')

class PhaseNet():
    def __init__(self,mseed_storage,json_path,
                out_dir):
        self.mseed_storage = mseed_storage
        self.json_path = json_path
        self.pick_storage = os.path.join(out_dir,"results")
        self.datadir = os.path.join(out_dir,"datadir")
        self.datalist = os.path.join(out_dir,"datalist",'fname.csv')
        self.msg_author = "Picker: PhaseNet"

    def mv_downloads2onefolder(self):
        """
        Move all downloads to folder to all_mseed folder.
        """
        tic = time.time()
        # ut.get_one_stream(self.mseed_storage,
        #                 self.datadir)
        ut.mv_mseed2onefolder(self.mseed_storage,self.datadir)
        # ut.mv_mseed2stationfolder(self.datadir,self.mseed_storage)
        # exit()
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Moving mseed2onefolder in {exetime.total_seconds()} seconds.")

    def mv_downloads2stationfolder(self):
        """
        Move all downloads to folder to all_mseed folder.
        """
        tic = time.time()
        ut.mv_mseed2stationfolder(self.datadir,self.mseed_storage)
        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Moving mseed2stationfolder in {exetime.total_seconds()} seconds.")

    def make_datalist(self):
        tic = time.time()
        printlog("info",self.msg_author,"Running to create datalist")
        
        datalist = ut.make_PhaseNet_datalist(datadir=self.datadir, 
                                json_path=self.json_path,
                                datalist_path=self.datalist)

        toc = time.time()
        exetime = dt.timedelta(seconds=toc-tic)
        printlog("info",self.msg_author,f"Datalist created in {exetime.total_seconds()} seconds.")

    def picker(self,pnet_obj):
        tf.compat.v1.reset_default_graph()
        tf.keras.backend.clear_session()

        tic = time.time()
        printlog("info",self.msg_author,"Running PhaseNet pretrained model.")

        args = pnet_obj
        args.data_dir = self.datadir
        args.data_list = self.datalist
        args.output_dir  = self.pick_storage

        # ut.phasenet_from_console(args,self.msg_author)

        datapicks = os.path.join(self.pick_storage,'picks.csv')
        seismonitor_datapicks = os.path.join(self.pick_storage,'seismonitor_picks.csv')
        ut.get_picks(datapicks=datapicks,
                            datalist=self.datalist,
                            min_p_prob=pnet_obj.tp_prob, 
                            min_s_prob=pnet_obj.ts_prob, 
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
    