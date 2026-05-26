"""
Phase picking interfaces for EQTransformer and PhaseNet.

This module provides wrapper classes and configuration containers for
running seismic phase picking workflows using:

- EQTransformer
- PhaseNet

The classes standardize:
    * model configuration
    * waveform preprocessing
    * prediction execution
    * result conversion to SeisMonitor format

The module also includes helper utilities for:
    * TensorFlow session cleanup
    * MiniSEED organization
    * datalist generation
    * duplicate pick removal

Notes
-----
TensorFlow warnings and verbose logging are suppressed to reduce
console noise during inference.

Examples
--------
Run EQTransformer::

    eqt_config = EQTransformerObj(
        model_path="model.h5",
        batch_size=4,
    )

    picker = EQTransformer(eqt_config)

    picks = picker.pick(
        mseed_storage="./mseed",
        metadata_dir="./metadata",
        out_dir="./results",
    )

Run PhaseNet::

    pnet_config = PhaseNetObj(
        pnet_path="./PhaseNet",
        model_path="./model",
    )

    picker = PhaseNet(pnet_config)

    picks = picker.pick(
        mseed_storage="./mseed",
        metadata_dir="./metadata",
        out_dir="./results",
    )
"""

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
from git import Repo
from SeisMonitor.utils import printlog, isfile
from EQTransformer.core.mseed_predictor import mseed_predictor

# Suppress warnings and TensorFlow logging
warnings.simplefilter(action='ignore', category=FutureWarning)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)


class EQTransformerObj:
    """
    Configuration container for EQTransformer predictions.

    Parameters
    ----------
    model_path : str
        Path to the trained EQTransformer model.

    n_processor : int, default=2
        Number of processors used during prediction.

    overlap : float, default=0.3
        Overlap fraction between prediction windows.

    detection_threshold : float, default=0.1
        Detection probability threshold.

    P_threshold : float, default=0.1
        P-phase probability threshold.

    S_threshold : float, default=0.1
        S-phase probability threshold.

    number_of_plots : int, default=1
        Number of prediction plots to generate.

    batch_size : int, default=1
        Number of waveforms processed per batch.

    plot_mode : int, default=1
        Plotting mode used by EQTransformer.

    overwrite : bool, default=False
        Whether existing outputs should be overwritten.

    rm_downloads : bool, default=False
        Remove MiniSEED files after prediction.

    Attributes
    ----------
    name : str
        Picker identifier.
    """

    def __init__(
        self,
        model_path,
        n_processor=2,
        overlap=0.3,
        detection_threshold=0.1,
        P_threshold=0.1,
        S_threshold=0.1,
        number_of_plots=1,
        batch_size=1,
        plot_mode=1,
        overwrite=False,
        rm_downloads=False
    ):
        """Initialize EQTransformer configuration."""
        self.model_path = model_path
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


class PhaseNetObj:
    """
    Configuration container for PhaseNet predictions.

    Parameters
    ----------
    pnet_path : str
        Path to PhaseNet repository.
    model_path : str
        Path to trained PhaseNet model.
    mode : str, default="pred"
        Operation mode ('pred' for prediction).
    P_threshold : float, default=0.3
        Minimum P-phase probability.
    S_threshold : float, default=0.3
        Minimum S-phase probability.
    batch_size : int, default=2
        Number of waveforms per batch.
    plot : bool, default=False
        Whether to generate plots.
    save_result : bool, default=False
        Whether to save prediction results.
    epochs : int, default=100
        Number of training epochs.
    learning_rate : float, default=0.01
        Initial learning rate.
    decay_step : int, default=100
        Steps for learning rate decay.
    decay_rate : float, default=0.9
        Learning rate decay rate.
    momentum : float, default=0.9
        Momentum for optimizer.
    filters_root : int, default=8
        Base number of filters.
    depth : int, default=5
        Network depth.
    kernel_size : list, default=[7, 5, 5, 5, 5]
        Kernel size for convolution.
    pool_size : list, default=[4, 2, 2, 2, 2]
        Pooling size.
    drop_rate : float, default=0.1
        Dropout rate.
    dilation_rate : list, default=[1, 2, 4, 8, 16]
        Dilation rate for convolution.
    loss_type : str, default="dice"
        Loss function type.
    weight_decay : float, default=1e-5
        Weight decay factor.
    optimizer : str, default="adam"
        Optimizer type.
    summary : bool, default=True
        Whether to show summary.
    class_weights : list, optional
        Weights for class balancing.
    log_dir : str, default="logs"
        Directory for logs.
    num_plots : int, default=10
        Number of plots to generate.
    input_length : int, optional
        Input length in samples.
    input_mseed : bool, default=True
        Whether input is MSEED format.
    filename_picks : str, default="picks"
        Base name for picks output file.
    one_single_sampling_rate : float, default=-1
        Override sampling rate (-1 for auto).
    data_dir : str
        Directory for prediction data.
    data_list : str
        Path to data list CSV.
    train_dir : str
        Directory for training data.
    train_list : str
        Path to training list CSV.
    valid_dir : str, optional
        Directory for validation data.
    valid_list : str, optional
        Path to validation list CSV.
    output_dir : str, optional
        Directory for output.
    rm_downloads : bool, default=False
        Whether to remove input MSEED files after processing.
    """

    def __init__(
        self,
        pnet_path,
        model_path,
        mode='pred',
        P_threshold=0.3,
        S_threshold=0.3,
        batch_size=2,
        plot=False,
        save_result=False,
        epochs=100,
        learning_rate=0.01,
        decay_step=-1,
        decay_rate=0.9,
        momentum=0.9,
        filters_root=8,
        depth=5,
        kernel_size=[7, 1],
        pool_size=[4, 1],
        drop_rate=0,
        dilation_rate=[1, 1],
        loss_type="cross_entropy",
        weight_decay=0,
        optimizer="adam",
        summary=True,
        class_weights=[1, 1, 1],
        log_dir="log",
        num_plots=10,
        input_length=None,
        input_mseed=True,
        filename_picks="picks",
        one_single_sampling_rate=-1,
        data_dir="./dataset/waveform_pred/",
        data_list="./dataset/waveform.csv",
        train_dir="./dataset/waveform_train/",
        train_list="./dataset/waveform.csv",
        valid_dir=None,
        valid_list=None,
        output_dir=None,
        rm_downloads=False
    ):
        """Initialize PhaseNetObj with prediction parameters."""
        self.phasenet_path = pnet_path
        self.model_path = model_path
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
        self.overlap = 0.5
        self.rm_downloads = rm_downloads
        self.one_single_sampling_rate = one_single_sampling_rate
        self.name = "PhaseNet"


class EQTransformer:
    """
    Wrapper class for EQTransformer predictions.

    Parameters
    ----------
    eqt_obj : EQTransformerObj
        EQTransformer configuration object.
    """

    def __init__(self, eqt_obj):
        """Initialize EQTransformer wrapper."""
        self.eqt_obj = eqt_obj

    def pick(self, mseed_storage, metadata_dir, out_dir):
        """
        Run EQTransformer phase picking.

        Parameters
        ----------
        mseed_storage : str
            Directory containing MiniSEED files.

        metadata_dir : str
            Directory containing station metadata.

        out_dir : str
            Output directory for prediction results.

        Returns
        -------
        pandas.DataFrame
            Picks converted to SeisMonitor format.

        Notes:
            Saves results to {out_dir}/results/seismonitor_picks.csv
        """
        self.mseed_storage = mseed_storage
        self.json_path = os.path.join(metadata_dir, "stations.json")
        self.pick_storage = os.path.join(out_dir, "results")
        self.msg_author = "Picker: EQTransformer"

        tic = time.time()
        printlog('info', self.msg_author, 'Running EQTransformer')
        
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
                overwrite=self.eqt_obj.overwrite
            )
        except ValueError as e:
            if "NoneType" in str(e) and "length 0" in str(e):
                printlog(
                    'info', self.msg_author,
                    f"Error in batch_size: {self.eqt_obj.batch_size} to the last mseed. "
                    "Check 'to_pick' parameter in DownloadRestrictions"
                )

        output_path = os.path.join(self.pick_storage, "seismonitor_picks.csv")
        result = ut.eqt_picks_2_seismonitor_fmt(
            self.pick_storage, self.mseed_storage, output_path
        )
        
        tf.keras.backend.clear_session()
        toc = time.time()
        exetime = dt.timedelta(seconds=toc - tic)
        printlog(
            'info', self.msg_author,
            f'Total time of execution: {exetime.total_seconds()} seconds'
        )

        if self.eqt_obj.rm_downloads:
            shutil.rmtree(self.mseed_storage, ignore_errors=True)
        
        return result


class PhaseNet:
    """
    Wrapper class for PhaseNet predictions.

    Parameters
    ----------
    pnet_obj : PhaseNetObj
        PhaseNet configuration object.
    """

    def __init__(self, pnet_obj):
        """Initialize PhaseNet with configuration object."""
        self.pnet_obj = pnet_obj
        self.msg_author = "Picker: PhaseNet"

    def __mv_downloads2onefolder(self):
        """Move MSEED files to a single folder for processing."""
        tic = time.time()
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)

        ut.mv_mseed2onefolder(self.mseed_storage, self.datadir)
        toc = time.time()
        exetime = dt.timedelta(seconds=toc - tic)
        printlog(
            "info", self.msg_author,
            f"Moving mseed2onefolder in {exetime.total_seconds()} seconds."
        )

    def __mv_downloads2stationfolder(self):
        """Move MSEED files back to station-specific folders."""
        tic = time.time()
        ut.mv_mseed2stationfolder(self.datadir, self.mseed_storage)
        toc = time.time()
        exetime = dt.timedelta(seconds=toc - tic)
        printlog(
            "info", self.msg_author,
            f"Moving mseed2stationfolder in {exetime.total_seconds()} seconds."
        )

    def __make_datalist(self):
        """Create datalist for PhaseNet processing."""
        tic = time.time()
        printlog("info", self.msg_author, "Running to create datalist")
        
        ut.make_phasenet_datalist(
            datadir=self.datadir,
            json_path=self.json_path,
            datalist_path=self.datalist
        )
        
        toc = time.time()
        exetime = dt.timedelta(seconds=toc - tic)
        printlog(
            "info", self.msg_author,
            f"Datalist created in {exetime.total_seconds()} seconds."
        )

    def pick(self, mseed_storage, metadata_dir, out_dir):
        """Run PhaseNet phase picking.

        Parameters
        ----------
        mseed_storage : str
            Directory containing MiniSEED files.

        metadata_dir : str
            Directory containing station metadata.

        out_dir : str
            Output directory for prediction results.

        Returns
        -------
        pandas.DataFrame
            Picks converted to SeisMonitor format.

        Notes:
            Saves results to {out_dir}/results/seismonitor_picks.csv
        """
        self.mseed_storage = mseed_storage
        self.json_path = os.path.join(metadata_dir, "stations.json")
        self.pick_storage = os.path.join(out_dir, "results")
        self.datadir = os.path.join(out_dir, "datadir")
        self.datalist = os.path.join(out_dir, "datalist", 'fname.csv')

        self.__mv_downloads2onefolder()
        self.__make_datalist()

        tf.compat.v1.reset_default_graph()
        tf.keras.backend.clear_session()

        tic = time.time()
        printlog("info", self.msg_author, "Running PhaseNet pretrained model.")

        args = self.pnet_obj
        args.data_dir = self.datadir
        args.data_list = self.datalist
        args.output_dir = self.pick_storage

        ut.phasenet_from_console(args, self.msg_author)
        datapicks = os.path.join(self.pick_storage, 'picks.csv')
        seismonitor_datapicks = os.path.join(self.pick_storage, 'seismonitor_picks.csv')
        
        result = ut.get_picks(
            datapicks=datapicks,
            datalist=self.datalist,
            min_p_prob=self.pnet_obj.tp_prob,
            min_s_prob=self.pnet_obj.ts_prob,
            one_single_sampling_rate=self.pnet_obj.one_single_sampling_rate,
            mode='df_obj',
            export=seismonitor_datapicks
        )

        try:
            ut.rm_phasenet_duplicate_picks(
                path=seismonitor_datapicks,
                output_path=seismonitor_datapicks
            )
            printlog(
                'info', self.msg_author,
                f'Duplicated picks had been removed. See your results in {seismonitor_datapicks}'
            )
        except Exception as e:
            printlog(
                'error', self.msg_author,
                f"Can't remove duplicated picks: {e}"
            )

        toc = time.time()
        exetime = dt.timedelta(seconds=toc - tic)
        printlog(
            "info", self.msg_author,
            f'Time of prediction: {exetime.total_seconds()} seconds'
        )

        self.__mv_downloads2stationfolder()
        if self.pnet_obj.rm_downloads:
            shutil.rmtree(self.mseed_storage, ignore_errors=True)

        return result