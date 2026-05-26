import logging
from pathlib import Path
from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
from SeisMonitor.monitor.picker import utils as piut

monitor_path = Path(__file__).parent / "sm"
downloads_path = monitor_path / "downloads"
stations_path = monitor_path / "stations"
picks_path = monitor_path / "picks"
sm_picks_path = picks_path / "seismonitor_picks.csv"

#check eqt models here: https://github.com/smousavi05/EQTransformer/tree/master/ModelsAndSampleData
eqt_model = "/groups/igonin/ecastillo/others/picking_models/EQTransformer_models/EqT_model.h5"

####### default configuration
logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] [%(name)s]  %(message)s',
        datefmt='%m-%d %H:%M'
    )

eqtobj = EQTransformerObj(model_path=eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.2,
            P_threshold = 0.1,
            S_threshold = 0.1,
            batch_size = 10, # This has to align with the batch size used in the downloader
            number_of_plots = 10,
            plot_mode = None ) 
eqt = EQTransformer(eqtobj)
eqt.pick(mseed_storage=str(downloads_path),
        metadata_dir= str(stations_path),
        out_dir=str(picks_path))
piut.eqt_picks_2_seismonitor_fmt(eqt_folder=str(picks_path),
                                mseed_folder=str(downloads_path),
                                out_path=str(sm_picks_path))