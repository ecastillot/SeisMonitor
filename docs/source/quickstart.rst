.. _quickstart-section:
Quickstart
============

This guide provides a **basic overview of SeisMonitor** and demonstrates its application using seismic data from Texas.

1. :ref:`download-section`
2. :ref:`detection-section`
3. :ref:`asso-section`
4. :ref:`location-section`
5. :ref:`magnitude-section`

.. warning::

   Make sure to use the latest stable version of SeisMonitor.:
   .. image:: https://img.shields.io/pypi/v/SeisMonitor?style=plastic
      :target: https://pypi.org/project/SeisMonitor/
      :alt: PyPI version

.. _download-section:
Download Earthquake Data
-------

:red:`Key Classes:` :class:`~SeisMonitor.core.objects.WaveformRestrictions`, :class:`~SeisMonitor.core.objects.Provider`, :class:`~SeisMonitor.monitor.downloader.seismonitor.MseedDownloader` 

:red:`Key Funcs:` :func:`~SeisMonitor.monitor.downloader.seismonitor.MseedDownloader.make_inv_and_json`, :func:`~SeisMonitor.monitor.downloader.seismonitor.MseedDownloader.download`

.. code-block:: python

   from pathlib import Path
   from obspy.core.utcdatetime import UTCDateTime
   from obspy.clients.fdsn import Client as FDSNClient
   from SeisMonitor.core.objects import WaveformRestrictions,Provider
   from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader


   monitor_path = Path(__file__).parent / "sm"


   downloads_path = monitor_path / "downloads"
   stations_path = monitor_path / "stations"


   client = FDSNClient("http://rtserve.beg.utexas.edu")
   chunklength_in_sec = 3600
   rest = WaveformRestrictions(network="*",
                     station="PB13,PB17,PB20,PB23,PB34,PB58,PB40,PB24,PB25,PB26,WB03",
                     location="*",
                     channel="*",
                     starttime=UTCDateTime("2022-11-16T21:00:00.000000Z"),
                     endtime=UTCDateTime("2022-11-16T23:00:00.000000Z"),
                     location_preferences=["","00","20","10","40"],
                     channel_preferences=["HH","BH","EH","HN","HL"],
                     filter_networks=[], 
                     filter_stations=[],
                     )


   ####### Default configuration
   provider = Provider(client,rest)
   md = MseedDownloader(providers=[provider])
   inv,json = md.make_inv_and_json(stations_path)

   mseed_storage = downloads_path / "{station}" / "{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
   md.download(str(mseed_storage),
               picker_args={"batch_size":10,
                           "overlap":0.3,
                           "length":60},
               chunklength_in_sec=chunklength_in_sec,
               n_processor=1)

.. _detection-section:
Earthquake Detection & Phase Picking
-----------------------------------

:red:`Key Classes:` :class:`~SeisMonitor.monitor.picker.ai.EQTransformerObj`, :class:`~SeisMonitor.monitor.picker.ai.EQTransformer`

:red:`Key Funcs:` :func:`~SeisMonitor.monitor.picker.ai.EQTransformer.pick`, :func:`~SeisMonitor.monitor.picker.utils.eqt_picks_2_seismonitor_fmt`

.. code-block:: python

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

.. _asso-section:
Phase Association
---------------------

:red:`Key Classes:` :class:`~SeisMonitor.monitor.associator.ai.GaMMAObj`, :class:`~SeisMonitor.monitor.associator.ai.GaMMA`

:red:`Key Funcs:` :func:`~SeisMonitor.monitor.associator.ai.GaMMA.associate`

.. code-block:: python

   import os
   from pathlib import Path
   from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
   from SeisMonitor.monitor.associator import utils as asut
   import matplotlib.pyplot as plt

   monitor_path = Path(__file__).parent / "sm"
   stations_path = monitor_path / "stations"
   picks_path = monitor_path / "picks"
   asso_path = monitor_path / "associations"

   # region = [min_lon, max_lon, min_lat, max_lat, min_depth_km, max_depth_km]
   # EPSG:3081 -> UTM zone 14N, suitable for Texas
   region = [-104.17816, -103.80355, 31.48921, 31.72133,0, 12] 
   epsg = "EPSG:3081" 


   ####### default configuration
   sm_picks_path = picks_path / "seismonitor_picks.csv"
   inv_stations_path = monitor_path / "stations" / "inv.xml"

   gc = GaMMAObj(region,epsg,
                  use_amplitude = False,
                  use_dbscan=False,
                  calculate_amp=False,
                  method="BGMM",
                  min_picks_per_eq=5,
                  oversample_factor=1,
                  max_sigma11=2.0,
                  vel = {"p": 6, "s": 6/ 1.70})


   g = GaMMA(gc)
   obspy_catalog, df_catalog,df_picks = g.associate(picks_csv=sm_picks_path,
                                       xml_path=inv_stations_path,
                                       out_dir=asso_path)
   print("Catalog\n",obspy_catalog)
   print("Events\n",obspy_catalog.to_df())
   print("Picks\n",obspy_catalog.picks_to_df())

.. _location-section:
Earthquake Location
-------

.. warning::

   The NonLinLoc (NLLoc) workflow in SeisMonitor is currently supported
   only on **Ubuntu Linux systems**. Other operating systems (Windows/macOS)
   are not officially tested and may fail due to system-level dependencies.

:red:`Key Classes:` :class:`~SeisMonitor.monitor.locator.nlloc.nlloc.NLLoc`, :class:`~SeisMonitor.monitor.locator.utils.VelModel`, :class:`~SeisMonitor.monitor.locator.utils.Stations`

:red:`Key Funcs:` :func:`~SeisMonitor.monitor.locator.nlloc.nlloc.NLLoc.compute_travel_times`, :func:`~SeisMonitor.monitor.locator.nlloc.nlloc.NLLoc.locate`

.. code-block:: python

   import os
   from pathlib import Path
   from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
   from SeisMonitor.monitor.locator import utils as lut

   main_nlloc_path = "/opt/ohpc/pub/apps/nonlinloc"
   vel_path = "/groups/igonin/ecastillo/SeisMonitor/examples/TX/sm/vel_model/DB_model.csv"

   #make sure the region covers all the stations and expected earthquake locations.  
   # Specially the elevation part and the stations, 0 respect to sea level, 
   # negative means below sea level, positive means above sea level.
   region = [-104.17816, -103.80355, 31.48921, 31.72133,-2, 12] 
   delta_in_km = 0.5 # grid spacing for nlloc, smaller means more precise but also more computationally expensive.

   monitor_path = Path(__file__).parent / "sm"
   stations_path = monitor_path / "stations"
   picks_path = monitor_path / "picks"
   asso_path = monitor_path / "associations"
   loc_path = monitor_path / "locations"

   tt_loc_folder = monitor_path / "tt" #travel time folder for nlloc, can be anywhere but it could consume significant disk space.

   ####### default configuration

   nlloc_output_path = loc_path / "nlloc"
   xml_stations_path = monitor_path / "stations" / "inv.xml"
   xml_asso_path = asso_path / "associations.xml"
   xml_nlloc_path = nlloc_output_path / "nlloc_catalog.xml"

   vel_model = lut.VelModel(str(vel_path),model_name="Velmodel1D")
   stations = lut.Stations(str(xml_stations_path))

   nlloc = NLLoc(
         core_path = str(main_nlloc_path), ### type your NLLoc path, 
         agency="SeisMonitor",
         region = region,
         vel_model = vel_model,
         stations = stations,
         delta_in_km = delta_in_km,
         tmp_folder=str(tt_loc_folder) 
         )

   # Use this to compute travel times only once, and then you can reuse the travel time files for future locations.
   nlloc.compute_travel_times()
   print(nlloc.tmp_folder)
   print("nlloc.tmp_folder dir: ",nlloc.tmp_folder)
   print("Important folders in nlloc.tmp_folder directory",os.listdir(nlloc.tmp_folder))

   eqt_nlloc_catalog = nlloc.locate(catalog=str(xml_asso_path),
                              nlloc_out_folder= str(nlloc_output_path),
                              out_filename = str(xml_nlloc_path.name),
                              out_format="QUAKEML" )
   print("Catalog\n",eqt_nlloc_catalog )
   print("Events\n",eqt_nlloc_catalog.to_df())
   print("Picks\n",eqt_nlloc_catalog.picks_to_df())

.. _magnitude-section:
Local Magnitude
-------

:red:`Key Classes:` :class:`~SeisMonitor.core.objects.WaveformRestrictions`, :class:`~SeisMonitor.core.objects.Provider`, :class:`~SeisMonitor.monitor.magnitude.mag.Magnitude`

:red:`Key Funcs:` :func:`~SeisMonitor.monitor.magnitude.mag.Magnitude.get_Ml`

.. code-block:: python

   import math
   from pathlib import Path
   from obspy.clients.fdsn import Client as FDSNClient
   from obspy.core.utcdatetime import UTCDateTime
   from SeisMonitor.monitor.magnitude.mag import Magnitude
   from SeisMonitor.core.objects import WaveformRestrictions,Provider

   monitor_path = Path(__file__).parent / "sm"

   loc_path = monitor_path / "locations"
   mag_path = monitor_path / "magnitude"

   nlloc_output_path = loc_path / "nlloc"
   xml_nlloc_path = nlloc_output_path / "nlloc_catalog.xml"

   client = FDSNClient("http://rtserve.beg.utexas.edu")
   rest = WaveformRestrictions(network="*",
                     station="PB13,PB17,PB20,PB23,PB34,PB58,PB40,PB24,PB25,PB26,WB03",
                     location="*",
                     channel="*",
                     starttime=UTCDateTime("2022-11-16T21:00:00.000000Z"),
                     endtime=UTCDateTime("2022-11-16T23:00:00.000000Z"),
                     location_preferences=["","00","20","10","40"],
                     channel_preferences=["HH","BH","EH","HN","HL"],
                     filter_networks=[], 
                     filter_stations=[],
                     )


   provider = Provider(client,rest)
   mag = Magnitude([provider],xml_nlloc_path,mag_path) #catalog,providers,out

   # Use your own Ml formula, here is an example for local magnitude, but make sure to check if the built-in formula is suitable for your region and data.
   # For Texas using Kavoura et al (2020)
   Ml = lambda ampl,epi_dist : math.log10(ampl*1e3 ) + 1.54 * math.log10(epi_dist) + 0.0002*(epi_dist-100) - 0.08

   cat = mag.get_Ml(mag_type=Ml ,
               trimmedtime=5, #seconds after pick S to trim the signal
               out_format="SC3ML")
   print(cat)
