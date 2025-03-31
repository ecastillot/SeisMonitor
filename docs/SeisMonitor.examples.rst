Examples
===================

WaveformRestrictions
--------------------

`WaveformRestrictions Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/core/objects.html#WaveformRestrictions>`_

.. code:: python
   
   sgc_rest = WaveformRestrictions(network="CM",
                     station="URMC,CLEJA,VILL,PRA,ORTC,GARC,FLO2,CHI,YOT",
                     location="*",
                     channel="*",
                     starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                     endtime=UTCDateTime("2019-12-25T01:00:00.000000Z"),
                     location_preferences=["","00","20","10","40"],
                     channel_preferences=["HH","BH","EH","HN","HL"],
                     filter_networks=[], 
                     filter_stations=[],
                     filter_domain= [-83.101,-64.549,-2.229,14.945],
                     )

Provider
-----------

`Provider Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/core/objects.html#Provider>`_

.. code:: python
  
   sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
   sgc_provider = Provider(sgc_client,sgc_rest)

Downloader
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/1.downloader.ipynb
   :alt: Open In Colab

.. code:: python
   
   mseed_storage = os.path.join(monitor_path,"downloads/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
   md = MseedDownloader(providers=[sgc_provider])
   md.download(mseed_storage,
               picker_args={"batch_size":100,"overlap":0.3,"length":60},
               chunklength_in_sec=7200,n_processor=None)

Picker
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/2.picker.ipynb
   :alt: Open In Colab

.. code:: python
   
   from SeisMonitor.monitor.picker.ai import EQTransformer,EQTransformerObj
   from SeisMonitor.monitor.picker import utils as piut

   # Eqtransformer object
   eqtobj = EQTransformerObj(model_path=eqt_model,
            n_processor = 6,
            overlap = 0.3,
            detection_threshold =0.1,
            P_threshold = 0.01,
            S_threshold = 0.01,
            batch_size = 20,
            number_of_plots = 0,
            plot_mode = 1 ) 

   out_dir = os.path.join(monitor_path ,"picks","eqt")
   result = os.path.join(monitor_path ,"picks","eqt","seismonitor_picks.csv")

   #picker
   eqt = EQTransformer(eqtobj)
   eqt.pick(downloads,stations,out_dir)

   # eqt to seismonitor format
   piut.eqt_picks_2_seismonitor_fmt(out_dir,downloads,result)
      

Associator
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/3.associator.ipynb
   :alt: Open In Colab

.. code:: python
   
   from SeisMonitor.monitor.associator.ai import GaMMA,GaMMAObj
   from SeisMonitor.monitor.associator import utils as asut
   
   # Region   lonw,lone,lats,latn, zmin_km,zmax_km
   region = [-76.729, -72.315,1.55, 5.314,0, 150]

   # gamma object
   gc = GaMMAObj(region,"EPSG:3116",
                  use_amplitude = False,
                  use_dbscan=False,
                  calculate_amp=False)

   inv = os.path.join(stations,"inv.xml")
   picks = os.path.join(picks,"eqt_seismonitor_picks.csv")
   out_dir = os.path.join(monitor_path,"gamma_asso","eqt")

   g = GaMMA(gc)
   obspy_catalog, df_catalog,df_picks = g.associate(picks,inv,out_dir)
   