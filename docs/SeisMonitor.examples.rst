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

`EQTransformerObj Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/picker/ai.html#EQTransformerObj>`_

`EQTransformer.pick Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/picker/ai.html#EQTransformer>`_

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

`GaMMAObj Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/associator/ai.html#GaMMAObj>`_

`GaMMA.associate Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/associator/ai.html#GaMMA>`_

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

Locator
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/4.locator.ipynb
   :alt: Open In Colab

`NLLoc Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/locator/nlloc/nlloc.html#NLLoc>`_

`VelModel Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/locator/utils.html#VelModel>`_

`Stations Documentation <https://seismonitor.readthedocs.io/en/latest/_modules/SeisMonitor/monitor/locator/utils.html#Stations>`_

Stations and velocity input data: `here <https://github.com/ecastillot/SeisMonitor/tree/others>`_ 

.. code:: python
   
   from SeisMonitor.monitor.locator.nlloc.nlloc import NLLoc
   from SeisMonitor.monitor.locator import utils as lut
   
   # vel model and stations
   vel_path = os.path.join(velmodel,"vel1d_col.csv")
   inv = os.path.join(stations,"inv.xml")

   vel_model = lut.VelModel(vel_path,model_name="Ojeda&Havskov(2004)")
   stations = lut.Stations(inv)

   # nlloc definition
   nlloc = NLLoc(
        core_path = nlloc_path, ### type your NLLoc path, 
        agency="SeisMonitor",
        region = [-85, -68,0, 15,-5, 205], #lonw,lone,#lats,latn,zmin_km,zmax_km
        vel_model = vel_model,
        stations = stations,
        delta_in_km = 2.5,
        tmp_folder=os.path.join(os.getcwd(),"NLLoc_grid") ### CHANGE PATH TO YOUR OWN PATH AND ALSO TAKE IN MIND THAT CONSUME DISK
        )

   # travel times (time consuming)
   nlloc.compute_travel_times()

   #earthquake location
   eqt_catalog = os.path.join(asso,"associations.xml")
   eqt_nlloc_catalog = nlloc.locate(catalog=eqt_catalog,
                              nlloc_out_folder= out_dir,
                              out_filename = "LOC.xml",
                              out_format="SC3ML" )

Locator
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/5.magnitude.ipynb
   :alt: Open In Colab

.. code:: python
   mag = Magnitude([sgc_provider],catalog,out_dir) #catalog,providers,out

   # parameter definitions for local magnitude
   ml_params = {"a":1.019,"b":0.0016,"r_ref":140} #ojeda
   k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
                     ml_params["b"]* (ml_params["r_ref"]-100) +3
   Ml = lambda ampl,epi_dist : math.log10(ampl * 1e3) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
                                 ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k

   cat = mag.get_Ml(mag_type=Ml ,
            trimmedtime=5, #seconds after pick S to trim the signal
            out_format="SC3ML")

All in one
-----------

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/monitor.ipynb
   :alt: Open In Colab