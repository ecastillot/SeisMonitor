Core
========================

client
------------------------------

.. automodule:: SeisMonitor.core.client
   :members:
   :undoc-members:
   :show-inheritance:



Example
-------

The data structure is organized as follows:

``{root_path}/{field_name}/seedfiles/{year}-{month:02d}/{year}-{month:02d}-{day:02d}/{network}.{station}.{location}.{channel}.{year}.{julianday:03d}``

``LocalClient`` is designed to upload local data. It inherits the SDS client functionalities.

.. code:: python

   root_path = "/home/emmanuel/myarchive"
   client = LocalClient(root_path, "/seedfiles/{year}-{month:02d}/{year}-{month:02d}-{day:02d}/{network}.{station}.{location}.{channel}.{year}.{julianday:03d}")
   st = client.get_waveforms("YY", "XXXX", "00",
                           channel="HHZ", starttime=UTCDateTime("20220102T000100"),
                           endtime=UTCDateTime("20220102T000200"))

objects
-------------------------------

.. automodule:: SeisMonitor.core.objects
   :members:
   :undoc-members:
   :show-inheritance:

utils
-----------------------------

.. automodule:: SeisMonitor.core.utils
   :members:
   :undoc-members:
   :show-inheritance:

