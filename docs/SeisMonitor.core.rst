Core
========================

client
------------------------------

.. automodule:: SeisMonitor.core.client
   :members:
   :undoc-members:
   :show-inheritance:

This script provides an example of creating a client class for a specific data structure archive on a local filesystem.

The data structure is organized as follows:

``{root_path}/{field_name}/seedfiles/{year}-{month:02d}/{year}-{month:02d}-{day:02d}/{network}.{station}.{location}.{channel}.{year}.{julianday:03d}``

To achieve this, I built the ``LocalClient`` class as a subclass of the SDS client. It inherits the SDS client functionalities. However, ``field_name``, ``month``, and ``day`` are not used in the format string (``fmt``) of the base SDS client instance. Therefore, it’s mandatory to override two private methods: ``_get_filenames`` and ``_get_filename``.

The mandatory parameters for the ``LocalClient`` class are ``root_path`` and ``field_name``.

Example
-------
Here’s how to use the ``LocalClient`` class::

   root_path = "/home/emmanuel/myarchive"
   client = LocalClient(root_path, "FIELD_1")
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

