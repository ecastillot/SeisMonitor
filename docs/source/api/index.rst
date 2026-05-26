.. _api-section:
API Reference
=============

SeisMonitor provides tools for seismic monitoring workflows, including
waveform handling, downloading, phase picking, association, event
location, and magnitude estimation.

.. grid:: 2

   .. grid-item-card:: Core
      :link: SeisMonitor.core.objects
      :link-type: doc

      Main objects.

   .. grid-item-card:: Downloader
      :link: SeisMonitor.monitor.downloader
      :link-type: doc

      Download waveform data from providers.

   .. grid-item-card:: Picker
      :link: SeisMonitor.monitor.picker
      :link-type: doc

      Automatic phase picking workflows.

   .. grid-item-card:: Associator
      :link: SeisMonitor.monitor.associator
      :link-type: doc

      Earthquake phase association tools.

   .. grid-item-card:: NonLinLoc Locator
      :link: SeisMonitor.monitor.locator.nlloc
      :link-type: doc

      Earthquake location tools using NonLinLoc.

   .. grid-item-card:: Magnitude
      :link: SeisMonitor.monitor.magnitude
      :link-type: doc

      Magnitude estimation methods.

.. toctree::
   :maxdepth: 3

   SeisMonitor.core.objects
   SeisMonitor.monitor.downloader
   SeisMonitor.monitor.picker
   SeisMonitor.monitor.associator
   SeisMonitor.monitor.locator.nlloc
   SeisMonitor.monitor.magnitude
