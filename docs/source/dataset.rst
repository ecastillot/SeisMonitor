.. _dataset-section:
SeisMonitor in Colombia 
============================

We applied it to Colombia using nearly 7 years (2016–2022) of seismic data, successfully revealing significant tectonic structures and crustal faults.


Paper
-------
See it here: `Colombian Seismic Monitoring Using Advanced Machine-Learning Algorithms <https://gaprieto.com/wp-content/uploads/2024/05/gprieto_24c.pdf>`_ 

Emmanuel Castillo, Daniel Siervo, Germán A. Prieto; Colombian Seismic Monitoring Using Advanced Machine‐Learning Algorithms. Seismological Research Letters 2024;; 95 (5): 2971–2985. doi: `https://doi.org/10.1785/0220240036 <https://doi.org/10.1785/0220240036>`_ 


Colombia Dataset
-------

.. list-table::
   :header-rows: 1
   :widths: 40 20

   * - Data
     - Link
   * - EQTransformer-picks
     - `Open <https://drive.google.com/file/d/1e3044OJBtFjg4HrrbawEJ-3GfqjN-Cbf/view?usp=sharing>`_
   * - GaMMA-catalog
     - `Open <https://drive.google.com/file/d/1OnMDNe4NZK98mNLjdrFhZNy14dJHZXhX/view?usp=drive_link>`_
   * - GaMMA-picks
     - `Open <https://drive.google.com/file/d/1qJhRHIYpFh8_TyWV0KyOvq9nj4y_zz5T/view?usp=drive_link>`_
   * - SeisMonitor-catalog
     - `Open <https://drive.google.com/file/d/1ZphhiOkZkeOZBwBVPmXmVsXp7WL0ZjEL/view?usp=sharing>`_
   * - SeisMonitor-picks
     - `Open <https://drive.google.com/file/d/1SS1OysOSk-9l-gFRpuCORRgXvwTKeniG/view?usp=sharing>`_

SeisMonitor Workflow
-------

Seismic Data
^^^^^^^^^^^^^^^^^^^^^
More than 100 stations across Colombia Seismic Network (`CM <https://www.fdsn.org/networks/detail/CM/>`_ )

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/MAP.png
   :alt: Colombia Seismic Network
   :width: 300px
   :align: center

Earthquake Detection & Phase Picking
^^^^^^^^^^^^^^^^^^^^^
We used `EQTransformer <https://eqtransformer.readthedocs.io/en/latest/>`_ to detect earthquakes and pick phases.

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/mesetas_examples.png
   :alt: Mesetas Earthquake Detection
   :width: 500px
   :align: center

Phase Association
^^^^^^^^^^^^^^^^^^^^^
We used `GaMMA <https://github.com/AI4EPS/GaMMA>`_ to associate phases and obtained our first catalog.

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/gamma_overview.PNG
   :alt: GaMMA Catalog
   :width: 700px
   :align: center


NLloc Catalog
^^^^^^^^^^^^^^^^^^^^^
Although preliminary GaMMA locations are good enough to allow us to visualize the main seismic activity, we use `NLLoc <http://alomax.free.fr/nlloc/>`_ to further improve earthquake locations 

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/sm_overview.png
   :alt: NLLoc Catalog
   :width: 700px
   :align: center


Colombian Seismicity
^^^^^^^^^^^^^^^^^^^^^
It stands out for its computational efficiency, as well as its earthquake detection and location performance, demonstrating a high level of quality that effectively highlights key tectonic features and seismicity trends in northern South America.

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/map_with_profiles.PNG
   :alt: Colombian Seismicity
   :width: 700px
   :align: center

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/profiles.PNG
   :alt: Cross plots of seismicity
   :width: 700px
   :align: center   