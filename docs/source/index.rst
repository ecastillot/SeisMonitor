.. SeisMonitor documentation master file, created by
   sphinx-quickstart on Fri Jan 30 13:39:44 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <h1>
   Welcome to 
   <span style="color:#335BFF;">S</span>eis<span style="color:#335BFF;">M</span>onitor
   </h1>

.. image:: https://img.shields.io/badge/GitHub-SeisMonitor-black?style=for-the-badge&logo=github
   :target: https://github.com/ecastillot/SeisMonitor/tree/master
   :alt: GitHub Repository
.. image:: https://img.shields.io/pypi/v/SeisMonitor?label=pypi
   :target: https://pypi.org/project/SeisMonitor/
   :alt: PyPI version

.. raw:: html

   <br><br>

**A Python open-source package for seismological monitoring**

Workflow
--------------------------

SeisMonitor monitors seismic activity that uses ready-made ML methods for event detection, phase picking and association, and other well-known methods for the rest of the steps

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/SeisMonitor_flow.png
   :width: 800px


The framework can be applied to any region worldwide. We applied it to Colombia using nearly 7 years (2016–2022) of seismic data, successfully revealing significant tectonic structures and crustal faults.



Example for Colombia
--------------------------
**Paper**: `Colombian Seismic Monitoring Using Advanced Machine-Learning Algorithms <https://gaprieto.com/wp-content/uploads/2024/05/gprieto_24c.pdf>`_ 

Emmanuel Castillo, Daniel Siervo, Germán A. Prieto; Colombian Seismic Monitoring Using Advanced Machine‐Learning Algorithms. Seismological Research Letters 2024;; 95 (5): 2971–2985. doi: `https://doi.org/10.1785/0220240036 <https://doi.org/10.1785/0220240036>`_ 

**Dataset**: :ref:`dataset-section`

.. image:: https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/sm_overview.png
   :alt: NLLoc Catalog
   :width: 700px
   :align: center


Content
--------------------------

.. toctree::
   :maxdepth: 2
   
   home
   installation
   quickstart
   dataset
   api/index
   authors
   contribution
