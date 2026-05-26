.. _installation-section:
Installation
============================

.. image:: https://img.shields.io/pypi/v/SeisMonitor?
   :target: https://pypi.org/project/SeisMonitor/
   :alt: PyPI version

Installing via Conda
----------------------------
You can install `SeisMonitor <https://pypi.org/project/seismonitor/>`_ directly from PyPI using `pip` in a conda environment. 
This will install the Python package and its dependencies.

.. code-block:: bash

   conda create --name seismonitor python=3.10
   conda activate seismonitor
   pip install SeisMonitor
   pip install git+https://github.com/ecastillot/EQTransformer.git@master
   pip install git+https://github.com/wayneweiqiang/GaMMA.git
