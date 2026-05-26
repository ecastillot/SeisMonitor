![seismonitor](docs/source/_static/logo/seismonitor.PNG)   
![PyPI](https://img.shields.io/pypi/v/SeisMonitor?style=plastic)
![Twitter Follow](https://img.shields.io/twitter/follow/manuavid?style=social)![GitHub followers](https://img.shields.io/github/followers/ecastillot?style=social)![GitHub stars](https://img.shields.io/github/stars/ecastillot/SeisMonitor?style=social)![GitHub forks](https://img.shields.io/github/forks/ecastillot/SeisMonitor?style=social)  

Python package for seismological monitoring. 

SeisMonitor monitors seismic activity that uses ready-made ML methods for event detection, phase picking and association, and other well-known methods for the rest of the steps.

![SeisMonitor Workflow](https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/SeisMonitor_flow.png)


# Documentation

 [![Read the Docs](https://img.shields.io/badge/Docs-Read%20the%20Docs-blue?style=for-the-badge&logo=read-the-docs)](https://seismonitor.readthedocs.io/en/latest/)


# Installation

## pip
```bash
pip install SeisMonitor
pip install git+https://github.com/ecastillot/EQTransformer.git@master
pip install git+https://github.com/wayneweiqiang/GaMMA.git
```
## conda
```bash
conda create --name seismonitor python=3.10
conda activate seismonitor
pip install SeisMonitor
pip install git+https://github.com/ecastillot/EQTransformer.git@master
pip install git+https://github.com/wayneweiqiang/GaMMA.git
```

# Colombian Catalog
![Colombia Catalog](https://raw.githubusercontent.com/ecastillot/SeisMonitor/master/figures/sm_overview.png)

Castillo, E., Prieto, G.A., Siervo, D. (2024). Colombian seismic monitoring using advanced machine-learning algorithms. Seism. Res. Lett. [PDF](https://gaprieto.com/wp-content/uploads/2024/05/gprieto_24c.pdf) 


| Data |  |
|---|---|
| EQTransformer-picks| [Open](https://drive.google.com/file/d/1e3044OJBtFjg4HrrbawEJ-3GfqjN-Cbf/view?usp=sharing) |
| GaMMA-catalog| [Open](https://drive.google.com/file/d/1OnMDNe4NZK98mNLjdrFhZNy14dJHZXhX/view?usp=drive_link) |
| GaMMA-picks| [Open](https://drive.google.com/file/d/1qJhRHIYpFh8_TyWV0KyOvq9nj4y_zz5T/view?usp=drive_link) |
| SeisMonitor-catalog| [Open](https://drive.google.com/file/d/1ZphhiOkZkeOZBwBVPmXmVsXp7WL0ZjEL/view?usp=sharing) |
| SeisMonitor-picks| [Open](https://drive.google.com/file/d/1SS1OysOSk-9l-gFRpuCORRgXvwTKeniG/view?usp=sharing) | 


# Authors

<table>
  <tr>
    <td valign="middle">
      <strong>Emmanuel Castillo</strong><br>
      University of Texas at Dallas<br>
      <a href="https://www.linkedin.com/in/ecastillot/" target="_blank">
        <img src="https://img.shields.io/badge/LinkedIn-Follow-blue?logo=linkedin&logoColor=white&style=flat" height="20">
      </a>
      <a href="mailto:castillo.280997@gmail.com" target="_blank">
        <img src="https://img.shields.io/badge/Email-Send-red?logo=gmail&logoColor=white&style=flat" height="20">
      </a>
    </td>
    <td width="100" align="center">
      <img src="docs/source/_static/photos/emmanuel.jpg" alt="Emmanuel" width="80" height="80" style="border-radius:50%; object-fit:cover;">
    </td>
  </tr>
  <tr>
    <td valign="middle">
      <strong>Daniel Siervo</strong><br>
      University of Texas at Dallas<br>
      <a href="https://www.linkedin.com/in/dsiervo/" target="_blank">
        <img src="https://img.shields.io/badge/LinkedIn-Follow-blue?logo=linkedin&logoColor=white&style=flat" height="20">
      </a>
      <a href="mailto:daniel.siervo@beg.utexas.edu" target="_blank">
        <img src="https://img.shields.io/badge/Email-Send-red?logo=gmail&logoColor=white&style=flat" height="20">
      </a>
    </td>
    <td width="100" align="center">
      <img src="docs/source/_static/photos/dsiervo.jpg" alt="Nadine" width="80" height="80" style="border-radius:50%; object-fit:cover;">
    </td>
  </tr>
  <tr>
    <td valign="middle">
      <strong>German Prieto</strong><br>
      Universidad Nacional de Colombia<br>
      <a href="https://www.linkedin.com/in/german-prieto-2798495/" target="_blank">
        <img src="https://img.shields.io/badge/LinkedIn-Follow-blue?logo=linkedin&logoColor=white&style=flat" height="20">
      </a>
      <a href="gaprietogo@unal.edu.co" target="_blank">
        <img src="https://img.shields.io/badge/Email-Send-red?logo=gmail&logoColor=white&style=flat" height="20">
      </a>
    </td>
    <td width="100" align="center">
      <img src="docs/source/_static/photos/gaprieto.jpg" alt="Marine" width="80" height="80" style="border-radius:50%; object-fit:cover;">






<!-- | Examples |  |
|---|---|
| 1.downloader| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/1.downloader.ipynb) |
| 2.picker| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/2.picker.ipynb) |
| 3.associator| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/3.associator.ipynb) |
| 4.locator| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ecastillot/SeisMonitor/blob/master/examples/4.locator.ipynb) | -->
