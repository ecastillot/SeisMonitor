# SeisMonitor

Set of routines to monitor seismic activity.

*We are still preparing the documentation for you.*

## Requirements

1. Picker: Clone the next repositories at home directory.
         
   - [EQTransformer](https://github.com/smousavi05/EQTransformer.git):
  
   ``` git clone https://github.com/smousavi05/EQTransformer.git``` 

   - [PhaseNet](https://github.com/wayneweiqiang/PhaseNet.git): 
  
   SeisMonitor only supports PhaseNet at specific commit. Then you could use the forked repository by me at this commit.

   ``` git clone https://github.com/ecastillot/PhaseNet.git```

   Otherwise, you could use the next code 

   ``` 
   git clone https://github.com/wayneweiqiang/PhaseNet.git
   git checkout f58d151e2fee0cbb78442e2b4f8072a9ab2095bd
   ```

2. Locator : (It is not ready to publish) 
   - Seisan software.

## Installation

### Conda 
```bash
conda env create -f seismonitor.yml
conda activate seismonitor
```
### Pip

```bash
pip install -r requirements.txt
```
For specific seismonitor environment
```bash
python3.7 -m pip install virtualenv
python3.7 -m virtualenv .seismonitor
source .seismonitor/bin/activate
pip install -r requirements.txt
```


## Examples
There are some examples to use SeisMonitor.

## Automatic Colombian Catalogs

### CM - Colombian Seismic Network
[Scanloc](https://docs.gempa.de/scanloc/current/apps/scanloc.html) is the method used to associate the picks and locate the earthquakes. [Scmag](https://docs.gempa.de/seiscomp3/current/apps/scmag.html) is used to compute the magnitudes. We get two automatic catalogs by each SDL pretrained-model.

- EQTransformer: [CM_eqt_catalog](https://drive.google.com/file/d/1EEBYM6ECLMApnrzeDAQzgDwrWGdeyt4s/view?usp=sharing)
- PhaseNet: [CM_pnet_catalog](https://drive.google.com/file/d/1ZgV0dbv9klnLnLeDksO30EdQp_iOz8PJ/view?usp=sharing)

### VMM - Middle Magdalena valley array
[Scanloc](https://docs.gempa.de/scanloc/current/apps/scanloc.html) is the method used to associate the picks and locate the earthquakes. [Scmag](https://docs.gempa.de/seiscomp3/current/apps/scmag.html) is used to compute the magnitudes. We only get one automatic catalog.

- EQTransformer: [VMM_eqt_catalog](https://drive.google.com/file/d/1BM32a4KUJ2DC-lOu0wdrtrq8gUx8NAUf/view?usp=sharing)

### YU - The Caribbean- Mérida Andes seismic array
[EQTransformer-associator](https://eqtransformer.readthedocs.io/en/latest/tutorial.html#phase-association) is the method used to associate the picks. [SEISAN](http://www.seisan.info/) is used to locate the earthquakes. We only get one automatic catalog.
