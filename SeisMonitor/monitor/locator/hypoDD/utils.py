import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from datetime import datetime
import numpy as np
import pandas as pd
from tqdm import tqdm
from glob import glob
from multiprocessing import Pool
import time
import SeisMonitor.utils as ut

"""
Some utils functions are taken from: https://github.com/wayneweiqiang/QuakeFlow/blob/master/HypoDD/gamma2hypodd.py
such us:
    - download_hypodd()
"""

CORE_HYPODD = os.path.join(os.path.dirname(__file__),"core")

def download_hypodd():
    '''
    HypoDD can be downloaded from https://www.ldeo.columbia.edu/~felixw/hypoDD.html
    Helpful compiling flags: FFLAGS = -O -I${INCLDIR} -mcmodel=large
    '''
    name = "HYPODD_1.3.tar.gz"
    gz_path = os.path.join(CORE_HYPODD,name)
    hypoDD_path = os.path.join(CORE_HYPODD,"HYPODD")
    make_path = os.path.join(hypoDD_path,"src")
    f77_path = os.path.join(CORE_HYPODD,"f77")
    g77_path = os.path.join(CORE_HYPODD,"g77")

    isfile = ut.isfile(gz_path)
    if not isfile:
        os.system(f"wget -O {gz_path} http://www.ldeo.columbia.edu/~felixw/HYPODD/HYPODD_1.3.tar.gz")
    
    if not os.path.isdir(hypoDD_path):
        os.system(f"tar -xf {gz_path} -C {CORE_HYPODD}")

    isfile = ut.isfile(f77_path)
    if not isfile:
        os.system(f"ln -s $(which gfortran) {f77_path}")

    isfile = ut.isfile(g77_path)
    if not isfile:
        os.system(f"ln -s $(which gfortran) {g77_path}")

    os.environ['PATH'] += os.pathsep + CORE_HYPODD


    # print(os.environ['PATH'])
    os.system(f"make -C {make_path}")

if __name__ == "__main__":
    download_hypodd()