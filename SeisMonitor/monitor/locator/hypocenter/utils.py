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
from obspy import read_inventory
import time
from obspy.core.event.catalog import Catalog, read_events
import pexpect
import subprocess
from SeisMonitor.utils import printlog, isfile
# import SeisMonitor.utils as ut

# if 'google.colab' in sys.modules:
#     __file__ =

libgfortran_path = os.path.join(os.path.dirname(__file__),"libgfortran.so.3.0.0")
CORE_SEISAN = os.path.join(os.path.dirname(__file__),"core")
SEISAN_path = os.path.join(CORE_SEISAN,"seismo")
COM_path = os.path.join(SEISAN_path,"COM")
PRO_path = os.path.join(SEISAN_path,"PRO")
DAT_path = os.path.join(SEISAN_path,"DAT")
STATION0_path = os.path.join(DAT_path,"STATION0.HYP")


def download_seisan(libgfortran_path):
    '''
    Seisan can be downloaded from https://www.geo.uib.no/seismo/SOFTWARE/SEISAN/seisan_v12.0_linux_64.tar.gz
    '''
    name = "seisan_v12.0_linux_64.tar.gz"
    gz_path = os.path.join(CORE_SEISAN,name)

    if not os.path.isdir(CORE_SEISAN):
      os.makedirs(CORE_SEISAN)
    if not os.path.isdir(gz_path):
      os.system(f"wget -O {gz_path} https://www.geo.uib.no/seismo/SOFTWARE/SEISAN/seisan_v12.0_linux_64.tar.gz")
    
    os.system("sudo apt-get update")
    os.system("sudo apt-get install gcc")
    os.system("sudo apt-get install g++")
    os.system("sudo apt-get install gfortran")

    if not os.path.isdir(SEISAN_path):
        os.makedirs(SEISAN_path)
        os.system(f"tar -xf {gz_path} -C {SEISAN_path}")
    
    os.system(f"cp {libgfortran_path} /usr/lib/x86_64-linux-gnu/libgfortran.so.3")

def split(sfile_folder,sfilename):
    """
    Parameters:
    -----------
    sfile_folder: str
        It is the folder where is located the sfile.
        In addition, it is the folder where will be saved all sfiles.

    sfilename: str
        Name of the sfile what will be splitted

    Returns:
    --------
        Each sfile produced by the split function

    """
    split_path = os.path.join(PRO_path,"split")
    Split = pexpect.spawn(split_path, cwd=sfile_folder)
    # Split = pexpect.spawn('split', cwd=sfile_folder)
    Split.expect(b' NAME')
    Split.sendline(f'{sfilename}\n')
    Split.expect(b'LOCAL DIRECTORY')
    Split.sendline('\n')
    Split.expect(b'CHARS')
    Split.sendline('usr\n')
    # print(Split.before)
    Split.interact()

def collect(sfile_folder):
    """
    Parameters:
    -----------
    sfile_folder: str
        It is the folder where is located all sfiles.
        In addition, it is the folder where will be saved 
        the collected sfile named as collect.out
    
    Returns:
        colelcted sfile in the next path
        {sfile_folder}/collect.out

    
    """
    collect_path = os.path.join(PRO_path,"collect")
    Collect = pexpect.spawn(collect_path, cwd=sfile_folder)
    Collect.expect(b' return for default')
    Collect.sendline(',,\n')
    Collect.expect(b":")
    Collect.sendline('\n')
    Collect.expect(b":")
    Collect.sendline('\n')
    Collect.expect(b"(Y/N=default)")
    # print(Collect.before)
    Collect.interact()

def update(sfile_folder):
    """
    Parameters:
    -----------
    sfile_folder: str
        It is the folder where is located all sfiles.
            
    Returns:
        It updates each sfile with the hypocenter location
    """
    update_path = os.path.join(PRO_path,"update")
    Update = pexpect.spawn(update_path, cwd=sfile_folder)
    Update.expect(b' return for default')
    Update.sendline(',,\n')
    Update.expect(b":")
    Update.sendline('usr\n')
    Update.expect(b":")
    Update.sendline('\n')
    # print(Update.before)
    Update.interact()

def select(sfile_folder):
    """
    Parameters:
    -----------
    sfile_folder: str
        It is the folder where is located all sfiles.
            
    Returns: 
        It executes select command from SEISAN.
        select.in: index
        select.out: Selected information
    """
    select_path = os.path.join(PRO_path,"select")
    Select = pexpect.spawn(select_path, cwd=sfile_folder)
    Select.expect('FILENAME FOR ONE FILE, MUST BE 6 OR MORE CHARACTERS OR HAVE A .')
    Select.sendline(',,\n')
    Select.expect(' yyyymmddhhmmss:')
    Select.sendline('\n')
    Select.expect(' :')
    Select.sendline('\n')
    Select.expect('RETURN TO SEARCH:')
    Select.sendline('\n')
    # print(Select.before)
    Select.interact()

    out = os.path.join(sfile_folder,"select.out ")
    os.system(f'head -1 {out}')

def norhead(sfile_folder):
    """
    Parameters:
    -----------
    sfile_folder: str
        It is the folder where is located all sfiles.
            
    Returns:
        It select the principal header information of the event in
        only one file
    """
    norhead_path = os.path.join(PRO_path,"norhead")
    select_out = os.path.join(sfile_folder,"select.out")
    hypo_out = os.path.join(sfile_folder,"hypo.out")
    msg = f"{norhead_path} {select_out} {hypo_out}"
    # print(msg)
    os.system(msg)

def cp_station0(station0,sfile_folder):
    """
    Parameters:
    -----------
    station0: str
        Path of the station0 file
    sfile_folder: str
        It is the folder where is located all sfiles.
            
    Returns:
        Copy the station0 file in the sfile_folder
    """
    station0_path = os.path.join(DAT_path,"STATION0.HYP")
    if os.path.isfile(station0) != True:

        msg = f"cp {station0} {sfile_folder}"
        # print(msg)
        os.system(msg)

def resp2df(resp):
    """
    Parameters:
    -----------
    resp: str
        RESP filepath

    Returns: DataFrame
        Dataframe with the next columns
        network,station,latitude,longitude,elevation
    """
    networks = []
    stations = []
    longitudes = []
    latitudes = []
    elevations = []
    inv = read_inventory(resp)
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
            elevations.append(sta.elevation)
            stations.append(sta.code)
            networks.append(net.code)

    df = {"network":networks,"station":stations,
        "latitude":latitudes,"longitude":longitudes,
        "elevation":elevations}
    df = pd.DataFrame(df)
    # print(df)
    return df

def sta2station0(df):
    """
    Parameters:
    ----------
    df: DataFrame
        Dataframe with the next columns
        network,station,latitude,longitude,elevation.
        Review resp2df function.

    Returns: 
    --------
    msgs: list
        List of each message by station.
        The station message contains the next information:
        
        station lat lon elevation
        ---example---
        CA01A 353.16N 7341.07W 436

    """
    df = df.sort_values(by="station")

    df["lon_dec"] = df["longitude"].apply(lambda x: round((abs(x) % 1)*60,2))
    df["lat_dec"] = df["latitude"].apply(lambda x: round((abs(x) % 1)*60,2))
    df["latitude"] = df["latitude"].apply(lambda x: round(int(x),2))
    df["longitude"] = df["longitude"].apply(lambda x: round(int(x),2))

    msgs = []
    for i,row in df.iterrows():
        sta = row.station
        lat = str(abs(row.latitude)).zfill(2)
        lon = str(abs(row.longitude)).zfill(2)

        if row.latitude < 0:
            latcord = "S"
        else:
            latcord = "N"

        if row.longitude < 0:
            loncord = "W"
        else:
            loncord = "E"

        fmt = "%6s%2i%5.2f%1s%3i%5.2f%1s%4i\n"
        msg = fmt % (sta,int(lat),float(row.lat_dec),
                    latcord,int(lon),float(row.lon_dec),
                    loncord,int(row.elevation))
        msgs.append(msg)
        # print(msg)
    return msgs

def vel2station0(df,only_vp=True):
    """
    Parameters:
    -----------
    df: DataFrame
        Dataframe with the next columns
        dep,vp,vs,disc

    Returns:
    --------
    msgs: list
        List of each message by discontinuity in the velocity model.
        The discontinuity message contains the next information:
        
        vp dep vs disc
        ---example---
        7.000 25.000  3.940B

    Warnings:
    ---------
        Discontinuities
        B: Conrad
        N: Moho
    """
    df = df.sort_values(by="vp")
    df = df.where(pd.notnull(df), None)
    msgs = []
    for i,row in df.iterrows():
        if row.disc != None:
            if only_vp:
                fmt = "%7.3f%7.3f%1s\n"
                msg = fmt % (row.vp,row.depth,row.disc)
            else:
                fmt = "%7.3f%7.3f%7.3f%1s\n"
                msg = fmt % (row.vp,row.depth,row.vs,row.disc)
        else: 
            if only_vp:
                fmt = "%7.3f%7.3f\n"
                msg = fmt % (row.vp,row.depth)
            else:
                fmt = "%7.3f%7.3f%7.3f\n"
                msg = fmt % (row.vp,row.depth,row.vs)

        # print(msg)
        msgs.append(msg)
    return msgs

def test2station0(test):
    """
    Parameters:
    -----------
    test: dict
        key: str
            Number of the SEISAN HYPOCENTER test line
        value: int
            Value of the respective key

    Returns:
    --------
    msgs: list
        List of each message by test line.

        The test message contains the next information:
        RESET TEST(key)=value

        ------example-------
        RESET TEST(02)=500.0
        
    """
    msgs = []
    for key,value in test.items():
        key = str(key).zfill(2)
        msg = f"RESET TEST({key})={value}\n"
        # print(msg)
        msgs.append(msg)

    return msgs

def vsp_ratio2station0(vsp_ratio):
    """
    Parameters:
    -----------
    vsp_ratio: ordered dict
        The key order of the dictionary must be the next:
        "starting_depth","xnear","xfar","vps"

        ---------example--------
        {"starting_depth":2,
        "xnear":100,
        "xfar": 800,
        "vps": 1.84},

    Returns:
    --------
    msg: str
        Represents the last line of the STATION0 file

    """
    fmts = []
    values = []
    for key,value in vsp_ratio.items():
        if key == "vps":
            fmt = "5.2f"   
        else:
            fmt = "5.0f"
            
        fmts.append(fmt)
        values.append(float(value))

    fmt = "%"+"%".join(fmts)
    # print(fmt)
    # print(tuple(values))
    msg = fmt % tuple(values) + "\n"
    # print(msg)
    return msg

def check_sfile_integrity(sfile_folder,rm_not_locatable=True):
    """
    It checks sfile integrity.
    We aim that obspy could be read it as catalog.

    Parameters:
    -----------
    sfile_folder: str
        It is the folderpath where is located all sfiles.
    rm_not_locatable: bool
        Removes not locatable events
    """
    for dp, dn, filenames in os.walk(sfile_folder):
        for f in filenames:
            # print(f)
            #check if is an sfile to collect
            if f.split(".")[-1][0] == "S":

                sfilepath = os.path.join(dp, f)

                with open(sfilepath,"r") as sf:
                    line0 = sf.readline()
                    # print(line0[0:10])
                    
                    latlondep_col = line0[23:45]

                    if latlondep_col.isspace():
                        print(sfilepath,"not_locatable")
                    
                        if rm_not_locatable:
                            os.remove(sfilepath)
                            print("Removed:",sfilepath,"not_locatable")

class STATION0():
    def __init__(self,xml_path, vel_path,
            test={"02":500.0,"11":99.0,
                "13":1.0,"41":20000.0,
                "43":5.0,"56":1.0,
                "85":0.1},
            vsp_ratio = {"starting_depth":2,
                        "xnear":100,
                        "xfar": 800,
                        "vps": 1.78},
            agency = "TES",
            only_vp = True
                ):
        """
        Parameters:
        -----------
        sta_df: DataFrame
            Dataframe with the next columns
            network,station,latitude,longitude,elevation.
            Review resp2df function.
        vel_df: DataFrame
            Dataframe with the next columns
            dep,vp,vs,disc
        test: dict
            key: str
                Number of the SEISAN HYPOCENTER test line
            value: int
                Value of the respective key
        vsp_ratio: ordered dict
            The key order of the dictionary must be the next:
            "starting_depth","xnear","xfar","vps"

            ---------example--------
            {"starting_depth":2,
            "xnear":100,
            "xfar": 800,
            "vps": 1.84},
        agency: str
            Agency 
        """

        self.test = test
        self.sta_df = resp2df(xml_path)
        self.vel_df = pd.read_csv(vel_path)
        self.vsp_ratio = vsp_ratio
        self.agency = agency
        self.only_vp = only_vp
    
    def _get_msgs(self):
        """
        Get the messages to write in the station0 file
        """
        test_msgs = test2station0(self.test) 
        sta_msgs = sta2station0(self.sta_df)
        vel_msgs = vel2station0(self.vel_df,only_vp=self.only_vp)
        vsp_msg = vsp_ratio2station0(self.vsp_ratio)


        msgs = test_msgs + ["\n"] +sta_msgs +\
                ["\n"] + vel_msgs + ["\n"] +\
                [vsp_msg] + [self.agency]
        # print(msgs)

        return msgs

    def write(self,out):
        """
        Parameters:
        -----------
        out: str
            Station0 path file
        """
        # ut.isfile(out)
        f = open(out, "w")
        for line in self._get_msgs():
            f.write(line)
        f.close()

class HypocenterTools():
    def __init__(self,sfile_folder,from_sfilename=None):
        """
        Parameters:
        -----------
        sfile_folder: str
            It is the folderpath where is located all sfiles.
        station0: str
            Station0 path
        from_sfilename: str
            Name of the only one sfile with specific filename
            that you want to locate.
        """
        self.sfile_folder = sfile_folder
        self.sfilename = from_sfilename

    def remodl_and_setbrn(self):
        """
        Copy station0 in sfiles folder
        """
        subprocess.Popen(["remodl","setbrn"],cwd=self.sfile_folder)

        # cp_station0(self.station0_base,self.sfile_folder)
    
    def update(self):
        """
        Locate the events with the update seisan command
        """
        update(self.sfile_folder)
    
    def collect(self):
        """
        Collect all sfiles in only one collect sfile and write it
        in {out} path with the {format} specified.

        Parameters:
        -----------
        out: str
            Output Path of the collected catalog
        format: str
            Format supported by OBspy catalogs
        """
        collect(self.sfile_folder)
        nordic = os.path.join(self.sfile_folder,"collect.out")
        ## check if nordic event is ok
        catalog = read_events(nordic,format="NORDIC")
        ## put network and channel in catalog

        # # ut.isfile(out)
        # catalog.write(out,format=format)
        return catalog

    def select(self):
        """
        Runs select seisan command
        """
        select(self.sfile_folder)

    def norhead(self):
        """
        Runs norhead seisan command
        """
        norhead(self.sfile_folder)

    def split(self):
        """
        Runs split seisan command
        """
        split(self.sfile_folder,self.sfilename)



if __name__ == "__main__":
    # download_seisan()

    # xml_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
    # vel_path = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"

    # sta0 = STATION0(xml_path,vel_path)
    # sta0.write("./STATION0.HYP")

    # seisan_folder= "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/hypocenter/test/"
    # sfilename = "12-2223-00L.S202106"
    # split(seisan_folder,sfilename )
    # update(seisan_folder )


    ####### cse
    xml_path = "/home/emmanuel/G-Ecopetrol/ecastillo/Avances/2022/cano_sur_datos_iniciales/CSE.xml"
    vel_path = "/home/emmanuel/G-Ecopetrol/ecastillo/Avances/2022/cano_sur_datos_iniciales/vel_model.csv"

    sta0 = STATION0(xml_path,vel_path)
    sta0.write("./STATION0.HYP")