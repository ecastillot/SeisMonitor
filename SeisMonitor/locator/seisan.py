# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:52:18
#  * @modify date 2021-12-22 09:52:18
#  * @desc [description]
#  */


import os
import pexpect
import subprocess
import pandas as pd
from obspy import read_inventory
from obspy.core.event.catalog import Catalog, read_events
import SeisMonitor.utils as ut

# https://seis.geus.net/software/seisan/node243.html : UNTERSTANDING NORDIC FILE


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
    Split = pexpect.spawn('split', cwd=sfile_folder)
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
    Collect = pexpect.spawn('collect', cwd=sfile_folder)
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
    Update = pexpect.spawn('update', cwd=sfile_folder)
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
    Select = pexpect.spawn('select', cwd=sfile_folder)
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
    select_out = os.path.join(sfile_folder,"select.out")
    hypo_out = os.path.join(sfile_folder,"hypo.out")
    msg = f"norhead {select_out} {hypo_out}"
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

def vel2station0(df):
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
            fmt = "%7.3f%7.3f%7.3f%1s\n"
            msg = fmt % (row.vp,row.dep,row.vs,row.disc)
        else: 
            fmt = "%7.3f%7.3f%7.3f\n"
            msg = fmt % (row.vp,row.dep,row.vs)

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
    def __init__(self,sta_df, vel_df,
            test={"02":500.0,"11":99.0,
                "13":1.0,"41":20000.0,
                "43":5.0,"56":1.0,
                "85":0.1},
            vsp_ratio = {"starting_depth":2,
                        "xnear":100,
                        "xfar": 800,
                        "vps": 1.84},
            agency = "TES"
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
        self.sta_df =sta_df
        self.vel_df = vel_df
        self.vsp_ratio = vsp_ratio
        self.agency = agency
    
    def _get_msgs(self):
        """
        Get the messages to write in the station0 file
        """
        test_msgs = test2station0(self.test) 
        sta_msgs = sta2station0(self.sta_df)
        vel_msgs = vel2station0(self.vel_df)
        vsp_msg = vsp_ratio2station0(self.vsp_ratio)


        msgs = test_msgs + ["\n"] +sta_msgs +\
                ["\n"] + vel_msgs + ["\n"] +\
                [vsp_msg] + [self.agency]
        # print(msgs)

        return msgs

    def write(self,out="./STATION0.HYP"):
        """
        Parameters:
        -----------
        out: str
            Station0 path file
        """
        ut.isfile(out)
        f = open(out, "w")
        for line in self._get_msgs():
            f.write(line)
        f.close()

class Hypocenter():
    def __init__(self,sfile_folder,station0_base,from_sfilename=None):
        """
        Parameters:
        -----------
        sfile_folder: str
            It is the folderpath where is located all sfiles.
        station0_base: str
            Station0 path
        from_sfilename: str
            Name of the only one sfile with specific filename
            that you want to locate.
        """
        self.sfile_folder = sfile_folder
        self.station0_base = station0_base
        self.sfilename = from_sfilename

    def station0(self):
        """
        Copy station0 in sfiles folder
        """
        subprocess.Popen(["remodl","setbrn"],cwd=self.sfile_folder)
        cp_station0(self.station0_base,self.sfile_folder)
    
    def update(self):
        """
        Locate the events with the update seisan command
        """
        update(self.sfile_folder)
    
    def collect(self,out,format="NORDIC"):
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

        ut.isfile(out)
        catalog.write(out,format=format)

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

