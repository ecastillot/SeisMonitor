import os
import shutil
import sys
from time import time
import pandas as pd
from subprocess import run
from obspy import read_inventory
from tqdm import tqdm
import glob
import SeisMonitor.utils as sut

CORE_NLLOC = os.path.join(os.path.dirname(__file__),"core")
NLLOC_path = os.path.join(CORE_NLLOC,"NonLinLoc-main")
src_path = os.path.join(NLLOC_path,"src")
bin_path = os.path.join(src_path,"bin")
vel2grid_exe_path = os.path.join(bin_path,"Vel2Grid")
grid2time_exe_path = os.path.join(bin_path,"Grid2Time")
nll_exe_path = os.path.join(bin_path,"NLLoc")

def run_nlloc(p_control_file_path,
            s_control_file_path):

    sut.printlog("info","NLLoc:Vel2Grid", "Running")
    vel2grid = os.system(f"{vel2grid_exe_path} {p_control_file_path} > /dev/null")
    sut.printlog("info","NLLoc:Grid2Time:P", "Running")
    grid2time = os.system(f"{grid2time_exe_path} {p_control_file_path} > /dev/null")
    sut.printlog("info","NLLoc:Grid2Time:S", "Running")
    grid2time = os.system(f"{grid2time_exe_path} {s_control_file_path} > /dev/null")
    sut.printlog("info","NLLoc:NLLoc", "Running")
    grid2time = os.system(f"{nll_exe_path} {s_control_file_path} > /dev/null")

def apt_install(pkgs):
    cmd = ['pkexec', 'apt-get', 'install', '-y'] + pkgs
    print('Running command: {}'.format(' '.join(cmd)))
    result = run(
        cmd,
        stdout=sys.stdout,
        stderr=sys.stderr,
        encoding='utf8',
        env={**os.environ, 'DEBIAN_FRONTEND': 'noninteractive'}
    )
    result.check_returncode()

def write_pref_origin_removing_phaselocinfo(catalog):
    events = []
    for ev in catalog:
        pref_origin = ev.preferred_origin()

        # replace none by NN
        for pick in ev.picks:
            pick.waveform_id.network_code= "NN"
            pick.waveform_id.location_code= "NN"
            pick.waveform_id.channel_code= "NNN"
            pick.evaluation_mode = "automatic"


        new_arrivals = []
        for arrival in pref_origin.arrivals:
            arrival.azimuth = None
            arrival.distance = None
            arrival.takeoff_angle=None
            arrival.time_residual=None
            arrival.time_weight=None
            new_arrivals.append(arrival)

        for i,origin in enumerate(ev.origins):
            if origin.resource_id.id == ev.preferred_origin_id:
                ev.origins[i].arrivals = new_arrivals 
            # else:
            #     continueagency_id=self.agency

        del ev.origins
        ev.origins = [pref_origin]
        events.append(ev)
    catalog.events = events
    return catalog

def download_nlloc(forced=False):
    name = "nll.zip"
    zip_path = os.path.join(CORE_NLLOC,name)
    cache_path = os.path.join(src_path,"CMakeCache.txt")
    # "https://github.com/alomax/NonLinLoc/archive/refs/heads/main.zip"

    if not forced:
        if os.path.isfile(nll_exe_path):
            return True
        else:
            pass

    isfile = sut.isfile(zip_path)
    if not isfile:
        os.system(f"wget https://github.com/alomax/NonLinLoc/archive/refs/heads/main.zip -O {zip_path}")

    if not os.path.isdir(NLLOC_path):
        os.system(f"unzip {zip_path} -d {CORE_NLLOC}")

    if os.path.isdir(bin_path):
        os.rmdir(bin_path)
        os.makedirs(bin_path)

    if os.path.isdir(cache_path):
        os.rmdir(cache_path)

    try:
        apt_install(["cmake"])
    except:
        sut.printlog("warning","Install Cmake","Could not install Cmake")

    try:
        os.system(f"cd {src_path} && cmake . && make")
    except:
        raise Exception("Could not compile NLLoc")

def write_1d_vel_model(vel_path,out,
                        compute_vs=True,
                        vp_vs_ratio=1.78):
    df = pd.read_csv(vel_path)

    vm = open(out, 'w')
    msg ="# model layers (LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad)\n"
    vm.write(msg)
    for i,row in df.iterrows():
        if compute_vs:
            vs = row.vp/vp_vs_ratio
        else:
            vs = row.vs

        if i == len(df)-1:
            enter = ""
        else:
            enter = "\n"

        msg = f"LAYER    {row.depth:<6.2f}    {row.vp:<.2f}    0.00    {vs:<.2f}    0.00    {row.rho:<.2f}    0.00{enter}"
        vm.write(msg)
    vm.close()

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
    return df 

def write_station_file(sta_path,out):
    df = resp2df(sta_path)

    vs = open(out, 'w')
    msg = "# GTSRCE label LATLON latSrce longSrce zSrce elev\n"
    vs.write(msg)
    for i,row in df.iterrows():
        elv = row.elevation/1e3
        if i == len(df)-1:
            enter = ""
        else:
            enter = "\n"

        msg = f"GTSRCE  {row.station:<5}  LATLON  {row.latitude:<7.3f}  {row.longitude:<8.3f}  0.000  {elv:<7.3f}{enter}"
        vs.write(msg)
    vs.close()

def join_args(args:str):
    return " ".join(list(map(lambda x: str(x),args)))
    
class GenericControlStatement():
    def __init__(self,
                trans:list,
                control:list=[1,54321]):
        self.control = join_args(control)
        self.trans = join_args(trans)
        sut.validate(self.__init__, locals())

    def get_msg(self):
        msg = (
            "#__________________START GENERIC CONTROL STATEMENTS\n\n"
            "CONTROL {}\n",
            "TRANS {}\n",
            "#__________________END\n"

            )
        msg = "".join(msg)
        msg = msg.format(self.control,
                    self.trans)
        sut.printlog("debug","NLLoc:GenericControlStatement",
                    "Message received")
        return msg

class Vel2Grid():
    def __init__(self,
                vel_path:str,
                grid_folder_out:str,
                grid:list,
                p_phase:str="P",
                s_phase:str="S") -> None:
        self.vel_path = vel_path
        self.grid = join_args(grid)
        self.grid_folder_out = grid_folder_out
        self.p_phase = p_phase
        self.s_phase = s_phase
        sut.validate(self.__init__, locals())

    def get_msg(self):
        msg = (
            "#__________________START VEL2GRID STATEMENTS\n\n"
            "VGOUT {}\n",
            "VGTYPE {}\n",
            "VGTYPE {}\n",
            "VGGRID {}\n",
            "INCLUDE {}\n",
            "#__________________END\n"

            )
        msg = "".join(msg)
        msg = msg.format(self.grid_folder_out,
                        self.p_phase,
                        self.s_phase,
                        self.grid,
                        self.vel_path)

        path = os.path.dirname(self.grid_folder_out)
        sut.isfile(path)

        sut.printlog("debug","NLLoc:Vel2grid",
                    "Message received")
        return msg

class Grid2Time():
    def __init__(self,
        station_path:str,
        grid_folder_out:str,
        time_folder_out:str,
        phase:str="P",
        mode:list=["GRID3D","ANGLES_YES"],
        plfd:list=[1e-3,0]) -> None:

        self.station_path = station_path 
        self.grid_folder_out = grid_folder_out 
        self.time_folder_out = time_folder_out 
        self.phase = phase
        self.mode = join_args(mode) 
        self.plfd = join_args(plfd) 
        sut.validate(self.__init__, locals())

    def get_msg(self):
        msg = (
            "#__________________START GRID2TIME STATEMENTS\n\n",
            "GTFILES {} {} {}\n",
            "GTMODE {}\n",
            "INCLUDE {}\n",
            "GT_PLFD {}\n"
            "#__________________END\n"
            )
        msg = "".join(msg)
        msg = msg.format(self.grid_folder_out,self.time_folder_out,self.phase,
                self.mode,
                self.station_path,
                self.plfd 
                )

        path = os.path.dirname(self.time_folder_out)
        sut.isfile(path)

        sut.printlog("debug","NLLoc:Grid2Time",
                    "Message received")
        return msg

class Time2Loc():
    def __init__(self,
        catalog:list,
        grid:list,
        time_folder_out:str,
        loc_folder_out:str,
        meth:list = ["GAU_ANALYTIC",9999,4,-1,-1,1.78,6],
        # meth:list = ["EDT_OT_WT",9999,4,-1,-1,1.78,6],
        search:list = ["OCT",37,58,7,1e-2,int(1e5),int(1e4)],
        sig:str = "SeisMonitor",
        com:str = "Comment",
        gau:list = [0.2,0.0],
        gau2:list = [0.05,0.05,2.0],
        p_phaseid:list = ["P","P","p","PN","PG","Pn","Pg"],
        s_phaseid:list = ["S","S","s","SN","SG","Sn","Sg"],
        qual2err:list = [0.1, 0.5, 1.0, 2.0, 99999.9],
        phstat:list = [9999.0,-1,9999.0,1.0,1.0,9999.0,-9999.0,9999.0],
        angles:list = ["ANGLES_YES",5],
        hypout:list = ["SAVE_NLLOC_ALL","SAVE_NLLOC_SUM",
                    "SAVE_HYPO71_SUM"],
        mag:list = ["ML_HB",1.0,1.110,0.00189]
        ):
        self.catalog = join_args(catalog)
        self.time_folder_out = time_folder_out
        self.loc_folder_out = loc_folder_out
        self.grid = join_args(grid)
        self.meth = join_args(meth)
        self.search = join_args(search)
        self.sig = sig
        self.com = com
        self.gau = join_args(gau)
        self.gau2 = join_args(gau2)
        self.p_phaseid = join_args(p_phaseid)
        self.s_phaseid = join_args(s_phaseid)
        self.qual2err = join_args(qual2err)
        self.phstat = join_args(phstat)
        self.angles = join_args(angles)
        self.hypout = join_args(hypout)
        self.mag = join_args(mag)
        sut.validate(self.__init__, locals())

    def get_msg(self):
        msg = (
            "#__________________START NLDIFFLOC STATEMENTS\n\n"
            "LOCSIG {}\n",
            "LOCCOM {}\n",
            "LOCFILES {} {} {}\n",
            "LOCHYPOUT {}\n",
            "LOCSEARCH {}\n",
            "LOCGRID {}\n",
            "LOCMETH {}\n",
            "LOCGAU {}\n",
            "LOCGAU2 {}\n",
            "LOCPHASEID {}\n",
            "LOCPHASEID {}\n",
            "LOCQUAL2ERR {}\n",
            "LOCPHSTAT {}\n",
            "LOCANGLES {}\n",
            "LOCMAG {}\n",
            "#__________________END\n")
            
        msg = "".join(msg)
        msg = msg.format(
                        self.sig,
                        self.com,
                        self.catalog,self.time_folder_out,self.loc_folder_out,
                        self.hypout,
                        self.search,
                        self.grid,
                        self.meth,
                        self.gau,
                        self.gau2,
                        self.p_phaseid,
                        self.s_phaseid,
                        self.qual2err,
                        self.phstat,
                        self.angles,
                        self.mag)

        path = os.path.dirname(self.loc_folder_out)
        sut.isfile(path)

        sut.printlog("debug","NLLoc:Time2Loc",
                    "Message received")
        return msg

class NLLocControlFile():
    def __init__(self,
                generic_control:GenericControlStatement,
                vel2grid:Vel2Grid,
                grid2time:Grid2Time,
                time2loc:Time2Loc):
        self.generic_control = generic_control
        self.vel2grid = vel2grid
        self.grid2time = grid2time
        self.time2loc = time2loc
        sut.validate(self.__init__, locals())

    def _validate_args(self,key:str,input:str,output:str):
        sut.printlog("debug","NLLoc:validate_control_file_args",
                    "validating control file args")

        if os.path.isfile(input):
            sut.printlog("debug","NLLoc:validate_control_file_args",
                f"{key} is ok.")
        else:
            raise Exception(f"NLLoc: {key} path doesn't exist. Check:{input}")

        sut.isfile(output)

        sut.printlog("debug","NLLoc:validate_control_file_args",
                    "control file args are ok")

    def get_msg(self,
            vel2grid:bool=True,
            grid2time:bool=True,
            time2loc:bool=True):

        msg = self.generic_control.get_msg()

        if vel2grid:
            msg += self.vel2grid.get_msg()
            input = self.vel2grid.vel_path
            output = self.vel2grid.grid_folder_out
            self._validate_args("vel_path",input,output)
        if grid2time:
            input = self.grid2time.station_path
            output = self.grid2time.time_folder_out
            self._validate_args("station_path",input,output)
            msg += self.grid2time.get_msg()
        if time2loc:
            input = self.time2loc.catalog.split(" ")[0]
            output = self.time2loc.loc_folder_out
            self._validate_args("catalog_path",input,output)
            msg += self.time2loc.get_msg()

        return msg

    def write(self,
            out:str,
            vel2grid:bool=True,
            grid2time:bool=True,
            time2loc:bool=True):

        sut.printlog("info","NLLoc:write_control_file",
                    "Running")
        sut.isfile(out)

        msg = self.get_msg(vel2grid,grid2time,time2loc)

        control_file_msg = open(out,"w")
        control_file_msg.write(msg)
        control_file_msg.close()
        sut.printlog("info","NLLoc:write_control_file",
                    f"Finished. Control file: {out} ")

if __name__=="__main__":
    catalog = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/select.out"
    station_path = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/station.dat"
    vel_path = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/model.dat"
    grid_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/model/layer"
    time_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/time/layer"
    loc_folder_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/loc/SeisMonitor"
    control_file_out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/test.in"

    gen_control = GenericControlStatement(trans=["SIMPLE",
                                             5.0,-73.0,0.0])
    vel2grid = Vel2Grid(vel_path,
                        grid_folder_out,
                        grid = [2,583,70,
                                -891.0,-891.0,-5.0,
                                3.0, 3.0, 3.0,
                                "SLOW_LEN"],
                        phase = "P")
    grid2time = Grid2Time(station_path,
                            grid_folder_out,
                            time_folder_out)
    time2loc = Time2Loc(catalog=[catalog,"SEISAN"],
                        grid = [374,583,70,
                                -891.0,-891.0,-5.0,
                                3.0,3.0,3.0,
                                "PROB_DENSITY","SAVE"],
                        time_folder_out=time_folder_out,
                        loc_folder_out=loc_folder_out)
    nlloc = NLLocObj(gen_control,vel2grid,
                    grid2time,time2loc)
    nlloc.write_control_file(control_file_out)
