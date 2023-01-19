import os
import glob
from posixpath import basename
import shutil
import numpy as np
from typing import Union
import concurrent.futures as cf
from obspy.core.event.base import CreationInfo
from obspy.core.event.catalog import Catalog, read_events
from obspy.geodetics.base import gps2dist_azimuth
from obspy.io.nlloc.core import read_nlloc_hyp
from obspy.core.event.base import CreationInfo
from obspy.core.event.resourceid import ResourceIdentifier
from obspy import UTCDateTime
import subprocess
from tqdm import tqdm
from . import utils as ut
from SeisMonitor import utils as sut
from SeisMonitor.core import utils as scut
from SeisMonitor.monitor.locator import utils as slut

CORE_NLLOC = os.path.join(os.path.dirname(__file__),"core")
NLLOC_path = os.path.join(CORE_NLLOC,"NonLinLoc-main")
src_path = os.path.join(NLLOC_path,"src")
bin_path = os.path.join(src_path,"bin")
vel2grid_exe_path = os.path.join(bin_path,"Vel2Grid")
grid2time_exe_path = os.path.join(bin_path,"Grid2Time")
nll_exe_path = os.path.join(bin_path,"NLLoc")


#https://github.com/SEISAN-EARTHQUAKE-ANALYSIS-SOFTWARE/SEISAN-manual/blob/master/appendix/nordic-format.tex seisan format
class NLLoc():
    def __init__(self,
                agency:str,
                region:list,
                vel_model:slut.VelModel,
                stations:slut.Stations,
                delta_in_km: float=2,
                kwargs_for_trans:dict={},
                kwargs_for_vel2grid:dict={},
                kwargs_for_grid2time:dict={},
                kwargs_for_time2loc:dict={},
                tmp_folder:str = os.getcwd(),
                exhaustively:bool = False,
                search_in_degrees:list = [],
                rm_attempts:bool = False
                ):

        self.agency = agency
        self.region = region
        self.vel_model = vel_model
        self.stations = stations
        self.basic_inputs = slut.LocatorBasicInputs(vel_model=vel_model,
                                                    stations=stations)
        self.delta_in_km = delta_in_km
        self.kwargs_for_trans = kwargs_for_trans
        self.kwargs_for_vel2grid = kwargs_for_vel2grid
        self.kwargs_for_grid2time = kwargs_for_grid2time
        self.kwargs_for_time2loc = kwargs_for_time2loc
        self.tmp_folder = tmp_folder
        self.rm_attempts = rm_attempts

        self.exhaustively = exhaustively
        if exhaustively:
            #     lon_distance_in_deg = abs(region[0]-region[1])
            #     lat_distance_in_deg = abs(region[2]-region[3])
            #     r = np.sqrt(lon_distance_in_deg**2 + lat_distance_in_deg**2 )
            #     self.search_in_degrees = np.linspace(1/111,r/2,5)[::-1]
            if not search_in_degrees:
                self.search_in_degrees = list(reversed(np.arange(1,4,0.5)))
            else:
                self.search_in_degrees = search_in_degrees

    def __initialize(self,write_nlloc_files=True):
        ### inputs
        self.vel_model_path = os.path.join(self.tmp_folder,"time_grid","vel_model.dat")
        self.station_path = os.path.join(self.tmp_folder,"time_grid","station.dat")
        sut.isfile(self.vel_model_path,overwrite=True)
        self.basic_inputs.vel_model.to_nlloc(self.vel_model_path)
        sut.isfile(self.station_path,overwrite=True)
        self.basic_inputs.stations.to_nlloc(self.station_path)

        self.grid_folder_out = os.path.join(self.tmp_folder,"time_grid","model","layer")
        self.time_folder_out = os.path.join(self.tmp_folder,"time_grid","time","layer")
        self.loc_folder_out = os.path.join(self.tmp_folder,"time_grid","loc","SeisMonitor")
        self.p_control_file_out = os.path.join(self.tmp_folder,"time_grid","p_nlloc.in")
        self.s_control_file_out = os.path.join(self.tmp_folder,"time_grid","s_nlloc.in")

        grid_args = self._prepare_grid_args()

        if "trans" in list(self.kwargs_for_trans.keys()):
            grid_args["trans"] = self.kwargs_for_trans["trans"]
        if "grid" in list(self.kwargs_for_vel2grid.keys()):
            grid_args["velgrid"] = self.kwargs_for_vel2grid["grid"]
        if "grid" in list(self.kwargs_for_grid2time.keys()):
            grid_args["locgrid"] = self.kwargs_for_grid2time["grid"]

        gen_control = ut.GenericControlStatement(trans=grid_args["trans"])
        gen_vel2grid = ut.Vel2Grid(vel_path=self.vel_model_path,
                    grid_folder_out=self.grid_folder_out,
                    grid=grid_args["velgrid"])
        p_grid2time = ut.Grid2Time(station_path=self.station_path,
                            grid_folder_out=self.grid_folder_out,
                            time_folder_out=self.time_folder_out,
                            phase="P")
        s_grid2time = ut.Grid2Time(station_path=self.station_path,
                            grid_folder_out=self.grid_folder_out,
                            time_folder_out=self.time_folder_out,
                            phase="S")

        self.catalog_path = os.path.join(self.tmp_folder,"time_grid","catalog.out")
        sut.isfile(self.catalog_path,overwrite=True)
        open(self.catalog_path, "w")
        gen_time2loc = ut.Time2Loc(catalog=[self.catalog_path,"SEISAN"],
                    grid = grid_args["locgrid"],
                    time_folder_out=self.time_folder_out,
                    loc_folder_out=self.loc_folder_out)

        def write_each_nlloc_file(nlloc_and_out):
            nlloc,out = nlloc_and_out
            nlloc.write(out,
                    vel2grid=True,
                    grid2time=True,
                    time2loc=True)
            return True


        p_nlloc = ut.NLLocControlFile(gen_control,gen_vel2grid,
                                    p_grid2time,gen_time2loc)
        s_nlloc = ut.NLLocControlFile(gen_control,gen_vel2grid,
                        s_grid2time,gen_time2loc)
        nlloc_control_files = [(p_nlloc,self.p_control_file_out),
                                (s_nlloc,self.s_control_file_out)]

        if write_nlloc_files:
            ## The parallelization doesn't work properly, therefore I only use max_workers=1. I am too lazy to change the code. 
            with cf.ThreadPoolExecutor(max_workers=1) as executor:
                executor.map(write_each_nlloc_file,nlloc_control_files)

        self.nll_control_file = p_nlloc

    def _prepare_grid_args(self):
        lonw,lone,lats,latn,zmin,zmax = self.region

        c_lat = (lats+latn)/2
        c_lon = (lonw+lone)/2

        x,_,_ = gps2dist_azimuth(lats,lonw,lats,c_lon)
        x = x/1e3
        x_num = int(x*2/self.delta_in_km)
        y,_,_ = gps2dist_azimuth(lats,lonw,c_lat,lonw)
        y = y/1e3
        y_num = int(y*2/self.delta_in_km)

        z_num = int((zmax-zmin)/self.delta_in_km)

        args = {
                "trans":["SIMPLE",c_lat,c_lon,0],
                # "trans":["LAMBERT","WGS-84",lats,lonw,0,15,0],
                "velgrid":[x_num,y_num,z_num,
                            -round(x,2),-round(y,2),round(zmin,2),
                            self.delta_in_km,self.delta_in_km,
                            self.delta_in_km,"SLOW_LEN"],
                "locgrid":[x_num,y_num,z_num,
                            -round(x,2),-round(y,2),round(zmin,2),
                            self.delta_in_km,self.delta_in_km,
                            self.delta_in_km,
                            "PROB_DENSITY","SAVE"]
                            }
        return args

    def download(self):
        if not os.path.isdir(ut.NLLOC_path):
            ut.download_nlloc()
        else:
            print(f"NonLinLoc is located in {ut.NLLOC_path}")

    def compute_travel_times(self):
        self.__initialize()

        sut.printlog("info","NLLoc:Vel2Grid", "Running")
        # vel2grid = subprocess.call(f"{vel2grid_exe_path} {self.p_control_file_out} > /dev/null")
        vel2grid = subprocess.call(f"{vel2grid_exe_path} {self.p_control_file_out}",shell=True)
        sut.printlog("info","NLLoc:Grid2Time:P", "Running")
        grid2time = subprocess.call(f"{grid2time_exe_path} {self.p_control_file_out}",shell=True)
        sut.printlog("info","NLLoc:Grid2Time:S", "Running")
        grid2time = subprocess.call(f"{grid2time_exe_path} {self.s_control_file_out}",shell=True)

    def _locate(self,
                catalog:Union[Catalog,str],
                nlloc_out_folder:str,
                out_filename:str = "locations.xml",
                out_format:str = "SC3ML"):

        if isinstance(catalog,Catalog):
            pass
        else:
            catalog = read_events(catalog)

        nlloc_inp = os.path.join(nlloc_out_folder,"catalog_input.inp")
        nlloc_folder = os.path.join(nlloc_out_folder,"nlloc","SeisMonitor")
        nlloc_out = os.path.join(nlloc_out_folder,out_filename)
        nlloc_control = os.path.join(nlloc_out_folder,"loc.in")

        sut.isfile(nlloc_inp,overwrite=True)
        catalog.write(nlloc_inp,format="NORDIC")

        picks_from_unproc_catalog = slut.get_picks(catalog)
        # arrivals_from_unproc_catalog = slut.get_picks(catalog)

        try:
            nlloc_control_file = self.nll_control_file
        except:
            self.__initialize(False)
            nlloc_control_file = self.nll_control_file

        nlloc_control_file.time2loc.catalog = " ".join((nlloc_inp,"SEISAN")) 
        nlloc_control_file.time2loc.loc_folder_out = nlloc_folder 
        nlloc_control_file.write(nlloc_control)

        sut.printlog("info","NLLoc:NLLoc", "Running")
        subprocess.call(f"{nll_exe_path} {nlloc_control}",shell=True)
        
        _nll_out = nlloc_folder+"*.hyp"
        all_events = []
        for path in tqdm(glob.glob(_nll_out)):
            basename = os.path.basename(path)
            date = basename.split(".")[1]
            if (date == "sum") or (date=="last.hyp"):
                continue
            else:
                try:
                    catalog = read_nlloc_hyp(path,format="NORDIC")
                except:
                    print(f"Unread: {path}")
                    continue

                for ev in catalog.events:
                    ori_pref  = ev.preferred_origin()
                    ori_pref = scut.add_aditional_origin_info(ori_pref,
                                    agency=self.agency,
                                    method_id="NLLOC",
                                    earth_model_id=self.vel_model.model_name)
                    ev.preferred_origin_id = ori_pref.resource_id.id

                    ev = slut.changing_picks_info(ev,picks_from_unproc_catalog)
                    ev = scut.add_aditional_event_info(ev,agency=self.agency)
                    all_events.append(ev)

        catalog = Catalog(events = all_events)
        catalog = scut.add_aditional_catalog_info(catalog,agency=self.agency)

        sut.isfile(nlloc_out)
        catalog.write(nlloc_out,format=out_format)
        sut.printlog("info","NLLoc:NLLoc", f"Finished. See your results in {nlloc_out}")
        return catalog

    def _iterlocate(self,
                catalog:Union[Catalog,str],
                nlloc_out_folder:str,
                out_filename:str = "locations.xml",
                out_format:str = "SC3ML",
                degrees=[4,3,2.5,2,1.5,1],
                rm_attempts=False
                ):     
        nlloc_out = os.path.join(nlloc_out_folder,out_filename)
        reloc_catalog = self._locate(catalog,nlloc_out_folder,
                    out_filename="base.xml",out_format="SC3ML")
        
        good_evs = []
        bad_catalog = reloc_catalog
        iter = 0
        while True:
            good_events,bad_events = slut.get_bad_and_good_events(bad_catalog)
            good_evs.append(good_events)
            degree = degrees[iter]
            if not bad_events:
                if iter == 0:
                    print(f"convergence without iterate")
                else:
                    print(f"convergence in {iter} iterations, only stations with distance < {degrees[iter-1]} degrees")
                break
            else:
                bad_events = slut.filter_arrivals_by_distance(bad_events,degree)
                
                if not bad_events:
                    break
                
                bad_catalog = Catalog(bad_events)
                tmp_path = os.path.join(nlloc_out_folder,"tmp",f"deg_{degree}")
                bad_catalog = self._locate(bad_catalog,tmp_path,
                    out_filename=f"{iter}.xml",out_format="SC3ML")
                iter +=1
        good_evs = [ y for x in good_evs for y in x]

        if not good_evs:
            pass
        else:
            reloc_catalog.events = good_evs

        sut.isfile(nlloc_out)
        reloc_catalog.write(nlloc_out,
                    format=out_format)
        sut.printlog("info","NLLoc:NLLoc", f"Finished. See your results in {nlloc_out}")
        
        if rm_attempts:
            shutil.rmtree(os.path.join(nlloc_out_folder,"tmp"))

        return catalog

    def locate(self,
                catalog:Union[Catalog,str],
                nlloc_out_folder:str,
                out_filename:str = "locations.xml",
                out_format:str = "SC3ML"
                ):
        if self.exhaustively:
            cat = self._iterlocate(catalog=catalog,
                            nlloc_out_folder=nlloc_out_folder,
                            out_filename=out_filename,
                            out_format=out_format,
                            degrees=self.search_in_degrees,
                            rm_attempts=self.rm_attempts)
        else:
            cat = self._locate(catalog=catalog,
                            nlloc_out_folder=nlloc_out_folder,
                            out_filename=out_filename,
                            out_format=out_format)
        return cat