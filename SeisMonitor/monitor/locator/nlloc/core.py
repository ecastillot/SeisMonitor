import os
import numpy as np
from typing import Union
import concurrent.futures as cf
from obspy.core.event.catalog import Catalog, read_events
from obspy.geodetics.base import gps2dist_azimuth
from obspy.io.nlloc.core import read_nlloc_hyp

from . import utils as ut
from SeisMonitor import utils as sut
from SeisMonitor.monitor.locator import utils as slut

CORE_NLLOC = os.path.join(os.path.dirname(__file__),"core")
NLLOC_path = os.path.join(CORE_NLLOC,"NonLinLoc-main")
src_path = os.path.join(NLLOC_path,"src")
bin_path = os.path.join(src_path,"bin")
vel2grid_exe_path = os.path.join(bin_path,"Vel2Grid")
grid2time_exe_path = os.path.join(bin_path,"Grid2Time")
nll_exe_path = os.path.join(bin_path,"NLLoc")

class NLLoc():
    def __init__(self,
                region:list,
                basic_inputs:slut.LocatorBasicInputs,
                delta_in_km: float=2,
                kwargs_for_trans:dict={},
                kwargs_for_vel2grid:dict={},
                kwargs_for_grid2time:dict={},
                kwargs_for_time2loc:dict={},
                tmp_folder:str = os.getcwd(),
                rm_tmp_folder:bool = False):

        self.region = region
        self.basic_inputs = basic_inputs
        self.delta_in_km = delta_in_km
        self.kwargs_for_trans = kwargs_for_trans
        self.kwargs_for_vel2grid = kwargs_for_vel2grid
        self.kwargs_for_grid2time = kwargs_for_grid2time
        self.kwargs_for_time2loc = kwargs_for_time2loc
        self.tmp_folder = tmp_folder
        self.tm_tmp_folder = rm_tmp_folder

    def __initialize(self):
        ### inputs
        self.vel_model_path = os.path.join(self.tmp_folder,"time_grid","vel_model.dat")
        self.station_path = os.path.join(self.tmp_folder,"time_grid","station.dat")
        sut.isfile(self.vel_model_path)
        self.basic_inputs.vel_model.to_nlloc(self.vel_model_path)
        sut.isfile(self.station_path)
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
        sut.isfile(self.catalog_path)
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

        ## The parallelization doesn't work, therefore I only use max_workers=1. I am too lazy to change the code. 
        with cf.ThreadPoolExecutor(max_workers=1) as executor:
            executor.map(write_each_nlloc_file,nlloc_control_files)

        self.nll_control_file = p_nlloc

    def _prepare_grid_args(self):
        lonw,lone,lats,latn,zmin,zmax = self.region

        x,_,_ = gps2dist_azimuth(lats,lonw,lats,lone)
        x = int(x/1e3/self.delta_in_km)
        y,_,_ = gps2dist_azimuth(lats,lonw,latn,lonw)
        y = int(y/1e3/self.delta_in_km)

        z = int((zmax-zmin)/self.delta_in_km)

        args = {"trans":["SIMPLE",lats,lonw,0],
                "velgrid":[2,y,z,
                            0,0,-5.0,
                            self.delta_in_km,self.delta_in_km,
                            self.delta_in_km,"SLOW_LEN"],
                "locgrid":[x,y,z,
                            0,0,-5.0,
                            self.delta_in_km,self.delta_in_km,
                            self.delta_in_km,
                            "PROB_DENSITY","SAVE"]
                            }
        return args

    def compute_travel_times(self):
        self.__initialize()

        sut.printlog("info","NLLoc:Vel2Grid", "Running")
        vel2grid = os.system(f"{vel2grid_exe_path} {self.p_control_file_out} > /dev/null")
        sut.printlog("info","NLLoc:Grid2Time:P", "Running")
        grid2time = os.system(f"{grid2time_exe_path} {self.p_control_file_out} > /dev/null")
        sut.printlog("info","NLLoc:Grid2Time:S", "Running")
        grid2time = os.system(f"{grid2time_exe_path} {self.s_control_file_out} > /dev/null")

    def relocate(self,
                nlloc_out_folder:str,
                out_format:str = "NORDIC"):
        

        nlloc_inp = os.path.join(nlloc_out_folder,"catalog_input.out")
        nlloc_folder = os.path.join(nlloc_out_folder,"nlloc","SeisMonitor")
        nlloc_out = os.path.join(nlloc_out_folder,"catalog_output.out")
        nlloc_control = os.path.join(nlloc_out_folder,"loc.in")

        sut.isfile(nlloc_inp)
        self.basic_inputs.catalog.write(nlloc_inp,
                                        format="NORDIC")
        
        
        nlloc_control_file = self.nll_control_file
        nlloc_control_file.time2loc.catalog = " ".join((nlloc_inp,"SEISAN")) 
        nlloc_control_file.time2loc.loc_folder_out = nlloc_folder 
        nlloc_control_file.write(nlloc_control)

        sut.printlog("info","NLLoc:NLLoc", "Running")
        os.system(f"{nll_exe_path} {nlloc_control} > /dev/null")

        _nll_out = os.path.join(os.path.dirname(nlloc_folder),"last.hyp")
        catalog = read_nlloc_hyp(_nll_out,format="NORDIC")

        sut.isfile(nlloc_out)
        catalog.write(nlloc_out,
                    format=out_format)
        sut.printlog("info","NLLoc:NLLoc", f"Finished. See your results in {nlloc_out}")
        return catalog
        