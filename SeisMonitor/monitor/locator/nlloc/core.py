import os
import numpy as np
from typing import Union
from obspy.core.event.catalog import Catalog, read_events
from obspy.geodetics.base import gps2dist_azimuth
from obspy.io.nlloc.core import read_nlloc_hyp

from . import utils as ut
from SeisMonitor import utils as sut
from SeisMonitor.monitor.locator import utils as slut

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


        self.vel_model_out = os.path.join(tmp_folder,"time_grid","vel_model.dat")
        self.station_out = os.path.join(tmp_folder,"time_grid","station.dat")
        self.catalog_out = os.path.join(tmp_folder,"catalog.out")
        self.p_control_file_out = os.path.join(tmp_folder,"time_grid","p_nlloc.in")
        self.s_control_file_out = os.path.join(tmp_folder,"time_grid","s_nlloc.in")
        self.grid_folder_out = os.path.join(tmp_folder,"time_grid","model","layer")
        self.time_folder_out = os.path.join(tmp_folder,"time_grid","time","layer")
        self.loc_folder_out = os.path.join(tmp_folder,"time_grid","loc","SeisMonitor")

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

    def _prepare_basic_inputs(self):
        
        
        sut.isfile(self.vel_model_out)
        self.basic_inputs.vel_model.to_nlloc(self.vel_model_out)
        sut.isfile(self.station_out)
        self.basic_inputs.stations.to_nlloc(self.station_out)
        sut.isfile(self.catalog_out)
        self.basic_inputs.catalog.write(self.catalog_out,
                                        format="NORDIC")
        return self.vel_model_out,self.station_out,self.catalog_out


    def _write_control_files(self):

        
        vel_model_path,station_path,catalog_path = self._prepare_basic_inputs()
        
        grid_args = self._prepare_grid_args()

        if "trans" in list(self.kwargs_for_trans.keys()):
            grid_args["trans"] = self.kwargs_for_trans["trans"]
        if "grid" in list(self.kwargs_for_vel2grid.keys()):
            grid_args["velgrid"] = self.kwargs_for_vel2grid["grid"]
        if "grid" in list(self.kwargs_for_grid2time.keys()):
            grid_args["locgrid"] = self.kwargs_for_grid2time["grid"]


        gen_control = ut.GenericControlStatement(trans=grid_args["trans"])
        vel2grid = ut.Vel2Grid(vel_path=vel_model_path,
                    grid_folder_out=self.grid_folder_out,
                    grid=grid_args["velgrid"])
        p_grid2time = ut.Grid2Time(station_path=station_path,
                            grid_folder_out=self.grid_folder_out,
                            time_folder_out=self.time_folder_out,
                            phase="P")
        s_grid2time = ut.Grid2Time(station_path=station_path,
                            grid_folder_out=self.grid_folder_out,
                            time_folder_out=self.time_folder_out,
                            phase="S")

        time2loc = ut.Time2Loc(catalog=[catalog_path,"SEISAN"],
                    grid = grid_args["locgrid"],
                    time_folder_out=time_folder_out,
                    loc_folder_out=loc_folder_out)

        p_nlloc = ut.NLLocControlFile(gen_control,vel2grid,
                                    p_grid2time,time2loc)
        p_nlloc.write(p_control_file_out)
        s_nlloc = ut.NLLocControlFile(gen_control,vel2grid,
                        s_grid2time,time2loc)
        s_nlloc.write(s_control_file_out)

        return p_control_file_out,s_control_file_out

    def compute_time_grids(self):


    def relocate(self,
                out:str = None,
                out_format:str = "NORDIC"):
        
        # p_control_file_path,s_control_file_path = self._write_control_files(tmp_folder)
        # ut.run_nlloc(p_control_file_path,
        #                 s_control_file_path)

        loc_folder_out = "/home/emmanuel/EDCT/test_nlloc/loc/last.hyp"
        catalog = read_nlloc_hyp(loc_folder_out,format="NORDIC")
        print(catalog[0])
        print(catalog)
        