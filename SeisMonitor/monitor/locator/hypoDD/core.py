import os
import pandas as pd
from obspy.core.event.catalog import Catalog, read_events
from SeisMonitor.utils import printlog, isfile
from . import utils as ut

class HypoDD():
    def __init__(self,catalog,xml_path,vel_path,out_dir):

        if isinstance(catalog,Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)
        
        self.xml_path = xml_path
        self.vel_path = vel_path
        self.out_dir = out_dir
        self.out_file = os.path.join(self.out_dir,"hypodd","hypodd_events.xml")
        self.paths = os.path.join(self.out_dir,"hypodd","hypodd_paths")
        if not os.path.isdir(self.paths):
            os.makedirs(self.paths)
        self.pha = os.path.join(self.paths,"hypoDD.pha")
        self.vel_df = pd.read_csv(vel_path)


    def locate(self,
                vp_vs_ratio = 1.84,
                out_format="QUAKEML",
                rm_not_locatable=True):
        df = ut.resp2df(self.xml_path)
        ut.write_hypoDDstation(df,self.paths)
        ut.write_pha(self.catalog,self.pha)
        ut.write_ph2dt_inp_file(self.paths)

        vel_layers = ut.get_vel_layers(self.vel_df)
        vel_model = ut.setup_velocity_model("layered_p_velocity_with_constant_vp_vs_ratio",
                        vp_vs_ratio=1.82,
                        layer_tops=vel_layers
                        )
        ut.write_hypoDD_inp_file(vel_model,self.paths)

        # print(ut.CORE_HYPODD)

        pos = "/home/emmanuel/QuakeFlow/HypoDD"
        ph2dt_path = os.path.join(pos,"HYPODD","src","ph2dt","ph2dt")
        hypodd_path = os.path.join(pos,"HYPODD","src","hypoDD","hypoDD")

        ph2dt_inp_path = os.path.join(self.paths,"ph2dt.inp")
        hypodd_inp_path = os.path.join(self.paths,"hypoDD.inp")
        PH2DT_CMD = f"cd {self.paths} && {ph2dt_path} ph2dt.inp"
        HYPODD_CMD = f"cd {self.paths} && {hypodd_path} hypoDD.inp"
        # print(PH2DT_CMD)
        # print(HYPODD_CMD )

        if os.system(PH2DT_CMD) != 0:
            raise("{PH2DT_CMD}" + " failed!")
        if os.system(HYPODD_CMD) != 0:
            raise("{HYPODD_CMD}" + " failed!")
