"""
NLL runner was mainly based on based on

https://github.com/saeedsltm/PyNLLRunner
"""

import sys
import os
import shutil
import glob
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)
from obspy import read_inventory
from subprocess import STDOUT, check_call
import SeisMonitor.utils as ut
from subprocess import run
import pandas as pd
CORE_NLLOC = os.path.join(os.path.dirname(__file__),"cor")

class LocObj():
    def __init__(self,
        nstd              = 1    ,    # number of standard deviation for plotting ellipses 1: %68 , 2: %95.
        lat_min           = 35.00,    # minimum latitude
        lat_max           = 37.00,    # maximum latitude
        lon_min           = 51.00,    # minimum longitude    
        lon_max           = 53.00,    # maximum longitude
        dep_max           = 25.00,    # maximum depth
        p_res_max         = 0.5  ,    # maximum p-residual 
        s_res_max         = 0.5  ,    # maximum s-residual
        max_dist          = 50   ,    # maximum distance for residual plot
        herr_max          = 10   ,    # maximum horizontal error
        zerr_max          = 10   ,    # maximum depth error
        rms_max           = 1.0  ,    # maximum rms 
        minds_max         = 20   ,    # maximum mindistance
        dep_hist_bw       = 2    ,    # depth histogram bin width
        her_hist_bw       = 0.25 ,    # herr histogram bin width
        zer_hist_bw       = 0.25 ,    # zerr histogram bin width
        rms_hist_bw       = 0.05 ,    # rms histogram bin width
    ):

        self.nstd              = nstd    # number of standard deviation for plotting ellipses 1: %68 , 2: %95.
        self.lat_min           = lat_min    # minimum latitude
        self.lat_max           = lat_max    # maximum latitude
        self.lon_min           = lon_min    # minimum longitude    
        self.lon_max           = lon_max    # maximum longitude
        self.dep_max           = dep_max    # maximum depth
        self.p_res_max         = p_res_max    # maximum p-residual 
        self.s_res_max         = s_res_max    # maximum s-residual
        self.max_dist          = max_dist    # maximum distance for residual plot
        self.herr_max          = herr_max    # maximum horizontal error
        self.zerr_max          = zerr_max    # maximum depth error
        self.rms_max           = rms_max    # maximum rms 
        self.minds_max         = minds_max    # maximum mindistance
        self.dep_hist_bw       = dep_hist_bw    # depth histogram bin width
        self.her_hist_bw       = her_hist_bw    # herr histogram bin width
        self.zer_hist_bw       = zer_hist_bw    # zerr histogram bin width
        self.rms_hist_bw       = rms_hist_bw    # rms histogram bin width
#  self.nlloc_cf.write('VGGRID  2 %s %s %s %s %s %s SLOW_LEN\n'%(self.vggrid_num_g_n_xy,
#                                                                       self.vggrid_num_g_n_z,
#                                                                       self.vggrid_grid_x,
#                                                                       self.vggrid_grid_y,
#                                                                       self.vggrid_grid_z,
#                                                                       self.vggrid_g_s_xyz))
class NLLocObj():
    def __init__(self,folder_inp,
        trans = "SIMPLE",
        trans_args = [5,-73,0],
        vggrid_num_nodes = [2,583,70],
        vggrid_coord_orig = [-891,-891,-5],
        vggrid_node_spacing =  [3,0,0],
        locsig            = "user"            ,   # signature.
        loccom            = "comment"         ,   # comments.
        locfiles_typ      = "SEISAN"         ,   # obs type: (nlloc_obs, hypo71, hypoellipse, renass_dep, seisan, hypodd_).
        lochypout         =  4              ,   # output type: (1:save_nlloc_all,2:save_nlloc_sum,3:save_hypo71_sum,4:all)  
        locsearch_min_xyz =  [10, 10, 4]        ,   # search type, initial number of cells in x/y/z directions. 
        locsearch_node    =  [0.01, 10000, 5000],   # search type, minimum node size, max number of nodes, number of scatters.
        locsearch_stp     =  [0, 1]            ,   #  1:use stations density or not  0 , 1:stop on min node size or not  0  .
        locgrid_num_nodes  = [2, 1001, 201]      ,   # number of nodes along x/y/z axis.
        locgrid_coord_orig = [-891,-891,-5],
        locgrid_node_spacing   = [1.0, 1.0, 1.0]     ,   # grid spacing along x/y/z axis.
        locmeth_args = [200,4,-1,-1,1.78,6,-1,0],
        locgau            =  [0.2, 0.0]      ,   # gaussian model error parameters (locgau sigma_t (s), corrlen (km)).
        locgau2           =  [0.01, 0.05, 2.0]  ,   # travel-time dependent gaussian model error parameters (%tt,tmax,tmin).
        locphaseid_p      =  ["P","p","PN","PG","Pn","Pg"],   # p-phase identifier mapping.
        locphaseid_s      =  ["S","s","SN","SG","Sn","Sg"],   # s-phase identifier mapping. to ignore s: ==>  $ .
        locphstat_args = [9999,-1,9999,1,1,9999,-9999,9999],

    ) -> None:
        self.folder_inp = folder_inp
        self.vggrid_vel_inp    = os.path.join(folder_inp,"model.dat")        # velocity model file.
        self.vggrid_sta_inp    = os.path.join(folder_inp,"station.dat")          # station file.
        self.locfiles_obs      = os.path.join(folder_inp,"select.out")        # observation file name  could be nordic. nlloc, ... .
        
        self.trans = trans
        self.trans_args = " ".join(list(map(lambda x: str(x),trans_args)))
        
        self.vggrid_num_nodes = " ".join(list(map(lambda x: str(x),vggrid_num_nodes)))
        self.vggrid_coord_orig = " ".join(list(map(lambda x: str(x),vggrid_coord_orig)))
        self.vggrid_node_spacing =  " ".join(list(map(lambda x: str(x),vggrid_node_spacing)))
        
        self.locsig            = locsig                # signature.
        self.loccom            = loccom                # comments.
        self.locfiles_typ      = locfiles_typ          # obs type: (nlloc_obs, hypo71, hypoellipse, renass_dep, seisan, hypodd_).
        self.lochypout         = lochypout             # output type: (1:save_nlloc_all,2:save_nlloc_sum,3:save_hypo71_sum,4:all)  
        self.locsearch_min_xyz = " ".join(list(map(lambda x: str(x),locsearch_min_xyz)))     # search type, initial number of cells in x/y/z directions. 
        self.locsearch_node    = " ".join(list(map(lambda x: str(x),locsearch_node)))        # search type, minimum node size, max number of nodes, number of scatters.
        self.locsearch_stp     = " ".join(list(map(lambda x: str(x),locsearch_stp)))         #  1:use stations density or not  0 , 1:stop on min node size or not  0  .
        self.locgrid_num_nodes  = " ".join(list(map(lambda x: str(x),locgrid_num_nodes)))      # number of nodes along x/y/z axis.
        self.locgrid_coord_orig = " ".join(list(map(lambda x: str(x),locgrid_coord_orig)))
        self.locgrid_node_spacing   = " ".join(list(map(lambda x: str(x),locgrid_node_spacing)))       # grid spacing along x/y/z axis.
        self.locmeth_args = " ".join(list(map(lambda x: str(x),locmeth_args))) 
        self.locgau            = " ".join(list(map(lambda x: str(x),locgau)))                # gaussian model error parameters (locgau sigma_t (s), corrlen (km)).
        self.locgau2           = " ".join(list(map(lambda x: str(x),locgau2)))               # travel-time dependent gaussian model error parameters (%tt,tmax,tmin).
        self.locphaseid_p      = " ".join(list(map(lambda x: str(x),locphaseid_p)))          # p-phase identifier mapping.
        self.locphaseid_s      = " ".join(list(map(lambda x: str(x),locphaseid_s)))          # s-phase identifier mapping. to ignore s: ==>  $ .
        self.locphstat_args = locphstat_args     # max hypocenter rms to include in ave residual.

        self.NLLOC_path = os.path.join(CORE_NLLOC,"NonLinLoc-main")
        self.src_path = os.path.join(self.NLLOC_path,"src")
        self.bin_path = os.path.join(self.src_path,"bin")
        self.nll_exe_path = os.path.join(self.bin_path,"NLLoc")
        self.loc_folder = os.path.join(self.folder_inp,"loc","SeisMonitor")
        self.model_folder = os.path.join(self.folder_inp,"model","layer")
        self.time_folder = os.path.join(self.folder_inp,"time","layer")

    def write_nlloc_cf(self, out='./nlloc.cf',
                        P_flag=True, 
                        S_flag=False,
                        rm_outs=False):

        """
        taken from https://github.com/saeedsltm/PyNLLRunner
        """

        self.nlloc_cf = open(out, 'w')

        self.nlloc_cf.write('#__________________START GENERIC CONTROL STATEMENTS\n\n')
        self.nlloc_cf.write('CONTROL 1 54321\n')
        self.nlloc_cf.write('TRANS  %s %s \n\n'%(self.trans,
                                            self.trans_args))
        self.nlloc_cf.write('#__________________END\n')
        self.nlloc_cf.write('#__________________START VEL2GRID STATEMENTS\n\n')
        self.nlloc_cf.write(f'VGOUT  {self.model_folder}\n')
        
        if P_flag: 
            self.nlloc_cf.write('VGTYPE P\n')
        if S_flag: 
            self.nlloc_cf.write('VGTYPE S\n')

        self.nlloc_cf.write('VGGRID  %s %s %s SLOW_LEN\n'%(self.vggrid_num_nodes,
                                                                self.vggrid_coord_orig,
                                                                    self.vggrid_node_spacing))
        self.nlloc_cf.write('INCLUDE %s\n\n'%(self.vggrid_vel_inp))
        self.nlloc_cf.write('#__________________END\n')
        self.nlloc_cf.write('#__________________START GRID2TIME STATEMENTS\n\n')

        if P_flag:
            self.nlloc_cf.write(f'GTFILES  {self.model_folder} {self.time_folder} P\n')
        
        if S_flag:
            self.nlloc_cf.write(f'GTFILES  {self.model_folder} {self.time_folder} S\n')

        self.nlloc_cf.write('GTMODE GRID2D ANGLES_YES\n')
        self.nlloc_cf.write('INCLUDE %s\n'%(self.vggrid_sta_inp))
        self.nlloc_cf.write('GT_PLFD  1.0e-3  0\n\n')
        self.nlloc_cf.write('#__________________END\n')


        self.nlloc_cf.write('#__________________START NLDIFFLOC STATEMENTS\n\n')
        self.nlloc_cf.write('LOCSIG %s \n'%(self.locsig))
        self.nlloc_cf.write('LOCCOM %s \n'%(self.loccom))

        self.nlloc_cf.write('LOCFILES %s %s %s %s\n'%(self.locfiles_obs,
                                                    self.locfiles_typ,
                                                    self.time_folder,
                                                    self.loc_folder))
            
        output = ['SAVE_NLLOC_ALL','SAVE_NLLOC_SUM','SAVE_HYPO71_SUM']
        output.append(' '.join(output))

        self.nlloc_cf.write('LOCHYPOUT %s\n'%(output[int(self.lochypout)-1]))
               
        self.nlloc_cf.write('LOCSEARCH OCT %s %s %s \n'%(self.locsearch,
                                                         self.locsearch_node,
                                                         self.locsearch_stp))

        self.nlloc_cf.write('LOCGRID %s %s %s %s %s PROB_DENSITY SAVE\n'%(self.locgrid_gnum_xyz ,
                                                                          self.locgrid_grid_x,
                                                                          self.locgrid_grid_y,
                                                                          self.locgrid_grid_z,
                                                                          self.locgrid_g_s_xyz))
       
        self.nlloc_cf.write('LOCMETH EDT_OT_WT %s %s %s %s %s %s %s %s\n'%(self.locmeth_max_st_d,
                                                                            self.locmeth_min_nm_ph,
                                                                            self.locmeth_max_nm_ph,
                                                                            self.locmeth_min_nm_s,
                                                                            self.locmeth_vp_vs,
                                                                            self.locmeth_max_g,
                                                                            self.locmeth_min_st_d,
                                                                            self.locmeth_dup))



        self.nlloc_cf.write('LOCGAU %s\n'%(self.locgau))
        self.nlloc_cf.write('LLOCGAU2 %s\n'%(self.locgau2))
        self.nlloc_cf.write('LOCPHASEID P %s\n'%(self.locphaseid_p))
        self.nlloc_cf.write('LOCPHASEID S %s\n'%(self.locphaseid_s))
        self.nlloc_cf.write('LOCQUAL2ERR 0.1 0.5 1.0 2.0 99999.9\n')
        self.nlloc_cf.write('LOCPHSTAT %s %s %s %s %s %s %s %s\n'%(self.locphstat_rms_max,
                                                                   self.locphstat_nr_min,
                                                                   self.locphstat_gap_max,
                                                                   self.locphstat_p_rmax,
                                                                   self.locphstat_s_rmax,
                                                                   self.locphstat_el3_max,
                                                                   self.locphstat_d_min,
                                                                   self.locphstat_d_max))


        self.nlloc_cf.write('LOCANGLES ANGLES_YES 5\n')
        self.nlloc_cf.write('LOCMAG ML_HB 1.0 1.110 0.00189\n')

        self.nlloc_cf.write('\n#__________________END')
            
        self.nlloc_cf.close()

        for _folder in [self.loc_folder,self.model_folder,self.time_folder]:
            folder = os.path.dirname(_folder)
            if os.path.isdir(folder):
                if rm_outs:
                    shutil.rmtree(folder)
                    os.makedirs(folder)
                else:
                    pass
            else:
                os.makedirs(folder)


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

def download_nlloc(forced=False):
    name = "nll.zip"
    zip_path = os.path.join(CORE_NLLOC,name)
    NLLOC_path = os.path.join(CORE_NLLOC,"NonLinLoc-main")
    src_path = os.path.join(NLLOC_path,"src")
    bin_path = os.path.join(src_path,"bin")
    nll_exe_path = os.path.join(bin_path,"NLLoc")
    cache_path = os.path.join(src_path,"CMakeCache.txt")
    # "https://github.com/alomax/NonLinLoc/archive/refs/heads/main.zip"

    if not forced:
        if os.path.isfile(nll_exe_path):
            return True
        else:
            pass

    isfile = ut.isfile(zip_path)
    if not isfile:
        os.system(f"wget https://github.com/alomax/NonLinLoc/archive/refs/heads/main.zip -O {zip_path}")

    if not os.path.isdir(NLLOC_path):
        os.system(f"unzip {zip_path} -d {CORE_NLLOC}")

    if os.path.isdir(bin_path):
        os.rmdir(bin_path)
        os.makedirs(bin_path)

    if os.path.isdir(cache_path):
        os.rmdir(cache_path)

    apt_install(["cmake"])

    os.system(f"cd {src_path} && cmake . && make")

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
    # print(df)
    return df 

def write_station_file(sta_path,out):
    df = resp2df(sta_path)
    print(df)

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



if __name__ == "__main__":
    # from obspy.core.event import read_events
    # cat_path = "/home/emmanuel/EDCT/test/associations/associations.xml"
    # out="/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test/select.out"
    # cat = read_events(cat_path)
    # cat.write(out,format="NORDIC")

    # download_nlloc()

    folder_inp = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3"
    out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/nlloc.cf"
    nll_obj = NLLocObj(folder_inp=folder_inp)
    # nll_obj.write_nlloc_cf(out=out,P_flag=True, S_flag=False,rm_outs=True)
    nll_obj.write_nlloc_cf(out=out,P_flag=False, S_flag=True,rm_outs=False)

    # vel_model = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/vel1d_col.csv"
    # out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/model.dat"
    # write_1d_vel_model(vel_model,out)

    # vel_model = "/home/emmanuel/EDCT/SeisMonitor/data/metadata/CM.xml"
    # out = "/home/emmanuel/EDCT/SeisMonitor/SeisMonitor/monitor/locator/nlloc/test3/station.dat"
    # write_station_file(vel_model,out)


# https://github.com/amaggi/waveloc/blob/master/nll/nlloc_sample.in