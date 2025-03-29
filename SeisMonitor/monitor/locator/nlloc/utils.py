import os
import shutil
import sys
from time import time
import pandas as pd
from subprocess import run
from obspy import read_inventory
from tqdm import tqdm
import glob
from git import Repo
import SeisMonitor.utils as sut

def get_nlloc_folders(core_path):
    """Generate dictionary of NonLinLoc folder paths.
    
    Args:
        core_path (str): Base path to NonLinLoc core directory
        
    Returns:
        dict: Dictionary containing all relevant NonLinLoc paths
    """
    pre_core_path = os.path.dirname(core_path)
    src_path = os.path.join(core_path, "src")
    bin_path = os.path.join(src_path, "bin")
    vel2grid_exe_path = os.path.join(bin_path, "Vel2Grid")
    grid2time_exe_path = os.path.join(bin_path, "Grid2Time")
    nll_exe_path = os.path.join(bin_path, "NLLoc")

    paths = {
        "pre_core_path": pre_core_path,
        "core_path": core_path,
        "src_path": src_path,
        "bin_path": bin_path,
        "vel2grid_exe_path": vel2grid_exe_path,
        "grid2time_exe_path": grid2time_exe_path,
        "nll_exe_path": nll_exe_path
    }
    return paths

def testing_nlloc_core_path(nlloc_core_path):
    """Validate existence of essential NonLinLoc directories.
    
    Args:
        nlloc_core_path (str): Path to NonLinLoc core directory
        
    Returns:
        dict: Validated paths dictionary
        
    Raises:
        Exception: If mandatory paths are missing
    """
    paths = get_nlloc_folders(nlloc_core_path)
    
    for key, path in paths.items():
        if key in ["pre_core_path", "core_path", "src_path", "bin_path"]:
            if not os.path.isdir(path):
                raise Exception(
                    f"Mandatory path was not found -> {path}.\n"
                    "There is no NLLoc core folder, or it could be corrupted. "
                    "If you are using Ubuntu, feel free to use NLLoc.download()"
                )
    return paths

def run_nlloc(nlloc_paths, p_control_file_path, s_control_file_path):
    """Execute NonLinLoc processing steps using system commands.
    
    Args:
        nlloc_paths (dict): Dictionary of NonLinLoc executable paths
        p_control_file_path (str): Path to P-wave control file
        s_control_file_path (str): Path to S-wave control file
    """
    sut.printlog("info", "NLLoc:Vel2Grid", "Running")
    os.system(f"{nlloc_paths['vel2grid_exe_path']} {p_control_file_path} > /dev/null")
    
    sut.printlog("info", "NLLoc:Grid2Time:P", "Running")
    os.system(f"{nlloc_paths['grid2time_exe_path']} {p_control_file_path} > /dev/null")
    
    sut.printlog("info", "NLLoc:Grid2Time:S", "Running")
    os.system(f"{nlloc_paths['grid2time_exe_path']} {s_control_file_path} > /dev/null")
    
    sut.printlog("info", "NLLoc:NLLoc", "Running")
    os.system(f"{nlloc_paths['nll_exe_path']} {s_control_file_path} > /dev/null")

def apt_install(pkgs):
    """Install packages using apt-get with sudo privileges.
    
    Args:
        pkgs (list): List of package names to install
    """
    cmd = ['pkexec', 'apt-get', 'install', '-y'] + pkgs
    print(f"Running command: {' '.join(cmd)}")
    result = run(
        cmd,
        stdout=sys.stdout,
        stderr=sys.stderr,
        encoding='utf8',
        env={**os.environ, 'DEBIAN_FRONTEND': 'noninteractive'}
    )
    result.check_returncode()

def write_pref_origin_removing_phaselocinfo(catalog):
    """Process seismic catalog to keep only preferred origins.
    
    Args:
        catalog (obspy.Catalog): Input seismic catalog
        
    Returns:
        obspy.Catalog: Modified catalog with cleaned origins
    """
    events = []
    for ev in catalog:
        pref_origin = ev.preferred_origin()
        
        # Standardize pick information
        for pick in ev.picks:
            pick.waveform_id.network_code = "NN"
            pick.waveform_id.location_code = "NN"
            pick.waveform_id.channel_code = "NNN"
            pick.evaluation_mode = "automatic"

        # Clear arrival parameters
        new_arrivals = []
        for arrival in pref_origin.arrivals:
            arrival.azimuth = None
            arrival.distance = None
            arrival.takeoff_angle = None
            arrival.time_residual = None
            arrival.time_weight = None
            new_arrivals.append(arrival)

        # Update origins with cleaned arrivals
        for i, origin in enumerate(ev.origins):
            if origin.resource_id.id == ev.preferred_origin_id:
                ev.origins[i].arrivals = new_arrivals
        
        ev.origins = [pref_origin]
        events.append(ev)
    
    catalog.events = events
    return catalog

def download_nlloc(nlloc_path, forced=False):
    """Download and compile NonLinLoc from Git repository.
    
    Args:
        nlloc_path (str): Destination path for NonLinLoc
        forced (bool): Force re-download even if exists (default: False)
        
    Returns:
        bool: True if already installed and not forced
    """
    nlloc_paths = get_nlloc_folders(nlloc_path)
    zip_path = os.path.join(nlloc_paths['pre_core_path'], "nll.zip")
    cache_path = os.path.join(nlloc_paths['src_path'], "CMakeCache.txt")

    if not forced and os.path.isfile(nlloc_paths['nll_exe_path']):
        return True

    git_url = "https://github.com/alomax/NonLinLoc.git"
    if not os.path.isdir(nlloc_paths['core_path']):
        Repo.clone_from(git_url, nlloc_paths['core_path'])

    if os.path.isdir(nlloc_paths['bin_path']):
        shutil.rmtree(nlloc_paths['bin_path'])
    os.makedirs(nlloc_paths['bin_path'])

    if os.path.isdir(cache_path):
        shutil.rmtree(cache_path)

    try:
        apt_install(["cmake"])
        os.system(f"cd {nlloc_paths['src_path']} && cmake . && make")
    except Exception as e:
        raise Exception(f"Could not compile NLLoc: {str(e)}")

def write_1d_vel_model(vel_path, out, compute_vs=True, vp_vs_ratio=1.78):
    """Write 1D velocity model file from CSV data.
    
    Args:
        vel_path (str): Path to input velocity CSV file
        out (str): Output file path
        compute_vs (bool): Calculate Vs from Vp if True (default: True)
        vp_vs_ratio (float): Vp/Vs ratio for Vs calculation (default: 1.78)
    """
    df = pd.read_csv(vel_path)
    
    with open(out, 'w') as vm:
        vm.write("# model layers (LAYER depth, Vp_top, Vp_grad, Vs_top, Vs_grad, p_top, p_grad)\n")
        for i, row in df.iterrows():
            vs = row.vp / vp_vs_ratio if compute_vs else row.vs
            enter = "" if i == len(df) - 1 else "\n"
            msg = f"LAYER    {row.depth:<6.2f}    {row.vp:<.2f}    0.00    {vs:<.2f}    0.00    {row.rho:<.2f}    0.00{enter}"
            vm.write(msg)

def resp2df(resp):
    """Convert RESP file to DataFrame with station information.
    
    Args:
        resp (str): Path to RESP file
        
    Returns:
        pd.DataFrame: DataFrame with network, station, lat, lon, elevation
    """
    networks, stations, longitudes, latitudes, elevations = [], [], [], [], []
    inv = read_inventory(resp)
    
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
            elevations.append(sta.elevation)
            stations.append(sta.code)
            networks.append(net.code)

    return pd.DataFrame({
        "network": networks,
        "station": stations,
        "latitude": latitudes,
        "longitude": longitudes,
        "elevation": elevations
    })

def write_station_file(sta_path, out):
    """Write station file from RESP data.
    
    Args:
        sta_path (str): Path to RESP file
        out (str): Output file path
    """
    df = resp2df(sta_path)
    
    with open(out, 'w') as vs:
        vs.write("# GTSRCE label LATLON latSrce longSrce zSrce elev\n")
        for i, row in df.iterrows():
            elv = row.elevation / 1e3
            enter = "" if i == len(df) - 1 else "\n"
            msg = f"GTSRCE  {row.station:<5}  LATLON  {row.latitude:<7.3f}  {row.longitude:<8.3f}  0.000  {elv:<7.3f}{enter}"
            vs.write(msg)

def join_args(args):
    """Convert list of arguments to space-separated string.
    
    Args:
        args (list): List of arguments
        
    Returns:
        str: Space-separated string of arguments
    """
    return " ".join(map(str, args))

class GenericControlStatement:
    """Class for generating generic control statements for NonLinLoc."""
    
    def __init__(self, trans, control=[1, 54321]):
        """
        Args:
            trans (list): Transformation parameters [type, lat, lon, depth]
            control (list): Control parameters [level, random_seed] (default: [1, 54321])
        """
        self.control = join_args(control)
        self.trans = join_args(trans)
        sut.validate(self.__init__, locals())

    def get_msg(self):
        """Generate formatted control statement string.
        
        Returns:
            str: Formatted control statement
        """
        msg = (
            "#__________________START GENERIC CONTROL STATEMENTS\n\n"
            "CONTROL {}\n"
            "TRANS {}\n"
            "#__________________END\n"
        )
        sut.printlog("debug", "NLLoc:GenericControlStatement", "Message received")
        return msg.format(self.control, self.trans)

class Vel2Grid:
    """Class for generating velocity grid statements for NonLinLoc."""
    
    def __init__(self, vel_path, grid_folder_out, grid, p_phase="P", s_phase="S"):
        """
        Args:
            vel_path (str): Path to velocity model file
            grid_folder_out (str): Output folder for grid files
            grid (list): Grid parameters [nx, ny, nz, x0, y0, z0, dx, dy, dz, type]
            p_phase (str): P-wave phase identifier (default: "P")
            s_phase (str): S-wave phase identifier (default: "S")
        """
        self.vel_path = vel_path
        self.grid = join_args(grid)
        self.grid_folder_out = grid_folder_out
        self.p_phase = p_phase
        self.s_phase = s_phase
        sut.validate(self.__init__, locals())

    def get_msg(self):
        """Generate formatted velocity grid statement string.
        
        Returns:
            str: Formatted velocity grid statement
        """
        msg = (
            "#__________________START VEL2GRID STATEMENTS\n\n"
            "VGOUT {}\n"
            "VGTYPE {}\n"
            "VGTYPE {}\n"
            "VGGRID {}\n"
            "INCLUDE {}\n"
            "#__________________END\n"
        )
        path = os.path.dirname(self.grid_folder_out)
        sut.isfile(path)
        sut.printlog("debug", "NLLoc:Vel2grid", "Message received")
        return msg.format(
            self.grid_folder_out,
            self.p_phase,
            self.s_phase,
            self.grid,
            self.vel_path
        )

class Grid2Time:
    """Class for generating travel time grid statements for NonLinLoc."""
    
    def __init__(self, station_path, grid_folder_out, time_folder_out, phase="P",
                 mode=["GRID3D", "ANGLES_YES"], plfd=[1e-3, 0]):
        """
        Args:
            station_path (str): Path to station file
            grid_folder_out (str): Input folder with grid files
            time_folder_out (str): Output folder for time files
            phase (str): Phase type (default: "P")
            mode (list): Mode parameters [grid_type, angles] (default: ["GRID3D", "ANGLES_YES"])
            plfd (list): Path length and first difference parameters (default: [1e-3, 0])
        """
        self.station_path = station_path
        self.grid_folder_out = grid_folder_out
        self.time_folder_out = time_folder_out
        self.phase = phase
        self.mode = join_args(mode)
        self.plfd = join_args(plfd)
        sut.validate(self.__init__, locals())

    def get_msg(self):
        """Generate formatted travel time grid statement string.
        
        Returns:
            str: Formatted travel time grid statement
        """
        msg = (
            "#__________________START GRID2TIME STATEMENTS\n\n"
            "GTFILES {} {} {}\n"
            "GTMODE {}\n"
            "INCLUDE {}\n"
            "GT_PLFD {}\n"
            "#__________________END\n"
        )
        path = os.path.dirname(self.time_folder_out)
        sut.isfile(path)
        sut.printlog("debug", "NLLoc:Grid2Time", "Message received")
        return msg.format(
            self.grid_folder_out,
            self.time_folder_out,
            self.phase,
            self.mode,
            self.station_path,
            self.plfd
        )

class Time2Loc:
    """Class for generating location statements for NonLinLoc."""
    
    def __init__(self, catalog, grid, time_folder_out, loc_folder_out,
                 meth=["GAU_ANALYTIC", 9999, 4, -1, -1, 1.78, 6],
                 search=["OCT", 37, 58, 7, 1e-2, int(1e5), int(1e4)],
                 sig="SeisMonitor", com="Comment", gau=[0.2, 0.0], gau2=[0.05, 0.05, 2.0],
                 p_phaseid=["P", "P", "p", "PN", "PG", "Pn", "Pg"],
                 s_phaseid=["S", "S", "s", "SN", "SG", "Sn", "Sg"],
                 qual2err=[0.1, 0.5, 1.0, 2.0, 99999.9],
                 phstat=[9999.0, -1, 9999.0, 1.0, 1.0, 9999.0, -9999.0, 9999.0],
                 angles=["ANGLES_YES", 5],
                 hypout=["SAVE_NLLOC_ALL", "SAVE_NLLOC_SUM", "SAVE_HYPO71_SUM"],
                 mag=["ML_HB", 1.0, 1.110, 0.00189]):
        """
        Args:
            catalog (list): [catalog_path, format]
            grid (list): Grid parameters for location
            time_folder_out (str): Input folder with time files
            loc_folder_out (str): Output folder for location results
            meth (list): Method parameters (default: GAU_ANALYTIC settings)
            search (list): Search parameters (default: OCT settings)
            sig (str): Signature (default: "SeisMonitor")
            com (str): Comment (default: "Comment")
            gau (list): Gaussian parameters (default: [0.2, 0.0])
            gau2 (list): Secondary Gaussian parameters (default: [0.05, 0.05, 2.0])
            p_phaseid (list): P-phase identifiers
            s_phaseid (list): S-phase identifiers
            qual2err (list): Quality to error mapping
            phstat (list): Phase statistics parameters
            angles (list): Angle parameters
            hypout (list): Hypocenter output options
            mag (list): Magnitude calculation parameters
        """
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
        """Generate formatted location statement string.
        
        Returns:
            str: Formatted location statement
        """
        msg = (
            "#__________________START NLDIFFLOC STATEMENTS\n\n"
            "LOCSIG {}\n"
            "LOCCOM {}\n"
            "LOCFILES {} {} {}\n"
            "LOCHYPOUT {}\n"
            "LOCSEARCH {}\n"
            "LOCGRID {}\n"
            "LOCMETH {}\n"
            "LOCGAU {}\n"
            "LOCGAU2 {}\n"
            "LOCPHASEID {}\n"
            "LOCPHASEID {}\n"
            "LOCQUAL2ERR {}\n"
            "LOCPHSTAT {}\n"
            "LOCANGLES {}\n"
            "LOCMAG {}\n"
            "#__________________END\n"
        )
        path = os.path.dirname(self.loc_folder_out)
        sut.isfile(path)
        sut.printlog("debug", "NLLoc:Time2Loc", "Message received")
        return msg.format(
            self.sig,
            self.com,
            self.catalog,
            self.time_folder_out,
            self.loc_folder_out,
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
            self.mag
        )

class NLLocControlFile:
    """Class for generating complete NonLinLoc control file."""
    
    def __init__(self, generic_control, vel2grid, grid2time, time2loc):
        """
        Args:
            generic_control (GenericControlStatement): Generic control instance
            vel2grid (Vel2Grid): Velocity grid instance
            grid2time (Grid2Time): Travel time grid instance
            time2loc (Time2Loc): Location instance
        """
        self.generic_control = generic_control
        self.vel2grid = vel2grid
        self.grid2time = grid2time
        self.time2loc = time2loc
        sut.validate(self.__init__, locals())

    def _validate_args(self, key, input_path, output_path):
        """Validate input and output paths for control file components.
        
        Args:
            key (str): Component identifier
            input_path (str): Input file path
            output_path (str): Output directory path
            
        Raises:
            Exception: If input path doesn't exist
        """
        sut.printlog("debug", "NLLoc:validate_control_file_args", "validating control file args")
        
        if os.path.isfile(input_path):
            sut.printlog("debug", "NLLoc:validate_control_file_args", f"{key} is ok.")
        else:
            raise Exception(f"NLLoc: {key} path doesn't exist. Check: {input_path}")
        
        sut.isfile(output_path)
        sut.printlog("debug", "NLLoc:validate_control_file_args", "control file args are ok")

    def get_msg(self, vel2grid=True, grid2time=True, time2loc=True):
        """Generate complete control file message.
        
        Args:
            vel2grid (bool): Include velocity grid section (default: True)
            grid2time (bool): Include travel time grid section (default: True)
            time2loc (bool): Include location section (default: True)
            
        Returns:
            str: Complete control file content
        """
        msg = self.generic_control.get_msg()

        if vel2grid:
            msg += self.vel2grid.get_msg()
            self._validate_args("vel_path", self.vel2grid.vel_path, self.vel2grid.grid_folder_out)
        if grid2time:
            msg += self.grid2time.get_msg()
            self._validate_args("station_path", self.grid2time.station_path, self.grid2time.time_folder_out)
        if time2loc:
            msg += self.time2loc.get_msg()
            self._validate_args("catalog_path", self.time2loc.catalog.split(" ")[0], self.time2loc.loc_folder_out)

        return msg

    def write(self, out, vel2grid=True, grid2time=True, time2loc=True):
        """Write control file to disk.
        
        Args:
            out (str): Output file path
            vel2grid (bool): Include velocity grid section (default: True)
            grid2time (bool): Include travel time grid section (default: True)
            time2loc (bool): Include location section (default: True)
        """
        sut.printlog("info", "NLLoc:write_control_file", "Running")
        sut.isfile(out)

        msg = self.get_msg(vel2grid, grid2time, time2loc)

        with open(out, "w") as control_file_msg:
            control_file_msg.write(msg)
        
        sut.printlog("info", "NLLoc:write_control_file", f"Finished. Control file: {out}")

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
