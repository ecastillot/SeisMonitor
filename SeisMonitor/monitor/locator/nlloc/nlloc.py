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
from obspy.core.event.resourceid import ResourceIdentifier
from obspy import UTCDateTime
import subprocess
from tqdm import tqdm
from . import utils as ut
from SeisMonitor import utils as sut
from SeisMonitor.core import utils as scut
from SeisMonitor.monitor.locator import utils as slut


class NLLoc:
    """NonLinLoc seismic event locator class.

    This class implements the NonLinLoc algorithm for locating seismic events
    using velocity models and station information.
    """

    def __init__(
        self,
        core_path: str,
        agency: str,
        region: list,
        vel_model: slut.VelModel,
        stations: slut.Stations,
        delta_in_km: float = 2,
        kwargs_for_trans: dict = {},
        kwargs_for_vel2grid: dict = {},
        kwargs_for_grid2time: dict = {},
        kwargs_for_time2loc: dict = {},
        tmp_folder: str = os.getcwd(),
        exhaustively: bool = False,
        search_in_degrees: list = [],
        rm_attempts: bool = False,
    ):
        """Initialize NLLoc locator instance.

        Args:
            core_path: Path to NonLinLoc core directory
            agency: Agency identifier string
            region: List of [lon_w, lon_e, lat_s, lat_n, z_min, z_max]
            vel_model: Velocity model object
            stations: Station information object
            delta_in_km: Grid spacing in kilometers
            kwargs_for_trans: Transformation parameters
            kwargs_for_vel2grid: Vel2Grid parameters
            kwargs_for_grid2time: Grid2Time parameters
            kwargs_for_time2loc: Time2Loc parameters
            tmp_folder: Temporary working directory
            exhaustively: Whether to perform exhaustive search
            search_in_degrees: Degrees for exhaustive search
            rm_attempts: Remove temporary attempt files
        """
        paths = ut.testing_nlloc_core_path(core_path)

        self.core_path = core_path
        self.nlloc_paths = paths
        self.agency = agency
        self.region = region
        self.vel_model = vel_model
        self.stations = stations
        self.basic_inputs = slut.LocatorBasicInputs(
            vel_model=vel_model,
            stations=stations
        )
        self.delta_in_km = delta_in_km
        self.kwargs_for_trans = kwargs_for_trans
        self.kwargs_for_vel2grid = kwargs_for_vel2grid
        self.kwargs_for_grid2time = kwargs_for_grid2time
        self.kwargs_for_time2loc = kwargs_for_time2loc
        self.tmp_folder = tmp_folder
        self.rm_attempts = rm_attempts

        self.exhaustively = exhaustively
        if exhaustively:
            if not search_in_degrees:
                self.search_in_degrees = list(reversed(np.arange(1, 4, 0.5)))
            else:
                self.search_in_degrees = search_in_degrees

    def __initialize(self, write_nlloc_files: bool = True):
        """Prepare input files and parameters for location process.

        Args:
            write_nlloc_files: Whether to write NLLoc control files
        """
        # Set up file paths
        self.vel_model_path = os.path.join(self.tmp_folder, "vel_model.dat")
        self.station_path = os.path.join(self.tmp_folder, "station.dat")
        
        # Create/overwrite input files
        sut.isfile(self.vel_model_path, overwrite=True)
        self.basic_inputs.vel_model.to_nlloc(self.vel_model_path)
        sut.isfile(self.station_path, overwrite=True)
        self.basic_inputs.stations.to_nlloc(self.station_path)

        # Define output directories
        self.grid_folder_out = os.path.join(self.tmp_folder, "model", "layer")
        self.time_folder_out = os.path.join(self.tmp_folder, "time", "layer")
        self.loc_folder_out = os.path.join(self.tmp_folder, "loc", "SeisMonitor")
        self.p_control_file_out = os.path.join(self.tmp_folder, "p_nlloc.in")
        self.s_control_file_out = os.path.join(self.tmp_folder, "s_nlloc.in")

        # Prepare grid parameters
        grid_args = self._prepare_grid_args()

        # Update grid arguments with kwargs if provided
        grid_args.update({
            "trans": self.kwargs_for_trans.get("trans", grid_args["trans"]),
            "velgrid": self.kwargs_for_vel2grid.get("grid", grid_args["velgrid"]),
            "locgrid": self.kwargs_for_grid2time.get("grid", grid_args["locgrid"])
        })

        # Create control objects
        gen_control = ut.GenericControlStatement(trans=grid_args["trans"])
        gen_vel2grid = ut.Vel2Grid(
            vel_path=self.vel_model_path,
            grid_folder_out=self.grid_folder_out,
            grid=grid_args["velgrid"]
        )
        p_grid2time = ut.Grid2Time(
            station_path=self.station_path,
            grid_folder_out=self.grid_folder_out,
            time_folder_out=self.time_folder_out,
            phase="P"
        )
        s_grid2time = ut.Grid2Time(
            station_path=self.station_path,
            grid_folder_out=self.grid_folder_out,
            time_folder_out=self.time_folder_out,
            phase="S"
        )

        # Prepare catalog output
        self.catalog_path = os.path.join(self.tmp_folder, "catalog.out")
        sut.isfile(self.catalog_path, overwrite=True)
        open(self.catalog_path, "w").close()
        gen_time2loc = ut.Time2Loc(
            catalog=[self.catalog_path, "SEISAN"],
            grid=grid_args["locgrid"],
            time_folder_out=self.time_folder_out,
            loc_folder_out=self.loc_folder_out
        )

        def write_each_nlloc_file(nlloc_and_out):
            """Helper function to write NLLoc control files."""
            nlloc, out = nlloc_and_out
            nlloc.write(
                out,
                vel2grid=True,
                grid2time=True,
                time2loc=True
            )
            return True

        # Create control files for P and S phases
        p_nlloc = ut.NLLocControlFile(gen_control, gen_vel2grid, p_grid2time, gen_time2loc)
        s_nlloc = ut.NLLocControlFile(gen_control, gen_vel2grid, s_grid2time, gen_time2loc)
        nlloc_control_files = [
            (p_nlloc, self.p_control_file_out),
            (s_nlloc, self.s_control_file_out)
        ]

        if write_nlloc_files:
            # Write control files using single-threaded executor
            with cf.ThreadPoolExecutor(max_workers=1) as executor:
                executor.map(write_each_nlloc_file, nlloc_control_files)

        self.nll_control_file = p_nlloc

    def _prepare_grid_args(self) -> dict:
        """Calculate grid parameters based on region and delta.

        Returns:
            Dictionary containing transformation, velocity grid, and location grid parameters
        """
        lon_w, lon_e, lat_s, lat_n, z_min, z_max = self.region

        # Calculate center point
        c_lat = (lat_s + lat_n) / 2
        c_lon = (lon_w + lon_e) / 2

        # Calculate grid dimensions in km
        x, _, _ = gps2dist_azimuth(lat_s, lon_w, lat_s, c_lon)
        x = x / 1e3  # Convert to km
        x_num = int(x * 2 / self.delta_in_km)
        
        y, _, _ = gps2dist_azimuth(lat_s, lon_w, c_lat, lon_w)
        y = y / 1e3  # Convert to km
        y_num = int(y * 2 / self.delta_in_km)

        z_num = int((z_max - z_min) / self.delta_in_km)

        return {
            "trans": ["SIMPLE", c_lat, c_lon, 0],
            "velgrid": [
                x_num, y_num, z_num,
                -round(x, 2), -round(y, 2), round(z_min, 2),
                self.delta_in_km, self.delta_in_km, self.delta_in_km,
                "SLOW_LEN"
            ],
            "locgrid": [
                x_num, y_num, z_num,
                -round(x, 2), -round(y, 2), round(z_min, 2),
                self.delta_in_km, self.delta_in_km, self.delta_in_km,
                "PROB_DENSITY", "SAVE"
            ]
        }

    def compute_travel_times(self):
        """Compute travel times using velocity model and station locations."""
        self.__initialize()
        
        sut.printlog("info", "NLLoc:Vel2Grid", "Running")
        subprocess.call(
            f"{self.nlloc_paths['vel2grid_exe_path']} {self.p_control_file_out}",
            shell=True
        )
        
        sut.printlog("info", "NLLoc:Grid2Time:P", "Running")
        subprocess.call(
            f"{self.nlloc_paths['grid2time_exe_path']} {self.p_control_file_out}",
            shell=True
        )
        
        sut.printlog("info", "NLLoc:Grid2Time:S", "Running")
        subprocess.call(
            f"{self.nlloc_paths['grid2time_exe_path']} {self.s_control_file_out}",
            shell=True
        )

    def _locate(
        self,
        catalog: Union[Catalog, str],
        nlloc_out_folder: str,
        out_filename: str = "locations.xml",
        out_format: str = "SC3ML"
    ) -> Catalog:
        """Perform single-pass event location.

        Args:
            catalog: Input catalog or path to catalog file
            nlloc_out_folder: Output directory
            out_filename: Output filename
            out_format: Output format

        Returns:
            Located event catalog
        """
        catalog = catalog if isinstance(catalog, Catalog) else read_events(catalog)

        # Set up file paths
        nlloc_inp = os.path.join(nlloc_out_folder, "catalog_input.inp")
        nlloc_folder = os.path.join(nlloc_out_folder, "nlloc", "SeisMonitor")
        nlloc_out = os.path.join(nlloc_out_folder, out_filename)
        nlloc_control = os.path.join(nlloc_out_folder, "loc.in")

        # Write input catalog
        sut.isfile(nlloc_inp, overwrite=True)
        catalog.write(nlloc_inp, format="NORDIC")

        picks_from_unproc_catalog = slut.get_picks(catalog)

        # Initialize if needed and update control file
        try:
            nlloc_control_file = self.nll_control_file
        except AttributeError:
            self.__initialize(False)
            nlloc_control_file = self.nll_control_file

        nlloc_control_file.time2loc.catalog = " ".join((nlloc_inp, "SEISAN"))
        nlloc_control_file.time2loc.loc_folder_out = nlloc_folder
        nlloc_control_file.write(nlloc_control)

        # Run location
        sut.printlog("info", "NLLoc:NLLoc", "Running")
        subprocess.call(
            f"{self.nlloc_paths['nll_exe_path']} {nlloc_control}",
            shell=True
        )

        # Process results
        all_events = []
        for path in tqdm(glob.glob(nlloc_folder + "*.hyp")):
            file_base = os.path.basename(path)
            date = file_base.split(".")[1]
            if date in ["sum", "last.hyp"]:
                continue
                
            try:
                catalog = read_nlloc_hyp(path, format="NORDIC")
            except Exception:
                print(f"Unread: {path}")
                continue

            for ev in catalog.events:
                ori_pref = ev.preferred_origin()
                ori_pref = scut.add_aditional_origin_info(
                    ori_pref,
                    agency=self.agency,
                    method_id="NLLOC",
                    earth_model_id=self.vel_model.model_name
                )
                ev.preferred_origin_id = ori_pref.resource_id.id
                ev = slut.changing_picks_info(ev, picks_from_unproc_catalog)
                ev = scut.add_aditional_event_info(ev, agency=self.agency)
                all_events.append(ev)

        catalog = Catalog(events=all_events)
        catalog = scut.add_aditional_catalog_info(catalog, agency=self.agency)
        
        sut.isfile(nlloc_out)
        catalog.write(nlloc_out, format=out_format)
        sut.printlog("info", "NLLoc:NLLoc", f"Finished. See your results in {nlloc_out}")
        
        return catalog

    def _iterlocate(
        self,
        catalog: Union[Catalog, str],
        nlloc_out_folder: str,
        out_filename: str = "locations.xml",
        out_format: str = "SC3ML",
        degrees: list = [4, 3, 2.5, 2, 1.5, 1],
        rm_attempts: bool = False
    ) -> Catalog:
        """Perform iterative event location with decreasing distance thresholds.

        Args:
            catalog: Input catalog or path
            nlloc_out_folder: Output directory
            out_filename: Output filename
            out_format: Output format
            degrees: List of distance thresholds
            rm_attempts: Remove temporary files

        Returns:
            Located event catalog
        """
        nlloc_out = os.path.join(nlloc_out_folder, out_filename)
        reloc_catalog = self._locate(
            catalog,
            nlloc_out_folder,
            out_filename="base.xml",
            out_format="SC3ML"
        )

        good_evs = []
        bad_catalog = reloc_catalog
        iter_count = 0
        
        while True:
            good_events, bad_events = slut.get_bad_and_good_events(bad_catalog)
            good_evs.append(good_events)
            degree = degrees[iter_count]
            
            if not bad_events:
                print(
                    f"convergence in {iter_count} iterations, "
                    f"only stations with distance < {degrees[iter_count-1]} degrees"
                    if iter_count > 0
                    else "convergence without iterate"
                )
                break
                
            bad_events = slut.filter_arrivals_by_distance(bad_events, degree)
            if not bad_events:
                break
                
            bad_catalog = Catalog(bad_events)
            tmp_path = os.path.join(nlloc_out_folder, "tmp", f"deg_{degree}")
            bad_catalog = self._locate(
                bad_catalog,
                tmp_path,
                out_filename=f"{iter_count}.xml",
                out_format="SC3ML"
            )
            iter_count += 1

        good_evs = [y for x in good_evs for y in x]
        if good_evs:
            reloc_catalog.events = good_evs

        sut.isfile(nlloc_out)
        reloc_catalog.write(nlloc_out, format=out_format)
        sut.printlog("info", "NLLoc:NLLoc", f"Finished. See your results in {nlloc_out}")
        
        if rm_attempts:
            shutil.rmtree(os.path.join(nlloc_out_folder, "tmp"))

        return reloc_catalog

    def locate(
        self,
        catalog: Union[Catalog, str],
        nlloc_out_folder: str,
        out_filename: str = "locations.xml",
        out_format: str = "SC3ML"
    ) -> Catalog:
        """Main location method with optional exhaustive search.

        Args:
            catalog: Input catalog or path
            nlloc_out_folder: Output directory
            out_filename: Output filename
            out_format: Output format

        Returns:
            Located event catalog
        """
        if self.exhaustively:
            return self._iterlocate(
                catalog=catalog,
                nlloc_out_folder=nlloc_out_folder,
                out_filename=out_filename,
                out_format=out_format,
                degrees=self.search_in_degrees,
                rm_attempts=self.rm_attempts
            )
        return self._locate(
            catalog=catalog,
            nlloc_out_folder=nlloc_out_folder,
            out_filename=out_filename,
            out_format=out_format
        )