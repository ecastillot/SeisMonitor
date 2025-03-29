import os
import pandas as pd
from obspy.core.event.catalog import Catalog, read_events
from SeisMonitor.utils import printlog, isfile
from . import utils as ut

class HypoDD:
    """
    A class to handle HypoDD earthquake relocation.
    
    Attributes:
        catalog (Catalog): ObsPy Catalog object containing earthquake events.
        xml_path (str): Path to the station metadata file in XML format.
        vel_path (str): Path to the velocity model file (CSV format).
        out_dir (str): Directory where output files will be stored.
        out_file (str): Path to the final relocated earthquake catalog.
        paths (str): Directory for HypoDD-related intermediate files.
        pha (str): Path to the phase file for HypoDD.
        vel_df (DataFrame): Pandas DataFrame containing velocity model information.
    """
    
    def __init__(self, catalog, xml_path, vel_path, out_dir):
        """
        Initializes the HypoDD class with the necessary input files and directories.
        
        Args:
            catalog (str or Catalog): Path to an earthquake catalog file or an ObsPy Catalog object.
            xml_path (str): Path to the station metadata XML file.
            vel_path (str): Path to the CSV file containing the velocity model.
            out_dir (str): Directory where output files will be stored.
        """
        
        # Load catalog from file if not already a Catalog object
        if isinstance(catalog, Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)
        
        self.xml_path = xml_path
        self.vel_path = vel_path
        self.out_dir = out_dir
        
        # Define paths for output files
        self.out_file = os.path.join(self.out_dir, "hypodd", "hypodd_events.xml")
        self.paths = os.path.join(self.out_dir, "hypodd", "hypodd_paths")
        
        # Create output directory if it doesn't exist
        if not os.path.isdir(self.paths):
            os.makedirs(self.paths)
        
        # Define paths for phase file and velocity model DataFrame
        self.pha = os.path.join(self.paths, "hypoDD.pha")
        self.vel_df = pd.read_csv(vel_path)

    def locate(self, vp_vs_ratio=1.84, out_format="QUAKEML", rm_not_locatable=True):
        """
        Runs the HypoDD relocation process.
        
        Args:
            vp_vs_ratio (float, optional): The Vp/Vs ratio used in velocity modeling. Default is 1.84.
            out_format (str, optional): Output format for the relocated catalog (e.g., "QUAKEML"). Default is "QUAKEML".
            rm_not_locatable (bool, optional): Whether to remove events that cannot be relocated. Default is True.
        """
        
        # Convert XML station metadata to a DataFrame
        df = ut.resp2df(self.xml_path)
        
        # Write station information required for HypoDD
        ut.write_hypoDDstation(df, self.paths)
        
        # Write phase file from event catalog
        ut.write_pha(self.catalog, self.pha)
        
        # Write phase-to-difference time input file for HypoDD
        ut.write_ph2dt_inp_file(self.paths)
        
        # Extract velocity layers from the velocity model
        vel_layers = ut.get_vel_layers(self.vel_df)
        
        # Setup velocity model for HypoDD
        vel_model = ut.setup_velocity_model(
            "layered_p_velocity_with_constant_vp_vs_ratio",
            vp_vs_ratio=vp_vs_ratio,
            layer_tops=vel_layers
        )
        
        # Write HypoDD input file
        ut.write_hypoDD_inp_file(vel_model, self.paths)
        
        # Define paths for HypoDD binaries
        hypodd_root = "/home/emmanuel/QuakeFlow/HypoDD"
        ph2dt_path = os.path.join(hypodd_root, "HYPODD", "src", "ph2dt", "ph2dt")
        hypodd_path = os.path.join(hypodd_root, "HYPODD", "src", "hypoDD", "hypoDD")
        
        # Define paths for input files
        ph2dt_inp_path = os.path.join(self.paths, "ph2dt.inp")
        hypodd_inp_path = os.path.join(self.paths, "hypoDD.inp")
        
        # Define shell commands for running HypoDD steps
        PH2DT_CMD = f"cd {self.paths} && {ph2dt_path} ph2dt.inp"
        HYPODD_CMD = f"cd {self.paths} && {hypodd_path} hypoDD.inp"
        
        # Execute phase-to-difference time conversion
        if os.system(PH2DT_CMD) != 0:
            raise RuntimeError(f"{PH2DT_CMD} failed!")
        
        # Execute HypoDD relocation
        if os.system(HYPODD_CMD) != 0:
            raise RuntimeError(f"{HYPODD_CMD} failed!")
        printlog("HypoDD relocation completed successfully.")