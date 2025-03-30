import os
from obspy.core.event.catalog import Catalog, read_events
from SeisMonitor.utils import printlog, isfile
from . import utils as ut

class Hypocenter:
    """
    A class for handling hypocenter localization using SEISAN format.
    
    Attributes:
        catalog (Catalog): ObsPy Catalog object containing earthquake events.
        xml_path (str): Path to the station metadata XML file.
        vel_path (str): Path to the velocity model CSV file.
        out_dir (str): Directory where output files will be stored.
        out_file (str): Path to the final relocated earthquake catalog.
        paths (str): Directory for Hypocenter-related intermediate files.
        sta0 (str): Path to the STATION0.HYP file used by SEISAN.
    """
    def __init__(self, catalog: str, xml_path: str, vel_path: str, out_dir: str):
        """
        Initializes the Hypocenter class with the necessary input files and directories.
        
        Args:
            catalog (str or Catalog): Path to an earthquake catalog file or an ObsPy Catalog object.
            xml_path (str): Path to the station metadata XML file.
            vel_path (str): Path to the CSV file containing the velocity model.
            out_dir (str): Directory where output files will be stored.
        """
        if isinstance(catalog, Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)
        
        self.xml_path = xml_path
        self.vel_path = vel_path
        self.out_dir = out_dir
        
        # Define paths for output files
        self.out_file = os.path.join(self.out_dir, "hypocenter", "hypocenter_events.xml")
        self.paths = os.path.join(self.out_dir, "hypocenter", "hypocenter_paths")
        self.sta0 = os.path.join(self.paths, "STATION0.HYP")

    def locate(self, sfilename: str = "catalog.sfile", out_format: str = "QUAKEML", rm_not_locatable: bool = True) -> Catalog:
        """
        Runs the hypocenter localization process.
        
        Args:
            sfilename (str, optional): Name of the S-file to be used. Defaults to "catalog.sfile".
            out_format (str, optional): Output format for the relocated catalog. Defaults to "QUAKEML".
            rm_not_locatable (bool, optional): Whether to remove events that cannot be relocated. Defaults to True.
        
        Returns:
            Catalog: The relocated earthquake catalog.
        """
        # Define the path for the S-file
        sfile = os.path.join(self.paths, sfilename)
        isfile(sfile)
        
        # Write the event catalog to S-file format
        self.catalog.write(sfile, format="NORDIC")
        
        # Generate the STATION0 file
        sta0 = ut.STATION0(self.xml_path, self.vel_path)
        sta0.write(self.sta0)
        
        # Run hypocenter tools
        hyp = ut.HypocenterTools(self.paths, sfilename)
        hyp.split()
        hyp.remodl_and_setbrn()
        hyp.update()
        
        # Check the integrity of the S-files
        ut.check_sfile_integrity(self.paths, rm_not_locatable)
        
        # Collect relocated events into a new catalog
        catalog = hyp.collect()
        
        # Write the output catalog if an output file path is specified
        if self.out_file is not None:
            print("Writing output file...")
            isfile(self.out_file)
            catalog.write(self.out_file, out_format)
        
        return catalog
