import warnings
import os
import glob
import shutil
from SeisMonitor.monitor.downloader.utils import get_chunktimes
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader
from SeisMonitor.monitor.downloader.utils import sanitize_provider_times
from SeisMonitor.monitor.picker import ai as ai_picker
from SeisMonitor.monitor.associator import ai as ai_asso
from SeisMonitor.monitor.locator.nlloc import nlloc
from SeisMonitor.monitor.magnitude.mag import Magnitude

# Suppress FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_preproc_providers(providers, chunklength_in_sec, out_folder):
    """Generate preprocessed providers split by time chunks.

    Args:
        providers (list): List of provider objects with waveform restrictions
        chunklength_in_sec (int): Length of each time chunk in seconds
        out_folder (str): Base output directory for processed data

    Returns:
        list: List of dictionaries containing providers and folder structure
              for each time chunk
    """
    sanitize_provider_times(providers)
    oneprovider = providers[0]
    restrictions = oneprovider.waveform_restrictions
    starttime = restrictions.starttime
    endtime = restrictions.endtime
    
    chunktimes = get_chunktimes(starttime, endtime, chunklength_in_sec)
    preproc_providers = []
    
    for starttime, endtime in chunktimes:
        folders = get_folders_by_chunk(out_folder, starttime, endtime)
        new_providers = []
        
        for provider in providers:
            provider.waveform_restrictions.starttime = starttime
            provider.waveform_restrictions.endtime = endtime
            new_providers.append(provider)
        
        chunk_provider = {
            "providers": new_providers,
            "folders": folders,
            "dates": {"starttime": starttime, "endtime": endtime}
        }
        preproc_providers.append(chunk_provider)
    
    return preproc_providers


def get_chunktimes_by_provider(providers, chunklength_in_sec):
    """Get time chunks for each provider.

    Args:
        providers (list): List of provider objects
        chunklength_in_sec (int): Length of each time chunk in seconds

    Returns:
        dict: Dictionary mapping providers to their time chunks

    Note:
        This function appears incomplete as it doesn't populate the dictionary.
    """
    times_by_provider = {}
    for provider in providers:
        restrictions = provider.waverform_restrictions  # Note: typo in original 'waverform'
        starttime = restrictions.starttime
        endtime = restrictions.endtime
        times = get_chunktimes(starttime, endtime, chunklength_in_sec)
        times_by_provider.append(times)  # Should be times_by_provider[provider] = times?
    return times_by_provider


def get_folders_by_chunk(out_folder, starttime, endtime):
    """Generate folder structure for a time chunk.

    Args:
        out_folder (str): Base output directory
        starttime (UTCDateTime): Start time of the chunk
        endtime (UTCDateTime): End time of the chunk

    Returns:
        dict: Dictionary of folder paths for different processing stages
    """
    st = starttime.strftime("%Y%m%dT%H%M%S")
    et = endtime.strftime("%Y%m%dT%H%M%S")
    chunk_name = f"{st}__{et}"
    chunk_dir = os.path.join(out_folder, chunk_name)

    return {
        "metadata": os.path.join(chunk_dir, 'metadata'),
        "downloads": os.path.join(chunk_dir, 'downloads'),
        "detections": os.path.join(chunk_dir, 'detections'),
        "associations": os.path.join(chunk_dir, 'associations'),
        "locations": os.path.join(chunk_dir, 'locations'),
        "magnitudes": os.path.join(chunk_dir, 'magnitudes')
    }


def sanitize_pick_batch_size(pickers, download_args):
    """Adjust download arguments based on picker parameters.

    Args:
        pickers (dict): Dictionary of picker objects and their arguments
        download_args (dict): Download configuration arguments

    Returns:
        dict: Updated download arguments with minimum overlap and batch size
    """
    overlaps = [args.overlap for args in pickers.values()]
    batch_sizes = [args.batch_size for args in pickers.values()]
    
    download_args["picker_args"]["overlap"] = min(overlaps)
    download_args["picker_args"]["batch_size"] = min(batch_sizes)
    return download_args


def sanitize_downloads(pickers):
    """Modify picker arguments to preserve downloads until last picker.

    Args:
        pickers (dict): Dictionary of picker objects and their arguments

    Returns:
        dict: Updated pickers dictionary with rm_download set to False except for last
    """
    new_pickers = {}
    for i, (picker, args) in enumerate(pickers.items()):
        if i != len(pickers) - 1:
            args.rm_download = False
        new_pickers[picker] = args
    return new_pickers


class SeisMonitor:
    """Main class for seismic monitoring pipeline."""

    def __init__(self, providers, out_folder, chunklength_in_sec=3600, overwrite=False):
        """Initialize SeisMonitor with basic configuration.

        Args:
            providers (list): List of data provider objects
            out_folder (str): Base output directory
            chunklength_in_sec (int): Length of time chunks in seconds
            overwrite (bool): Whether to overwrite existing files
        """
        self.providers = providers
        self.out_folder = out_folder
        self.chunklength_in_sec = chunklength_in_sec
        self.overwrite = overwrite
        self.process = {}

    def add_downloader(
        self,
        threshold=60,
        overlap_in_sec=0,
        picker_args={},
        groupby='{network}.{station}.{channel}',
        n_processor=None
    ):
        """Add downloader process to the pipeline.

        Args:
            threshold (int): Minimum length in seconds for download
            overlap_in_sec (int): Overlap between chunks in seconds
            picker_args (dict): Arguments for picker compatibility
            groupby (str): Pattern for grouping traces
            n_processor (int, optional): Number of processors to use
        """
        dld_args = locals().copy()
        dld_args["chunklength_in_sec"] = self.chunklength_in_sec
        dld_args.pop("self")
        self.process["downloader"] = dld_args

    def add_picker(self, pickers={}):
        """Add picker processes to the pipeline.

        Args:
            pickers (dict): Dictionary of picker names and their argument objects

        Returns:
            list: List of picker names added
        """
        if pickers:
            pickers = sanitize_downloads(pickers)
            self.process["picker"] = pickers
            if "downloader" in self.process:
                self.process["downloader"] = sanitize_pick_batch_size(
                    pickers, self.process["downloader"]
                )
        
        self.picker_output = list(pickers.keys())
        return self.picker_output

    def add_associator(self, input, associators={}):
        """Add associator processes to the pipeline.

        Args:
            input (list): List of picker names to associate
            associators (dict): Dictionary of associator names and their arguments

        Returns:
            dict: Dictionary of association output paths
        """
        self.associator_input = input
        if associators:
            self.process["associator"] = associators

        out = {}
        for associator in associators:
            for picker in self.associator_input:
                name = f"{associator}_{picker}"
                out[name] = os.path.join(associator, picker)
        self.associator_output = out
        return out

    def add_locator(self, input, locators={}):
        """Add locator processes to the pipeline.

        Args:
            input (dict): Dictionary of task types and project tuples
            locators (dict): Dictionary of locator names and their arguments

        Returns:
            dict: Dictionary of location output paths
        """
        self.locator_input = input
        if locators:
            self.process["locator"] = locators

        out = {}
        for locator in locators:
            for task, project in self.locator_input.items():
                assert task in ["associations", "locations", "magnitudes"]
                assert isinstance(project, tuple) and len(project) == 2
                out_name = f"{locator}_{task}_{project[0]}_{project[1]}"
                out[out_name] = os.path.join(locator, task, *project)
        
        self.locator_output = out
        return out

    def add_magnitude(self, input, magnitudes={}):
        """Add magnitude calculation processes to the pipeline.

        Args:
            input (dict): Dictionary of task types and project tuples
            magnitudes (dict): Dictionary of magnitude methods and their arguments

        Returns:
            dict: Dictionary of magnitude output paths
        """
        self.magnitude_input = input
        if magnitudes:
            self.process["magnitude"] = magnitudes

        out = {}
        for magnitude in magnitudes:
            for task, project in self.magnitude_input.items():
                assert task in ["associations", "locations", "magnitudes"]
                assert isinstance(project, tuple) and len(project) == 2
                out_name = f"{magnitude}_{task}_{project[0]}_{project[1]}"
                out[out_name] = os.path.join(magnitude, task, *project)
        
        self.magnitude_output = out
        return out

    def run(self):
        """Execute the configured seismic monitoring pipeline."""
        preproc_providers = get_preproc_providers(
            self.providers, self.chunklength_in_sec, self.out_folder
        )

        for chunk_provider in preproc_providers:
            print(
                f"chunk: {chunk_provider['dates']['starttime']} "
                f"-- {chunk_provider['dates']['endtime']}"
            )
            providers = chunk_provider["providers"]
            folders = chunk_provider["folders"]

            for process, process_args in self.process.items():
                if process == "downloader":
                    structure = "{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
                    download_path = os.path.join(folders["downloads"], structure)
                    md = MseedDownloader(providers)
                    md.make_inv_and_json(folders["metadata"])
                    md.download(download_path, **process_args)
                    del md

                elif process == "picker":
                    for picker, picker_args in process_args.items():
                        out_path = os.path.join(folders["detections"], picker)
                        if picker == "EQTransformer":
                            _picker = ai_picker.EQTransformer(picker_args)
                        elif picker == "PhaseNet":
                            _picker = ai_picker.PhaseNet(picker_args)
                        else:
                            continue
                        
                        result = _picker.pick(folders["downloads"], folders["metadata"], out_path)
                        if result.empty:
                            print("No picks")
                            continue
                        del _picker
                        del result

                elif process == "associator":
                    inv = os.path.join(folders["metadata"], "inv.xml")
                    for picker in self.associator_input:
                        picker_path = os.path.join(folders["detections"], picker)
                        picks_path = os.path.join(picker_path, "results", "seismonitor_picks.csv")
                        for associator, associator_args in process_args.items():
                            out_name = self.associator_output[f"{associator_args.name}_{picker}"]
                            out_folder = os.path.join(folders["associations"], out_name)
                            if associator == "GaMMA":
                                _associator = ai_asso.GaMMA(associator_args)
                                _, result, _ = _associator.associate(picks_path, inv, out_folder)
                                if result.empty:
                                    print("No associated picks")
                                    continue

                elif process == "locator":
                    for locator, locator_args in process_args.items():
                        for task, project in self.locator_input.items():
                            catalog = os.path.join(folders[task], *project, f"{task}.xml")
                            nlloc_folder = os.path.join(folders["locations"], locator, *project)
                            locator_args.locate(catalog, nlloc_folder)

                elif process == "magnitude":
                    for magnitude, magnitude_args in process_args.items():
                        for task, project in self.magnitude_input.items():
                            catalog = os.path.join(folders[task], *project, f"{task}.xml")
                            mag_folder = os.path.join(folders["magnitudes"], magnitude, *project)
                            mag = Magnitude(providers=self.providers, catalog=catalog, out_dir=mag_folder)
                            if magnitude == "Ml":
                                mag.get_Ml(**magnitude_args)