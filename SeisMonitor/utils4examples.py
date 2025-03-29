import os
import shutil
from git import Repo
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

def clone_seismonitor_data(output_folder,branch):
    """
    Clones the SeisMonitor repository from GitHub into a specified folder.

    Parameters
    ----------
    output_folder : str
        The path to the folder where the repository will be cloned.
    branch : str
        The branch of the repository to clone.

    Returns
    -------
    bool
        Returns True if the cloning operation was successful.

    Notes
    -----
    If the target folder already exists, it will be removed before cloning the repository.
    """
    git_url = "https://github.com/ecastillot/SeisMonitor.git"

    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder)

    Repo.clone_from(git_url, output_folder,branch=branch)
    return True
    

def quick_download(out_download_folder):
    """
    Downloads waveform data from the Colombian Seismological Network (SGC) for a predefined time period and location range.

    Parameters
    ----------
    out_download_folder : str
        The path to the folder where the downloaded data will be stored.

    Returns
    -------
    None
        This function does not return any values but saves the downloaded data in the specified folder.

    Notes
    -----
    This function defines a specific time window (from 2019-12-24T19:00 to 2019-12-25T01:00) and a set of stations.
    The data is downloaded in MiniSEED format, and metadata is also saved as a JSON file.

    Side Effects
    ------------
    Downloads waveform data into the given folder, including creating directories if needed.
    """

    sgc_rest = WaveformRestrictions(network="CM",
                    station="URMC,VILL,PRA,ORTC,GARC,FLO2,CHI,YOT",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2019-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2019-12-25T01:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    filter_domain= [-83.101,-64.549,-2.229,14.945],
                    )
    sgc_client = FDSNClient('http://sismo.sgc.gov.co:8080')
    sgc_provider = Provider(sgc_client,sgc_rest)
    md = MseedDownloader(providers=[sgc_provider])
    json_path = os.path.join(out_download_folder,"stations")
    inv,json = md.make_inv_and_json(json_path)
    mseed_storage = os.path.join(out_download_folder,"archive","{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed")
    md.download(mseed_storage,
            picker_args={"batch_size":100,"overlap":0.3,"length":60},
            chunklength_in_sec=7200,n_processor=None)

if __name__ == "__main__":
    clone_seismonitor_data("/home/emmanuel/EDCT/seismonitor_dataset")