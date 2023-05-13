import os
from git import Repo
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader

def clone_seismonitor_data(output_folder,branch):
    git_url = "https://github.com/ecastillot/SeisMonitor.git"
    Repo.clone_from(git_url, output_folder,branch=branch)
    return True
    

def quick_download(out_download_folder):

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