from pathlib import Path
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient
from SeisMonitor.core.objects import WaveformRestrictions,Provider
from SeisMonitor.monitor.downloader.seismonitor import MseedDownloader


monitor_path = Path(__file__).parent / "sm2"


downloads_path = monitor_path / "downloads"
stations_path = monitor_path / "stations"


client = FDSNClient("http://rtserve.beg.utexas.edu")
chunklength_in_sec = 3600
rest = WaveformRestrictions(network="*",
                    station="PB13,PB17,PB20,PB23,PB34,PB58,PB40,PB24,PB25,PB26,WB03",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2022-11-16T21:00:00.000000Z"),
                    endtime=UTCDateTime("2022-11-16T23:00:00.000000Z"),
                    location_preferences=["","00","20","10","40"],
                    channel_preferences=["HH","BH","EH","HN","HL"],
                    filter_networks=[], 
                    filter_stations=[],
                    )


####### Default configuration
provider = Provider(client,rest)
md = MseedDownloader(providers=[provider])
inv,json = md.make_inv_and_json(stations_path)

mseed_storage = downloads_path / "{station}" / "{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed"
md.download(str(mseed_storage),
            picker_args={"batch_size":10,
                        "overlap":0.3,
                        "length":60},
            chunklength_in_sec=chunklength_in_sec,
            n_processor=1)