import os 
import glob
from datetime import timedelta
from obspy.clients.filesystem.sds import Client
from obspy.core.util.misc import BAND_CODE

class LocalClient(Client):

    def __init__(self, root, fmt, **kwargs):
        """
        Parameters:
        -----------
        root: str 
            Path where is located the Local structure
        
        fmt: str 
            The parameter should name the corresponding keys of the stats object, e.g. "{year}-{month:02d}/{year}-{month:02d}-{day:02d}/{network}.{station}.{location}.{channel}.{year}.{julday:03d}"
        
        kwargs: 
            SDS client additional args
        """
        self.root = root
        self.fmt = fmt
        super().__init__(root,**kwargs)

    def _get_filenames(self, network, station, location, channel, starttime,
                       endtime, sds_type=None):
        """
        Get a list of filenames for a given waveform and time span.

        Parameters
        ----------
        network: str
            Network code of requested data (e.g., "IU").
        
        station: str
            Station code of requested data (e.g., "ANMO").
        
        location: str
            Location code of requested data (e.g., "").
        
        channel: str
            Channel code of requested data (e.g., "HHZ").
        
        time: obspy.core.utcdatetime.UTCDateTime
            Time of interest.
        
        sds_type: str
            SDS type (description not provided).
        """
        sds_type = sds_type or self.sds_type
        # SDS has data sometimes in adjacent days, so also try to read the
        # requested data from those files. Usually this is only a few seconds
        # of data after midnight, but for now we play safe here to catch all
        # requested data (and with MiniSEED - the usual SDS file format - we
        # can use starttime/endtime kwargs anyway to read only desired parts).
        year_doy = set()
        # determine how far before starttime/after endtime we should check
        # other dayfiles for the data
        t_buffer = self.fileborder_samples / BAND_CODE.get(channel[:1], 20.0)
        t_buffer = max(t_buffer, self.fileborder_seconds)
        t = starttime - t_buffer
        t_max = endtime + t_buffer
        # make a list of year/doy combinations that covers the whole requested
        # time window (plus day before and day after)
        while t < t_max:
            year_doy.add((t.year,t.month,t.day, t.julday))
            t += timedelta(days=1)
        year_doy.add((t_max.year,t_max.month,t_max.day, t_max.julday))

        full_paths = set()
        for year,month,day,doy in year_doy:
            filename = self.fmt.format(
                            network=network, station=station, location=location,
                            channel=channel, year=year, month=month, 
                            day=day, julday=doy,sds_type=sds_type)
            full_path = os.path.join(self.sds_root, filename)
            full_paths = full_paths.union(glob.glob(full_path))
        
        return full_paths

    def _get_filename(self, network, station, location, channel, time, sds_type=None):
        """
        Get the filename for a given waveform.

        Parameters
        ----------
        network: str
            Network code of requested data (e.g., "IU").
        
        station: str
            Station code of requested data (e.g., "ANMO").
        
        location: str
            Location code of requested data (e.g., "").
        
        channel: str
            Channel code of requested data (e.g., "HHZ").
        
        time: obspy.core.utcdatetime.UTCDateTime
            Time of interest.

        Returns
        -------
        str
            The filename corresponding to the given parameters.
        """
        sds_type = sds_type or self.sds_type
        filename = self.fmt.format(
                    network=network, station=station, location=location,
                    channel=channel, year=time.year, month=time.month, 
                    day=time.day, doy=time.julday,sds_type=sds_type)
        return os.path.join(self.sds_root, filename)
