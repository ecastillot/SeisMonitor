from . import utils as ut
from obspy import read_inventory
import copy

class Provider:
    """
    Handles the retrieval and storage of station metadata (inventory) 
    from either an XML file or a client service.
    """
    def __init__(self, client, waveform_restrictions, processing=None, xml=None) -> None:
        """
        Initializes the Provider class.

        Parameters:
        -----------
        client : object
            An ObsPy-compatible client to fetch station metadata.
        waveform_restrictions : WaveformRestrictions
            Object containing selection criteria for waveform data.
        processing : Processing, optional
            Processing steps for waveform data (default is None).
        xml : str, optional
            Path to an XML file containing station metadata (default is None).
        """
        self.client = client
        self.waveform_restrictions = waveform_restrictions
        self.processing = processing
        self.xml = xml

        # Load station inventory from an XML file if provided, otherwise fetch from client
        if xml is not None:
            self.inventory = read_inventory(xml)
        else:
            self.inventory = client.get_stations(
                network=waveform_restrictions.network,
                station=waveform_restrictions.station,
                location="*",
                channel="*",
                level='response'
            )
    
    def copy(self):
        """
        Creates a deep copy of the Provider instance.

        Returns:
        --------
        Provider
            A new instance of Provider with copied attributes.
        """
        return copy.copy(self) # be careful, maybe deepcopy is needed


class WaveformRestrictions:
    """
    Defines criteria for selecting waveform data, including network, station,
    location, and time range constraints.
    """
    def __init__(self, network, station, location, channel, starttime, endtime,
                 location_preferences=None, channel_preferences=None,
                 filter_networks=None, filter_stations=None,
                 filter_domain=None):
        """
        Initializes the waveform selection criteria.
        
        Parameters:
        -----------
        network : str
            Comma-separated list of network codes (wildcards allowed).
        station : str
            Comma-separated list of station codes (wildcards allowed).
        location : str
            Comma-separated list of location identifiers (wildcards allowed).
        channel : str
            Comma-separated list of channel codes (e.g., "BHZ,HHZ").
        starttime : obspy.UTCDateTime
            Start time for waveform selection.
        endtime : obspy.UTCDateTime
            End time for waveform selection.
        location_preferences : list, optional
            Ordered list of preferred locations (default is empty list).
        channel_preferences : list, optional
            Ordered list of preferred channels (default is empty list).
        filter_networks : list, optional
            List of networks to filter out (default is empty list).
        filter_stations : list, optional
            List of stations to filter out (default is empty list).
        filter_domain : list, optional
            Geographic bounding box [lon_west, lon_east, lat_south, lat_north]
            (default is global coverage).
        """
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.starttime = starttime
        self.endtime = endtime
        self.location_preferences = location_preferences or []
        self.channel_preferences = channel_preferences or []
        self.filter_networks = filter_networks or []
        self.filter_stations = filter_stations or []
        self.filter_domain = filter_domain or [-180, 180, -90, 90]


class Processing:
    """
    Defines processing steps to be applied to waveform data.
    """
    def __init__(self, order=None, decimate=None, detrend=None, filter=None, merge=None,
                 normalize=None, resample=None, taper=None,
                 select_networks=None, select_stations=None,
                 filter_networks=None, filter_stations=None):
        """
        Initializes processing steps for waveform data.

        Parameters:
        -----------
        order : list of str, optional
            Order of preprocessing steps (default includes 'normalize', 'merge', etc.).
        decimate : dict, optional
            Parameters for the decimate method.
        detrend : dict, optional
            Parameters for the detrend method.
        filter : dict, optional
            Parameters for the filter method.
        merge : dict, optional
            Parameters for the merge method.
        normalize : dict, optional
            Parameters for the normalize method.
        resample : dict, optional
            Parameters for the resample method.
        taper : dict, optional
            Parameters for the taper method.
        select_networks : list, optional
            List of networks to select (default is empty list).
        select_stations : list, optional
            List of stations to select (default is empty list).
        filter_networks : list, optional
            List of networks to filter out (default is empty list).
        filter_stations : list, optional
            List of stations to filter out (default is empty list).
        """
        self.order = order or ['normalize', 'merge', 'detrend', 'taper', "filter"]
        self.decimate = decimate or {"factor": 2}
        self.detrend = detrend or {"type": "demean"}
        self.filter = filter or {"type": 'bandpass', "freqmin": 1, "freqmax": 45, "corners": 2, "zerophase": True}
        self.merge = merge or {"method": 0, "fill_value": 'latest'}
        self.normalize = normalize or {"global_max": False}
        self.resample = resample or {"sampling_rate": 200}
        self.taper = taper or {"max_percentage": 0.001, "type": "cosine", "max_length": 2}
        self.select_networks = select_networks or []
        self.select_stations = select_stations or []
        self.filter_networks = filter_networks or []
        self.filter_stations = filter_stations or []

    def run(self, st):
        """
        Applies the defined processing steps to a given waveform stream.

        Parameters:
        -----------
        st : obspy.Stream
            The waveform stream to process.
        
        Returns:
        --------
        obspy.Stream
            The processed waveform stream.
        """
        return ut.preproc_stream(
            st, self.order, self.decimate, self.detrend, self.filter,
            self.merge, self.normalize, self.resample, self.taper,
            self.select_networks, self.select_stations,
            self.filter_networks, self.filter_stations
        )