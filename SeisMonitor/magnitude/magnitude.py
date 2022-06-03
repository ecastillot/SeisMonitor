
from obspy.core.inventory.inventory import Inventory
from obspy.io.xseed.parser import Parser
from obspy.core.event import Magnitude as Mag
from obspy.core.event import Catalog
from obspy.core.event import read_events, Comment
import mtspec
import numpy as np
import scipy

class MwPhysicalMagParams():
    def __init__(self,vp=4800.0,
                vsp_factor=1.84,
                density=2700.0,
                waterlevel=10,
                p_radiation_pattern=0.52,
                s_radiation_pattern=0.63):
        self.vp = vp
        self.vsp_factor = vsp_factor
        self.vs = vp/vsp_factor
        self.density = density
        self.waterlevel = waterlevel
        self.p_radiation_pattern = p_radiation_pattern
        self.s_radiation_pattern = s_radiation_pattern

class MwProcessingMagParams():
    def __init__(self,
                time_before_pick = 0.2,
                time_after_pick = 0.8,
                padding = 20):
        self.time_before_pick = time_before_pick
        self.time_after_pick = time_after_pick
        self.padding = padding

def calculate_source_spectrum(frequencies, omega_0, corner_frequency, Q,
    traveltime):
    """
    After Abercrombie (1995) and Boatwright (1980).

    Abercrombie, R. E. (1995). Earthquake locations using single-station deep
    borehole recordings: Implications for microseismicity on the San Andreas
    fault in southern California. Journal of Geophysical Research, 100,
    24003–24013.

    Boatwright, J. (1980). A spectral theory for circular seismic sources,
    simple estimates of source dimension, dynamic stress drop, and radiated
    energy. Bulletin of the Seismological Society of America, 70(1).

    The used formula is:
        Omega(f) = (Omege(0) * e^(-pi * f * T / Q)) / (1 + (f/f_c)^4) ^ 0.5

    :param frequencies: Input array to perform the calculation on.
    :param omega_0: Low frequency amplitude in [meter x second].
    :param corner_frequency: Corner frequency in [Hz].
    :param Q: Quality factor.
    :param traveltime: Traveltime in [s].
    """
    num = omega_0 * np.exp(-np.pi * frequencies * traveltime / Q)
    denom = (1 + (frequencies / corner_frequency) ** 4) ** 0.5
    return num / denom

def fit_spectrum(spectrum, frequencies, 
    traveltime, initial_omega_0,
    initial_f_c,quality_factor=1000):
    """
    Fit a theoretical source spectrum to a measured source spectrum.

    Uses a Levenburg-Marquardt algorithm.

    :param spectrum: The measured source spectrum.
    :param frequencies: The corresponding frequencies.
    :para traveltime: Event traveltime in [s].
    :param initial_omega_0: Initial guess for Omega_0.
    :param initial_f_c: Initial guess for the corner frequency.


    :returns: Best fits and standard deviations.
        (Omega_0, f_c, Omega_0_std, f_c_std)
        Returns None, if the fit failed.
    """
    def f(frequencies, omega_0, f_c):
        return calculate_source_spectrum(frequencies, omega_0, f_c,
                quality_factor, traveltime)
    popt, pcov = scipy.optimize.curve_fit(f, frequencies, spectrum, \
        p0=list([initial_omega_0, initial_f_c]), maxfev=100000)
    if popt is None:
        return None
    return popt[0], popt[1], pcov[0, 0], pcov[1, 1]

def get_M0_magnitude_by_pick(st,picktime,traveltime,
                                phasehint,
                                physparams,
                                procparams):


    if st is None or len(st) != 3:
        return (None,None)

    if phasehint.upper() == "P":
        velocity = physparams.vp
        radiation_pattern = physparams.p_radiation_pattern
    elif phasehint.upper() == "S":
        velocity = physparams.vs
        radiation_pattern = physparams.s_radiation_pattern
    else:
        return (None,None)

    distance = traveltime * velocity
    if distance <= 0.0:
        return None

    omegas = []
    corner_freqs = []
    for trace in st:
        pick_index = int(round((picktime - trace.stats.starttime) / \
                                trace.stats.delta))
        # Choose date window 0.5 seconds before and 1 second after pick.
        data_window = trace.data[pick_index - \
                    int(procparams.time_before_pick * trace.stats.sampling_rate): \
                    pick_index + int(procparams.time_after_pick * trace.stats.sampling_rate)]

        spec, freq = mtspec.mtspec(data_window, trace.stats.delta, 2)
        try:
            fit = fit_spectrum(spec, freq, traveltime,
                    spec.max(), 10.0)
        except:
            continue

        if fit is None:
            continue
        Omega_0, f_c, err, _ = fit
        Omega_0 = np.sqrt(Omega_0)
        omegas.append(Omega_0)
        corner_freqs.append(f_c)

    if omegas:

        M_0 = 4.0 * np.pi * physparams.density * velocity ** 3 * distance * \
            np.sqrt(omegas[0] ** 2 + omegas[1] ** 2 + omegas[2] ** 2) / \
            radiation_pattern

        return (M_0,corner_freqs)
    else:
        return (None,None)

def write_magnitude_values(value,uncertainty,station_count,mag_type,
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           origin_id=None,
                           comments=None):
    mag = Mag()
    mag.mag = value
    mag.mag_errors.uncertainty = uncertainty
    mag.magnitude_type = mag_type
    mag.origin_id = origin_id
    mag.station_count = station_count
    mag.evaluation_mode = evaluation_mode 
    mag.evaluation_status = evaluation_status
    mag.comments.append(Comment( \
        "Magnitude=%e Nm; uncertainty=%e" % (value,
        uncertainty)))
    if comments != None:
        mag.comments.append(Comment(comments))
    return mag

            

class Magnitude():
    def __init__(self,client,catalog,response,
                jsonfile=None) -> None:
        self.client = client
        # super().__init__(sds_archive,sds_type,format)
        self.response = response
        self.jsonfile = jsonfile

        if isinstance(catalog,Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)

    def get_Mw(self,physparams,procparams,
                outfile=None,out_format="QUAKEML"):

        Mws = []
        Mws_std = []

        for event in self.catalog:
            print("event:",event.resource_id)

            if not event.origins:
                print ("No origin for event %s" % event.resource_id)
                continue

            origin_time = event.origins[0].time
            latitude = event.origins[0].latitude
            longitude = event.origins[0].longitude
            depth = event.origins[0].depth

            moments = []
            corner_frequencies = []
            for pick in event.picks:

                
                st = self._get_corresponding_stream(pick.waveform_id, pick.time,
                                                procparams.padding,
                                                physparams.waterlevel)
                if st == None:
                    continue

                traveltime = pick.time - origin_time
                M_0,corner_freqs = get_M0_magnitude_by_pick(st,pick.time,traveltime,
                                                    pick.phase_hint,
                                                    physparams,
                                                    procparams)

                if M_0 != None:
                    staname = ".".join((pick.waveform_id.network_code, pick.waveform_id.station_code,
                              pick.waveform_id.location_code ,pick.waveform_id.channel_code))
                    print(f"\t-> M0 | {staname} | {M_0}")

                    moments.append(M_0)
                    corner_frequencies.extend(corner_freqs)
            
            if not len(moments):
                print (f"No moments could be calculated for event: {event.resource_id}") 
                continue

            # Calculate the seismic moment via basic statistics.
            moments = np.array(moments)
            moment = moments.mean()
            moment_std = moments.std()

            corner_frequencies = np.array(corner_frequencies)
            corner_frequency = corner_frequencies.mean()
            corner_frequency_std = corner_frequencies.std()

            Mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
            Mw_std = 2.0 / 3.0 * moment_std / (moment * np.log(10))
            Mws_std.append(Mw_std)
            Mws.append(Mw)

            print(f"Mw | {Mw} | {Mw_std} | {event.resource_id}")

            mag = write_magnitude_values(Mw,Mw_std,len(moments),"Mw",
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           origin_id=event.origins[0].resource_id,
                           comments=None)

            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

        if outfile != None:
            print ("Writing output file...")
            self.catalog.write(outfile, out_format)
        return self.catalog
                

    def _get_paz_from_response(self,seed_id,
                                datetime=None):
        if isinstance(self.response,Parser):
            try:
                paz = self.response.get_paz(seed_id,datetime)
            except:
                return None

        elif isinstance(self.response,Inventory):
            try:
                response = self.response.get_response(seed_id,
                                        datetime)
            
                paz_stage = response.get_paz()
            except:
                return None

            sensitivity = response.instrument_sensitivity.value
            gain = paz_stage.normalization_factor
            poles = paz_stage.poles
            zeros = paz_stage.zeros

            paz = {'poles': poles,
                    'zeros': zeros,
                    'gain': gain,
                    'sensitivity': sensitivity}
        return paz


    def _get_corresponding_stream(self,waveform_id, 
                                pick_time,
                                padding=1.0,
                                waterlevel=10):
        """
        Helper function to find a requested waveform in the previously created
        waveform_index file.
        Also performs the instrument correction.
        Returns None if the file could not be found.
        """
        
        start = pick_time - padding
        end = pick_time + padding
        try:
            st= self.client.get_waveforms(waveform_id.network_code,
                                    waveform_id.station_code,
                                    "*",
                                    waveform_id.channel_code[0:2]+"*",start,end
                                    )
            # print(f"\t->OK: {waveform_id.network_code}-{waveform_id.station_code}"+\
            #         f"-*-{waveform_id.channel_code[0:2]}* |"+\
            #             f"{start} - {end}  ")
        except:
            print(f"\t->Not found: {waveform_id.network_code}-{waveform_id.station_code}"+\
                    f"-*-{waveform_id.channel_code[0:2]}* |"+\
                        f"{start} - {end}  ")
            return None

        for trace in st:
            paz = self._get_paz_from_response(trace.id, start)

            if paz == None:
                print(f"\t->No response found: {trace.id}-{start}")
                return None

            # PAZ in SEED correct to m/s. Add a zero to correct to m.
            paz["zeros"].append(0 + 0j)
            trace.detrend()
            trace.simulate(paz_remove=paz, water_level=waterlevel)
        return st

if __name__ =="__main__":
    from obspy.clients.fdsn import Client
    from obspy.core.inventory.inventory import read_inventory
    client = 'http://sismo.sgc.gov.co:8080'
    client = Client(client)

    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022knomqj.xml"
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022cvykgs.xml"
    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022ktsrvi.xml"
    catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022krnhiu.xml"
    # resp = "/home/emmanuel/Ecopetrol/SeisMonitor/data/metadata/CM.dataless"
    resp = "/home/emmanuel/EDCT/SeisMonitor/data/events/public_CM.xml"
    resp = read_inventory(resp)
    
    out="./test_magnitude,xml"
    mag = Magnitude(client,catalog,resp,
                jsonfile=None) 

    physparams = MwPhysicalMagParams()
    procparams = MwProcessingMagParams()
    mag.get_Mw(physparams,procparams,out)