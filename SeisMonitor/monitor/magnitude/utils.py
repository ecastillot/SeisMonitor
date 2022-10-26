from obspy.core.inventory.inventory import Inventory
from obspy.io.xseed.parser import Parser
from obspy.core.event import Magnitude as Mag
from obspy.core.event.magnitude import StationMagnitude
from obspy.core.event import Catalog
from obspy.core.event import read_events, Comment
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import estimate_wood_anderson_amplitude_using_response
from obspy import UTCDateTime
from obspy.core.event.base import CreationInfo,WaveformStreamID,TimeWindow
from obspy.core.event.resourceid import ResourceIdentifier

from obspy.core.event.magnitude import Amplitude
import mtspec
import numpy as np
import scipy
import math 
import types
paz_wa = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}

# SED_Ml_params = {"a":0.018,"b":2.17}
Ml_params = {"RSNC":{"a":1.019,"b":0.0016,"r_ref":140},
            "RSNC_by_zone":{"Ml_1":{"a":1.2488,"b":0.0024,"c":-2.05},
                            "Ml_2":{"a":1.0563,"b":0.002,"c":-1.760},
                            "Ml_3":{"a":1.0705,"b":0.0013,"c":-1.531},
                            "Ml_4":{"a":1.2399,"b":0.0015,"c":-2.178},
                            "Ml_5":{"a":0.7096,"b":0.0009,"c":-0.690},
                            }
            }

def get_paz_from_response(seed_id,response,
                            datetime=None):
    if isinstance(response,Parser):
        try:
            paz = response.get_paz(seed_id,datetime)
        except:
            return None

    elif isinstance(response,Inventory):
        try:
            response = response.get_response(seed_id,
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

def get_Ml(ampl,epi_dist,mag_type,zone=None):

    if isinstance(mag_type,str):
        if mag_type not in list(Ml_params.keys()):
            raise Exception(f"{mag_type} not in Ml_params")

        if mag_type.upper() == "RSNC":

            if zone != None:
                
                ml_params = Ml_params["RSNC_by_zone"][f"Ml_{zone}"]

                return math.log10(ampl * 1e9) +  ml_params["a"] * math.log10(epi_dist) +\
                       ml_params["b"] * epi_dist  + ml_params["c"]




            else:
                ml_params = Ml_params["RSNC"]

                k = ml_params["a"]*math.log10(ml_params["r_ref"]/100) +\
                    ml_params["b"]* (ml_params["r_ref"]-100) +3
                    
                return math.log10(ampl * 1e3) + ml_params["a"] * math.log10(epi_dist/ml_params["r_ref"]) +\
                                ml_params["b"] * (epi_dist-ml_params["r_ref"]) + k


    elif isinstance(mag_type,types.LambdaType):
        return mag_type(ampl,epi_dist)

def get_Ml_magparams_by_station(st,response,
                                ev_params,
                                trimmedtime=50,
                                waterlevel=10):
    picktime = ev_params["picktime"]
    latitude = ev_params["latitude"]
    longitude = ev_params["longitude"]

    # if st is None or len(st) != 3:
    if st is None or len(st)==0:
        return (None,None,None)

    for tr in st:
        paz = get_paz_from_response(tr.id,response,
                                tr.stats.starttime)
        if paz == None:
            print(f"\t->No response found: {tr.id}-{tr.stats.starttime}")
            return (None,None,None)

        coords = response.get_coordinates(tr.id,tr.stats.starttime)
        # paz_wa = estimate_wood_anderson_amplitude_using_response(response, amplitude,
        #                                             timespan)
        tr.simulate(paz_remove=paz, 
                    paz_simulate=paz_wa, 
                    water_level=waterlevel)


    st.trim(picktime,picktime+trimmedtime)

    # common_channels_info = list(st._get_common_channels_info().keys())
    # channels = list(map(lambda x: x[-1] ,common_channels_info))
    components = [ tr.stats.channel[-1] for tr in st]
    components = list(set(components))

    if len(components) == 1:
        tr= st.select(component=components[0])[0]
        ampl = max(abs(tr.data))
        tr_id = tr.get_id()
    else:
        if ("N" in components) or ("E" in components):
            if "N" in components:
                tr_n = st.select(component="N")[0]
                tr_n_id = tr_n.get_id()
                ampl_n = max(abs(tr_n.data))
            else:
                ampl_n = None
                tr_n_id = None

            if "E" in components:
                tr_e = st.select(component="E")[0]
                tr_e_id = tr_e.get_id()
                ampl_e = max(abs(tr_e.data))
            else:
                ampl_e = None
                tr_e_id = None

            if ampl_n == None:
                ampl = ampl_e
                tr_id = tr_e_id
            elif ampl_e == None:
                ampl = ampl_n
                tr_id = tr_n_id
            else:
                ampls = [ampl_n,ampl_e]
                trs_id = [tr_n_id,tr_e_id]
                ampl = max(ampls)
                tr_id = trs_id[ampls.index(ampl)]
        elif "Z" in components:
            tr= st.select(component="Z")[0]
            ampl = max(abs(tr.data))
            tr_id = tr.get_id()
        else:
            return (None,None,None)


    sta_lat = coords["latitude"]
    sta_lon = coords["longitude"]

    epi_dist, az, baz = gps2dist_azimuth(latitude, longitude,
                                         sta_lat, sta_lon)
    epi_dist = epi_dist / 1e3

    return ampl,epi_dist,tr_id

def Mw_st_processing(st,response,waterlevel,datetime):
    for trace in st:
        paz = get_paz_from_response(trace.id, response,datetime)

        if paz == None:
            print(f"\t->No response found: {trace.id}-{datetime}")
            return None

        # PAZ in SEED correct to m/s. Add a zero to correct to m.
        paz["zeros"].append(0 + 0j)
        trace.detrend()
        trace.simulate(paz_remove=paz, water_level=waterlevel)
    # st.plot()
    return st

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
    initial_f_c,quality_factor=500):
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

def write_magsta_values(value,
                        uncertainty,
                        mag_type,
                        origin_id = ResourceIdentifier(),
                        amplitude_id = ResourceIdentifier(),
                        method_id = ResourceIdentifier(),
                        waveform_id = WaveformStreamID(),
                        agency=None ):
    stamag = StationMagnitude()
    stamag.resource_id = ResourceIdentifier()
    stamag.mag = value
    stamag.mag_errors.uncertainty = uncertainty
    stamag.station_magnitude_type = mag_type
    stamag.origin_id = origin_id
    stamag.method_id = method_id
    stamag.amplitude_id = amplitude_id 
    stamag.waveform_id = waveform_id
    stamag.creation_info = CreationInfo(agency_id=agency,
										agency_uri=ResourceIdentifier(id=agency),
										author="SeisMonitor",
										author_uri=ResourceIdentifier(id="SeisMonitor"),
										creation_time=UTCDateTime.now())
    return stamag

def write_amplitude_values(value,amp_type="A",category="duration",
                           unit="m/s",
                           method_id = "https://docs.obspy.org/tutorial/advanced_exercise/advanced_exercise_solution_3b.html",
                           time_window=TimeWindow(),
                           pick_id = ResourceIdentifier(),
                           waveform_id = WaveformStreamID(),
                           magnitude_hint = "M",
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           agency=None ):
    amp = Amplitude()
    amp.resource_id = ResourceIdentifier()
    amp.generic_amplitude = value
    amp.type = amp_type
    amp.category = category
    amp.unit = unit
    amp.method_id = method_id
    amp.time_window = time_window
    amp.pick_id = pick_id 
    amp.waveform_id = waveform_id
    amp.magnitude_hint  = magnitude_hint
    amp.evaluation_mode = evaluation_mode
    amp.evaluation_status = evaluation_status
    amp.creation_info = CreationInfo(agency_id=agency,
										agency_uri=ResourceIdentifier(id=agency),
										author="SeisMonitor",
										author_uri=ResourceIdentifier(id="SeisMonitor"),
										creation_time=UTCDateTime.now())
    return amp

def write_magnitude_values(value,uncertainty,station_count,mag_type,
                           evaluation_mode = "automatic",
                           evaluation_status = "preliminary",
                           method_id=ResourceIdentifier(),
                           agency=None,
                           origin_id=ResourceIdentifier(),
                           comments=None):
    mag = Mag()
    mag.mag = value
    mag.mag_errors.uncertainty = uncertainty
    mag.magnitude_type = mag_type
    mag.origin_id = origin_id
    mag.station_count = station_count
    mag.evaluation_mode = evaluation_mode 
    mag.evaluation_status = evaluation_status
    mag.method_id = method_id
    mag.comments.append(Comment( \
        "Magnitude=%e Nm; uncertainty=%e" % (value,
        uncertainty)))
    if comments != None:
        mag.comments.append(Comment(comments))
    mag.creation_info = CreationInfo(agency_id=agency,
										agency_uri=ResourceIdentifier(id=agency),
										author="SeisMonitor",
										author_uri=ResourceIdentifier(id="SeisMonitor"),
										creation_time=UTCDateTime.now())
    
    return mag