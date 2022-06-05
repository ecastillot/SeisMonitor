# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-10-12 10:58:25
#  * @modify date 2021-10-12 10:58:25
#  * @desc [description]
#  */
######## 
# Moment magnitude is adapted from 
# https://github.com/krischer/moment_magnitude_calculator/tree/master/scripts
# :copyright:
# Lion Krischer (krischer@geophysik.uni-muenchen.de), 2012

# import colorama
import matplotlib.pylab as plt
import mtspec
import numpy as np
# from obspy.core.event import Magnitude as Mag, catalog
from obspy.core.event import Magnitude as Mag
from obspy.core.event import Catalog
from obspy.core.event import read_events, Comment
from obspy.io.xseed.parser import Parser
from obspy.core.event import WaveformStreamID
# import progressbar
import scipy
import scipy.optimize
# import warnings
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client
from obspy.io.xseed.parser import Parser
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import estimate_magnitude
from math import log10
import concurrent.futures as cf
import numpy as np
import time
import json

paz_wa = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}

def fit_spectrum(spectrum, frequencies, traveltime, initial_omega_0,
    initial_f_c,quality_factor=10000):
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

def fit_moment_magnitude_relation_curve(Mls, Mws, Mw_stds):
    """
    Fits a quadratic curve to
        Mw = a + b * Ml + c * Ml ** 2
    Returns the best fitting [a, b, c]
    """
    def y(x, a, b, c):
        return a + b * x + c * x ** 2
    # Use a straight line as starting point.
    Mls = np.ma.masked_invalid(Mls)
    Mws = np.ma.masked_invalid(Mws)
    inds = ~(Mls.mask | Mws.mask | np.isnan(Mw_stds) | (Mw_stds <= 0))
    popt, pcov = scipy.optimize.curve_fit(y, Mls[inds], Mws[inds], \
        p0=[0.0, 1.0, 0.0], sigma=Mw_stds[inds], maxfev=100000)
    return popt[0], popt[1], popt[2]

def plot_ml_vs_mw(catalog):
    moment_magnitudes = []
    moment_magnitudes_std = []
    local_magnitudes = []
    local_magnitudes_std = []
    for event in catalog:
        Mw = None
        Mw_std = None
        Ml = None
        Ml_std = None
        for mag in event.magnitudes:
            if Mw is not None and Ml is not None:
                break
            mag_type = mag.magnitude_type.lower()
            if mag_type == "mw":
                if Mw is not None:
                    continue
                Mw = mag.mag
                Mw_std = mag.mag_errors.uncertainty
            elif mag_type == "ml":
                if Ml is not None:
                    continue
                Ml = mag.mag
                Ml_std = mag.mag_errors.uncertainty
        moment_magnitudes.append(Mw)
        moment_magnitudes_std.append(Mw_std)
        local_magnitudes.append(Ml)
        local_magnitudes_std.append(Ml_std)
    moment_magnitudes = np.array(moment_magnitudes, dtype="float64")
    moment_magnitudes_std = np.array(moment_magnitudes_std, dtype="float64")
    local_magnitudes = np.array(local_magnitudes, dtype="float64")
    local_magnitudes_std = np.array(local_magnitudes_std, dtype="float64")

    # Fit a curve through the data.
    a, b, c = fit_moment_magnitude_relation_curve(local_magnitudes,
            moment_magnitudes, moment_magnitudes_std)
    x_values = np.linspace(-2.0, 4.0, 10000)
    fit_curve = a + b * x_values + c * x_values ** 2

    plt.figure(figsize=(10, 8))
    # Show the data values as dots.
    plt.scatter(local_magnitudes, moment_magnitudes, color="blue",
        edgecolor="black")

    # Plot the Ml=Mw line.
    plt.plot(x_values, x_values, label="$Mw=Ml$", color="k", alpha=0.8)
    plt.plot(x_values, 0.67 + 0.56 * x_values + 0.046 * x_values ** 2,
        label="$Mw=0.67 + 0.56Ml + 0.046Ml^2 (gruenthal etal 2003)$", color="green", ls="--")
    plt.plot(x_values, 0.53 + 0.646 * x_values + 0.0376 * x_values ** 2,
        label="$Mw=0.53 + 0.646Ml + 0.0376Ml^2 (gruenthal etal 2009)$", color="green")
    plt.plot(x_values, 0.594 * x_values + 0.985,
        label="$Mw=0.985 + 0.594Ml (goertz-allmann etal 2011)$", color="orange")
    plt.plot(x_values, (x_values + 1.21) / 1.58,
        label="$Mw=(Ml + 1.21) / 1.58 (bethmann etal 2011)$", color="red")
    plt.plot(x_values, fit_curve, color="blue",
        label="$Data$ $fit$ $with$ $Mw=%.2f + %.2fMl + %.3fMl^2$" % (a, b, c))
    # Set limits and labels.
    plt.xlim(-2, 4)
    plt.ylim(-2, 4)
    plt.xlabel("Ml", fontsize="x-large")
    plt.ylabel("Mw", fontsize="x-large")
    # Show grid and legend.
    plt.grid()
    plt.legend(loc="lower right")
    plt.savefig("moment_mag_automatic.pdf")

def plot_source_radius(cat):
    mw = []
    mw_std = []
    source_radius = []
    source_radius_std = []

    plt.figure(figsize=(10, 4.5))

    # Read the source radius.
    for event in cat:
        mag = event.magnitudes[1]
        if len(mag.comments) != 2:
            continue
        mw.append(mag.mag)
        mw_std.append(mag.mag_errors.uncertainty)
        sr, std = mag.comments[1].text.split(";")
        _, sr = sr.split("=")
        _, std = std.split("=")
        sr = float(sr[:-1])
        std = float(std)
        source_radius.append(sr)
        source_radius_std.append(std)
    plt.errorbar(mw, source_radius, yerr=source_radius_std,
        fmt="o", linestyle="None")
    plt.xlabel("Mw", fontsize="x-large")
    plt.ylabel("Source Radius [m]", fontsize="x-large")
    plt.grid()
    plt.savefig("/Users/lion/Desktop/SourceRadius.pdf")

class Magnitude():
    def __init__(self,client,cat,resp_file,
                jsonfile=None) -> None:
        self.client = client
        # super().__init__(sds_archive,sds_type,format)
        self.resp_file = resp_file
        self.cat = cat
        self.jsonfile = jsonfile
    
    @property
    def info_from_json(self):
        if self.jsonfile == None:
            return None
        else:
            f = open (self.jsonfile, "r")
            return json.load(f)

    @property
    def catalog(self):
        if isinstance(self.cat,Catalog):
            return self.cat
        else:
            return read_events(self.cat)

    @property
    def resp(self):
        return Parser(self.resp_file)

    def _get_channel_resp(self):
        parser = self.resp
        channels = [c['channel_id'] for c in parser.get_inventory()['channels']]
        # print(channels)
        # exit()
        parsers = dict.fromkeys(channels, parser)
        print(parsers)
        exit()
        return parsers

    def _get_paz_wa(self,seed_id):
        try:
            parse = self.resp
        except Exception as e:
            print(e)
        
        return parse.get_paz(seed_id)
    
    def _get_coordinates(self,seed_id):
        try:
            parse = self.resp
        except Exception as e:
            print(e)
        
        return parse.get_coordinates(seed_id)

    def _estimate_local_CC_magnitude(self,event_coords:tuple,times2trim:dict,
                            origin_time:UTCDateTime):
        event_lat,event_lon,event_depth = event_coords

        mags = []
        def compute_ml(t2t):
        # for info,times in times2trim.items():
            info,times = t2t
            net,sta,loc,cha = info.split(".")
            starttime,endtime = times
            try: 
                st = self.client.get_waveforms(net,sta,loc,cha, starttime,
                        endtime, metadata=True)
                assert(len(st) == 3)
            except Exception:
                print(info, "---")
                return

            amps = []
            pazs = []
            for tr in st:
                single_cha = ".".join((tr.stats.network,
                                        tr.stats.station,
                                        tr.stats.location,
                                        tr.stats.channel ))
                paz = self._get_paz_wa(single_cha)
                coords = self._get_coordinates(single_cha)

                # tr.simulate(paz_remove=paz, paz_simulate=paz_wa, 
                #             water_level=10)
                if tr.stats.channel[-1] != "Z":
                    amp =  max(abs(tr.data))
                    amps.append(amp)
                    pazs.append(paz)

            epi_dist, az, baz = gps2dist_azimuth(event_lat, event_lon,
                                             coords["latitude"], coords["longitude"])
            epi_dist = epi_dist / 1000

            timespan = endtime-starttime
            ml = estimate_magnitude(pazs,amps,[timespan,timespan],epi_dist)

            # print(sta, ml)
            mags.append(ml)
        
        with cf.ThreadPoolExecutor() as executor:
            executor.map(compute_ml,times2trim.items())

        net_mag = np.median(mags)
        print(f"{origin_time}: ({event_lat},{event_lon},{event_depth}) | Ml:{round(net_mag,3)}")
        return round(net_mag,3)

    def _estimate_local_SED_magnitude(self,event_coords:tuple,times2trim:dict,waterlevel=10):
        event_lat,event_lon = event_coords

        mags = []
        for info,times in times2trim.items():
            net,sta,loc,cha = info.split(".")
            starttime,endtime = times
            try: 
                st = self.client.get_waveforms(net,sta,loc,cha, starttime,
                        endtime, metadata=True)
                assert(len(st) == 3)
            except Exception:
                print(info, "---")
                continue

            amps = []
            for tr in st:
                single_cha = ".".join((tr.stats.network,
                                        tr.stats.station,
                                        tr.stats.location,
                                        tr.stats.channel ))
                paz = self._get_paz_wa(single_cha)
                coords = self._get_coordinates(single_cha)

                tr.simulate(paz_remove=paz, paz_simulate=paz_wa, 
                            water_level=waterlevel)
                if tr.stats.channel[-1] != "Z":
                    amp =  max(abs(tr.data))
                    amps.append(amp)

            ampl = max(amps)

            epi_dist, az, baz = gps2dist_azimuth(event_lat, event_lon,
                                             coords["latitude"], coords["longitude"])
            epi_dist = epi_dist / 1000


            a = 0.018
            b = 2.17
            ml = log10(ampl * 1000) + a * epi_dist + b
            print(sta, ml)
            mags.append(ml)
        
        net_mag = np.median(mags)
        print("Network magnitude:", net_mag)
        return net_mag

    def _estimate_moment_magnitudes(self, outfile=None,
                                vp=4266.96,vsp_factor=1.84,
                                density=2650.0,
                                time_before_pick = 0.2,
                                time_after_pick = 0.8,
                                padding = 20,waterlevel=10,
                                out_format="QUAKEML"):
        """
        :param cat: obspy.core.event.Catalog object.
        :param outfile: str

        # # Rock density in km/m^3.
        # # Velocities in m/s.

        # # time_before/after_pick: How many seconds before and after the pick
        #  to choose for calculating the spectra.

        # # Fixed quality factor. Very unstable inversion for it. Has almost no influence
        # # on the final seismic moment estimations but has some influence on the corner
        # # frequency estimation and therefore on the source radius estimation.
        """
        cat = self.catalog
        vs = vp/vsp_factor
        Mws = []
        Mls = []
        Mws_std = []
        print("wait",catalog)

        for event in cat:
            print(event)
            if not event.origins:
                print ("No origin for event %s" % event.resource_id)
                continue
            # if not event.magnitudes:
            #     print ("No magnitude for event %s" % event.resource_id)
            #     continue

            origin_time = event.origins[0].time
            latitude = event.origins[0].latitude
            longitude = event.origins[0].longitude
            depth = event.origins[0].depth
            # local_magnitude = event.magnitudes[0].mag

            #if local_magnitude < 1.0:
                #continue
            moments = []
            source_radii = []
            corner_frequencies = []
            # print(event.picks)
            # exit()
            for pick in event.picks:
                # Only p phase picks.
                if pick.phase_hint.lower() == "p":
                    radiation_pattern = 0.52
                    velocity = vp
                    k = 0.32
                elif pick.phase_hint.lower() == "s":
                    radiation_pattern = 0.63
                    velocity = vs
                    k = 0.21
                else:
                    continue
                distance = (pick.time - origin_time) * velocity
                if distance <= 0.0:
                    continue

                # wav_id = pick.waveform_id
                if self.jsonfile != None:
                    wav_id = pick.waveform_id
                    sta = wav_id.station_code
                    net = self.info_from_json[sta]["network"]
                    cha = self.info_from_json[sta]["channels"][0][:2]+"*"
                    pick.waveform_id = WaveformStreamID(network_code=net,
                                                    station_code=sta,
                                                    location_code="*",
                                                    channel_code=cha)
                # print(pick.waveform_id, pick.time,
                #                                 padding,waterlevel)

                stream = self.get_corresponding_stream(pick.waveform_id, pick.time,
                                                padding,waterlevel)
                # exit()
                # print(pick.waveform_id)
                # print(stream)
                if stream is None or len(stream) != 3:
                    continue
                omegas = []
                corner_freqs = []
                for trace in stream:
                    # print(pick.time)
                    # print(trace.stats.starttime)
                    # print(trace.stats.delta)
                    # Get the index of the pick.
                    pick_index = int(round((pick.time - trace.stats.starttime) / \
                        trace.stats.delta))
                    # print(len(trace.data))
                    # print(pick_index - \
                    #     int(time_before_pick * trace.stats.sampling_rate))
                    # print(pick_index + int(time_after_pick * trace.stats.sampling_rate))
                    # Choose date window 0.5 seconds before and 1 second after pick.
                    print(trace.data)
                    data_window = trace.data[pick_index - \
                        int(time_before_pick * trace.stats.sampling_rate): \
                        pick_index + int(time_after_pick * trace.stats.sampling_rate)]
                    # print(data_window)
                    # Calculate the spectrum.
                    spec, freq = mtspec.mtspec(data_window, trace.stats.delta, 2)
                    try:
                        fit = fit_spectrum(spec, freq, pick.time - origin_time,
                                spec.max(), 10.0)
                    except:
                        continue
                    if fit is None:
                        continue
                    # print(fit)
                    Omega_0, f_c, err, _ = fit
                    Omega_0 = np.sqrt(Omega_0)
                    omegas.append(Omega_0)
                    corner_freqs.append(f_c)
                
                if omegas:
                    M_0 = 4.0 * np.pi * density * velocity ** 3 * distance * \
                        np.sqrt(omegas[0] ** 2 + omegas[1] ** 2 + omegas[2] ** 2) / \
                        radiation_pattern
                    r = 3 * k * vs / sum(corner_freqs)
                    moments.append(M_0)
                    source_radii.append(r)
                    corner_frequencies.extend(corner_freqs)
                else:
                    pass
            if not len(moments):
                print ("No moments could be calculated for event") 
                # print ("No moments could be calculated for event %s") % \
                #     event.resource_id.resource_id
                continue

            # Calculate the seismic moment via basic statistics.
            moments = np.array(moments)
            moment = moments.mean()
            moment_std = moments.std()

            corner_frequencies = np.array(corner_frequencies)
            corner_frequency = corner_frequencies.mean()
            corner_frequency_std = corner_frequencies.std()

            # Calculate the source radius.tation_code
            source_radii = np.array(source_radii)
            source_radius = source_radii.mean()
            source_radius_std = source_radii.std()

            # Calculate the stress drop of the event based on the average moment and
            # source radii.
            stress_drop = (7 * moment) / (16 * source_radius ** 3)
            stress_drop_std = np.sqrt((stress_drop ** 2) * \
                (((moment_std ** 2) / (moment ** 2)) + \
                (9 * source_radius * source_radius_std ** 2)))
            # if source_radius > 0 and source_radius_std < source_radius:
            #     print ("Source radius:", source_radius, " Std:", source_radius_std)
            #     print ("Stress drop:", stress_drop / 1E5, " Std:", stress_drop_std / 1E5)

            Mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
            Mw_std = 2.0 / 3.0 * moment_std / (moment * np.log(10))
            Mws_std.append(Mw_std)
            Mws.append(Mw)
            # Mls.append(local_magnitude)
            # calc_diff = abs(Mw - local_magnitude)
            # Mw = ("%.3f" % Mw).rjust(7)
            # Ml = ("%.3f" % local_magnitude).rjust(7)
            # diff = ("%.3e" % calc_diff).rjust(7)
            # ret_string = colorama.Fore.GREEN + \
            #     "For event %s: Ml=%s | Mw=%s | " % (event.resource_id.resource_id,
            #     Ml, Mw)
            # if calc_diff >= 1.0:
                # ret_string += colorama.Fore.RED
            # ret_string += "Diff=%s" % diff
            # ret_string += colorama.Fore.GREEN
            # ret_string += " | Determined at %i stations" % len(moments)
            # ret_string += colorama.Style.RESET_ALL
            # print (ret_string)
            # print("Mw", Mw)
            print(f"{origin_time}: ({latitude},{longitude},{depth}) | Mw:{round(Mw,3)}")
            mag = Mag()
            mag.mag = Mw
            mag.mag_errors.uncertainty = Mw_std
            mag.magnitude_type = "Mw"
            mag.origin_id = event.origins[0].resource_id
            mag.method_id = "smi:com.github/krischer/moment_magnitude_calculator/automatic/1"
            mag.station_count = len(moments)
            mag.evaluation_mode = "automatic"
            mag.evaluation_status = "preliminary"
            mag.comments.append(Comment( \
                "Seismic Moment=%e Nm; standard deviation=%e" % (moment,
                moment_std)))
            mag.comments.append(Comment("Custom fit to Boatwright spectrum"))
            if source_radius > 0 and source_radius_std < source_radius:
                mag.comments.append(Comment( \
                    "Source radius=%.2fm; standard deviation=%.2f" % (source_radius,
                    source_radius_std)))
            
            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

        if outfile != None:
            print ("Writing output file...")
            cat.write(outfile, out_format)
        return cat

    def get_corresponding_stream(self,waveform_id, pick_time, padding=1.0,
                                waterlevel=10):
        """
        Helper function to find a requested waveform in the previously created
        waveform_index file.
        Also performs the instrument correction.
        Returns None if the file could not be found.
        """
        
        # from obspy.clients.filesystem.sds import Client as SDS_Client
        # sds_archive = "/home/emmanuel/archive/sds"
        # client = SDS_Client(sds_archive,
        #             sds_type='D', format='MSEED')
        start = pick_time - padding
        end = pick_time + padding
        try:
            st= self.client.get_waveforms(waveform_id.network_code,
                                    waveform_id.station_code,
                                    "*",
                                    waveform_id.channel_code[0:2]+"*",start,end
                                    )
            print(f"OK: {waveform_id.network_code}-{waveform_id.station_code}"+\
                    f"-*-{waveform_id.channel_code[0:2]}* |"+\
                     f"{start} - {end}  ")
        except:
            print(f"Not found: {waveform_id.network_code}-{waveform_id.station_code}"+\
                    f"-*-{waveform_id.channel_code[0:2]}* |"+\
                     f"{start} - {end}  ")
            return None

        for trace in st:
            parsers = self._get_channel_resp()
            # print(trace,parsers)
            paz = parsers[trace.id].get_paz(trace.id, start)
            print("paz",paz)
            # PAZ in SEED correct to m/s. Add a zero to correct to m.
            paz["zeros"].append(0 + 0j)
            trace.detrend()
            trace.simulate(paz_remove=paz, water_level=waterlevel)
        return st

    def estimate_local_magnitude(self,outfile=None,padding=0.5,
                        out_format="QUAKEML",n_processor=None):
        cat = self.catalog
        
        tic = time.time()

        
        # def local_main(event):
        for event in cat:
            if not event.origins:
                print ("No origin for event %s" % event.resource_id)
                continue
            # if not event.magnitudes:
            #     print ("No magnitude for event %s" % event.resource_id)
            #     continue

            origin = event.preferred_origin()
            lat = origin.latitude
            lon = origin.longitude
            depth = origin.depth
            if lat ==None or lon==None:
                continue
            # local_magnitude = event.magnitudes[0].mag

            #if local_magnitude < 1.0:
                #continue
            times2trim={}
            for pick in event.picks:
                if pick.phase_hint.lower() == "s":
                    wav_id = pick.waveform_id

                    if self.jsonfile != None:
                        net = self.info_from_json[wav_id.station_code]["network"]
                        cha = self.info_from_json[wav_id.station_code]["channels"][0][:2]+"*"
                    else: 
                        net = wav_id.network_code
                        cha = wav_id.channel_code[:2]+"*"


                    key = ".".join((net,
                                    wav_id.station_code,"*",
                                    cha))
                    start = pick.time - padding
                    end = pick.time + padding
                    times2trim[key]=(start,end)
            if not times2trim:
                continue

            Ml = self._estimate_local_CC_magnitude((lat,lon,depth),times2trim,origin.time)
                # print(key)

            mag = Mag()
            mag.mag = Ml
        #     mag.mag_errors.uncertainty = Mw_std
            mag.magnitude_type = "Ml"
            mag.origin_id = event.origins[0].resource_id
            mag.method_id = "smi:com.github/local_mag"
            mag.station_count = len(times2trim)
            mag.evaluation_mode = "automatic"
            mag.evaluation_status = "preliminary"
            mag.comments.append(Comment( \
                "local Moment=%e Nm" % (Ml)))
            event.magnitudes.append(mag)
            event.preferred_magnitude_id = mag.resource_id

        # with cf.ThreadPoolExecutor(n_processor) as executor:
        #     executor.map(local_main,cat)

        toc = time.time()
        print(toc-tic,"exe")
        if outfile != None:
            print ("Writing output file...")
            cat.write(outfile, out_format)
        return cat

    def estimate_moment_magnitude(self,output_file=None,density=2650.0,
                        vp=4266.96,
                        vsp_factor=1.84,time_before_pick=0.2,
                        time_after_pick=0.8,padding=20,
                        waterlevel=10.0,out_format="QUAKEML"):

        cat = self._estimate_moment_magnitudes(output_file,
                        vp,vsp_factor,density,time_before_pick,
                        time_after_pick,padding,
                        waterlevel,out_format)
        return cat


if __name__ =="__main__":
    from obspy.clients.fdsn import Client
    from obspy.core.inventory.inventory import read_inventory
    client = 'http://sismo.sgc.gov.co:8080'
    client = Client(client)

    resp = "/home/emmanuel/Ecopetrol/SeisMonitor/data/metadata/CM.dataless"
    p = Parser(resp)
    print(p)
    seed_id = "CM.TAPM.00.HHN"
    psi = p.get_paz(seed_id)
    print(psi)

    # resp = "/home/emmanuel/EDCT/SeisMonitor/data/events/public_CM.xml"
    # p = read_inventory(resp)
    # seed_id = "CM.BAR2.00.HHN"
    # psi = p.get_response(seed_id,"2022-05-01T00:00:00")
    # paz = psi.get_paz()
    # sensitivity = psi.instrument_sensitivity.value
    # gain = paz.normalization_factor
    # poles = paz.poles
    # zeros = paz.zeros

    # paz = {'poles': poles,
    #         'zeros': zeros,
    #         'gain': gain,
    #         'sensitivity': sensitivity}

    # st = client.get_waveforms("CM","BAR2","00","HHN",
    #                         UTCDateTime(2022,5,1,0,0,0),
    #                         UTCDateTime(2022,5,1,0,3,0))
    # tr = st[0]
    # print(tr.data)
    # x = tr.simulate(paz_remove=paz, water_level=10)
    # print(x.data)

    # catalog = "/home/emmanuel/EDCT/SeisMonitor/data/events/SGC2022knomqj.xml"
    # resp = "/home/emmanuel/Ecopetrol/SeisMonitor/data/metadata/CM.dataless"
    # # resp = "/home/emmanuel/EDCT/SeisMonitor/data/events/public_CM.xml"
    # out="./test_magnitude,xml"
    # mag = Magnitude(client,catalog,resp,
    #             jsonfile=None) 
    # mag.estimate_moment_magnitude(out)
    # print(mag)
