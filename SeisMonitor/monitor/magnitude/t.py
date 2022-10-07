from obspy import read_events

cat = read_events("/home/emmanuel/nlloc_iter/mag_out/Ml_magnitude.xml")
print(cat[0].station_magnitudes )