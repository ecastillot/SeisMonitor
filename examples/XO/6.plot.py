from obspy import read_events

catalog_path = "/home/emmanuel/XO_monitor_results/magnitude/nlloc/Ml/Ml_magnitude.xml"
catalog = read_events(catalog_path)
print(catalog)


### pip install Cartopy first
catalog.plot()
catalog.plot(projection="local")

