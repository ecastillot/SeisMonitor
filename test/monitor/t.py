from obspy.core.event.catalog import read_events


# cat = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/2016/20160101T000000__20160102T000000/associations/GaMMA/EQTransformer/associations.xml"
cat = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/2016/20160101T000000__20160102T000000/magnitudes/Ml/GaMMA/EQTransformer/Ml_magnitude.xml"

cat = read_events(cat)
for ev in cat:
    ev.origins[0].depth *= 1e3
    print(ev.origins[0].depth)
cat.write("/home/emmanuel/t.xml",format="SC3ML")