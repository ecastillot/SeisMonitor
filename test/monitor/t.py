from obspy.core.event.catalog import read_events
import os
import glob
import shutil
# # cat = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/2016/20160101T000000__20160102T000000/associations/GaMMA/EQTransformer/associations.xml"
# cat = "/media/emmanuel/TOSHIBA EXT/ColSeismicity/2016/20160101T000000__20160102T000000/magnitudes/Ml/GaMMA/EQTransformer/Ml_magnitude.xml"

# cat = read_events(cat)
# for ev in cat:
#     ev.origins[0].depth *= 1e3
#     print(ev.origins[0].depth)
# cat.write("/home/emmanuel/t.xml",format="SC3ML")


folder = "/home/emmanuel/E_ColSeismicity/ColSeismicity/2018"

# for x in glob.glob(os.path.join(folder,"**","NLLOC"),recursive=True):
#     print(x)

dfs = []
for dp, dn, filenames in os.walk(folder):
    for f in filenames:
        # search_path = os.path.join(dp, f)
        if "NLLOC" in dp:
            shutil.rmtree(dp)
            print(dp)
