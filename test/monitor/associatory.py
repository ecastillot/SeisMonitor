import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from SeisMonitor.monitor.associator import utils as ut
from obspy.core.utcdatetime import UTCDateTime

sp = "/home/emmanuel/EDCT/test/picks/picks/seismonitor_picks.csv"
eqt_f = "/home/emmanuel/EDCT/test/picks/picks/eqt_fmt"
ut.seismonitor_picks_to_eqt_fmt(sp,eqt_f,
                            tt=30,st=2,et=4)