import sys
import os
seismopath = "/home/emmanuel/EDCT"
seismonitor = os.path.join(seismopath,"SeisMonitor")
sys.path.insert(0,seismonitor)

from obspy.clients.fdsn import Client
from obspy.core.inventory.inventory import read_inventory
from SeisMonitor.monitor.magnitude.mag import Magnitude

client = 'http://sismo.sgc.gov.co:8080'
client = Client(client)
catalog = "/home/emmanuel/EDCT/test/associations/associations.xml"

resp = "/home/emmanuel/EDCT/SeisMonitor/data/events/public_CM.xml"
resp = read_inventory(resp)

out = "/home/emmanuel/EDCT/test"

mag = Magnitude(client,catalog,resp,out) 
cat = mag.get_Ml(mag_type="RSNC",trimmedtime=5,out_format="SC3ML")
print(cat)