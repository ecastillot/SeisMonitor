import os
from obspy.core.event.catalog import Catalog, read_events
from SeisMonitor.utils import printlog, isfile
import utils as ut

class Hypocenter():
    def __init__(self,catalog,xml_path,vel_path,out_dir):

        if isinstance(catalog,Catalog):
            self.catalog = catalog
        else:
            self.catalog = read_events(catalog)
        
        self.xml_path = xml_path
        self.vel_path = vel_path
        self.out_dir = out_dir
        self.out_file = os.path.join(self.out_dir,"hypocenter","hypocenter_events.xml")
        self.paths = os.path.join(self.out_dir,"hypocenter","hypocenter_paths")
        self.sta0 = os.path.join(self.paths,"STATION0.HYP")

    def locate(self,sfilename = "catalog.sfile",
                out_format="QUAKEML",
                rm_not_locatable=True):
        sfile = os.path.join(self.paths,sfilename)
        isfile(sfile)
        self.catalog.write(sfile,format="NORDIC")

        sta0 = ut.STATION0(self.xml_path,self.vel_path)
        sta0.write(self.sta0)

        hyp = ut.HypocenterTools(self.paths,sfilename)
        hyp.split()
        hyp.remodl_and_setbrn()
        hyp.update()

        ut.check_sfile_integrity(self.paths,rm_not_locatable)
        catalog = hyp.collect()
        if self.out_file != None:
            print ("Writing output file...")
            isfile(self.out_file)
            catalog.write(self.out_file, out_format)
        return catalog