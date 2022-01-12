# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-21 12:03:16
#  * @modify date 2021-12-21 12:03:16
#  * @desc [description]
#  */

from obspy.clients.filesystem.sds import Client as SDS_Client
from obspy.core.utcdatetime import UTCDateTime
import argparse

DEFAULT_args = {'sds':"/home/emmanuel/archive/sds"}

def read_args():
    """
    read arguments from console
    """

    prefix = "+"

    parser = argparse.ArgumentParser("get seismic stream. ",
                                prefix_chars=prefix,
                                usage=f'get seismic stream.')

    parser.add_argument(prefix+"net",prefix*2+"network",
                        type=str,
                        metavar='',
                        help= "network | examples: 'EY'", 
                        required = True)
    
    parser.add_argument(prefix+"sta",prefix*2+"station",
                        type=str,
                        metavar='',
                        help= "station | examples: 'CA01A'",
                        required = True)

    parser.add_argument(prefix+"loc",prefix*2+"location",
                        type=str,
                        metavar='',
                        help= "location | examples: '00'",
                        required = True)
    
    parser.add_argument(prefix+"cha",prefix*2+"channel",
                        type=str,
                        metavar='',
                        help= "channel | examples: 'HHZ'",
                        required = True)

    parser.add_argument(prefix+"sds",prefix*2+"sds",
                        default=DEFAULT_args['sds'],
                        type=str,
                        metavar='',
                        help= "Root path of the seiscomp data structure",
    )
    
    parser.add_argument(prefix+"f",prefix*2+"filter",
                        action="store_true",
                        help= "Filter the data of all traces in the Stream.",
    )

    parser.add_argument(prefix+"n",prefix*2+"normalize",
                        action="store_true",
                        help= "Normalize all Traces in the Stream.",
    )

    parser.add_argument(prefix+"d",prefix*2+"detrend",
                        action="store_true",
                        help= "Remove a trend from the trace",
    )

    parser.add_argument(prefix+"p",prefix*2+"plot",
                        action="store_true",
                        help= "plot traces",
    )

    parser.add_argument(prefix+"st",prefix*2+"starttime",
                        type=str,
                        metavar='',
                        help="date | format: 'yyyymmddThhmmss'",
                         required = True)

    parser.add_argument(prefix+"et",prefix*2+"endtime",
                        type=str,
                        metavar='',
                        help="date | format: 'yyyymmddThhmmss'", 
                        required = True)

    args = parser.parse_args()
    args.starttime = UTCDateTime(args.starttime)
    args.endtime = UTCDateTime(args.endtime)
    vars_args = vars(args)
    return vars_args


def extraction(args):
    """
    args: dict
        -------example-------
        {'network': 'EY', 
        'station': 'CA05A', 
        'location': '00', 
        'channel': 'HHZ', 
        'sds': '/home/emmanuel/archive/sds',
        'filter': True, 
        'normalize': True, 
        'detrend': True,
        'plot': True,
        'starttime': UTCDateTime(2020, 11, 1, 3, 1, 58), 
        'endtime': UTCDateTime(2020, 11, 1, 3, 3)}

    returns: st
        Obspy stream
    """
    client = SDS_Client(args["sds"],
                        sds_type='D',
                         format='MSEED')
    wav_args = args.copy()
    if ('sds' in wav_args) or\
        ('detrend' in wav_args) or\
        ('filter' in wav_args) or\
        ('plot' in wav_args) or\
        ('normalize' in wav_args):   
        del wav_args['sds']
        del wav_args['detrend']
        del wav_args['filter']
        del wav_args['normalize' ]
        del wav_args['plot' ]

    st = client.get_waveforms(**wav_args)
    print(st)

    if args['detrend'] != False:
        st.detrend('demean')
    if args['filter'] != False:
        st.filter(type='bandpass', freqmin = 5.0, 
                freqmax = 45, corners=2, 
                zerophase=True)
    if args['normalize'] != False:
        st.normalize() 
    if args['plot'] != False:
        st.plot()
    return st



if __name__ == "__main__":
    args = read_args()
    print(args)
    extraction(args)