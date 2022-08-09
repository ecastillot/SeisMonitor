#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 20:00:00 2020
@author: Emmanuel_Castillo
last update: 14-11-2020 
"""


from obspy import UTCDateTime
from obspy.clients.fdsn import Client as FDSNClient

sgc_client = FDSNClient('http://10.100.100.232:8091')
st = sgc_client.get_waveforms(network="CM",
                    station="*",
                    location="*",
                    channel="*",
                    starttime=UTCDateTime("2017-12-24T19:00:00.000000Z"),
                    endtime=UTCDateTime("2017-12-24T19:01:00.000000Z"))
print(st)