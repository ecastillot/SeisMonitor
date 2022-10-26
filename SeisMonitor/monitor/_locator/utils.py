#!/home/dsiervo/anaconda3/envs/pnet/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mar 25 11:25:33 2020

@author: Daniel Siervo, emetdan@gmail.com
"""
import csv
from obspy import UTCDateTime
import datetime
import os
import numpy as np
import pandas as pd
import glob
import sys

def prepare_eqt(df):
    """

    Parameters
    ----------
    df: DataFrame
        Dataframe that contains EQTransformer outputs
    min_prob : float
        Mimimum probability to consider a pick
    """
    interest_col = ['network', 'station', 'instrument_type', 'detection_probability',
                    'p_arrival_time', 'p_probability', 'p_snr',
                    's_arrival_time', 's_probability', 's_snr']
    # Keeping with columns of interest
    df = df[interest_col]
    # Deleting rows without P phase
    df = df.dropna(subset=['p_probability'])
    # Changing NaN for "no pick"
    df = df.fillna('no pick')
    # Removing white spaces arround station name
    df['station'] = df['station'].str.strip()

    # Create a list of picks
    p_list = read_eqt_picks(df['network'].tolist(), df['station'].tolist(),
                            df['instrument_type'].tolist(),
                            df['p_arrival_time'].tolist(),
                            df['p_probability'].tolist(),
                            df['s_arrival_time'].tolist(),
                            df['s_probability'].tolist(),
                            p_snrs=df['p_snr'].tolist(),
                            s_snrs=df['s_snr'].tolist(),
                            ev_probs=df['detection_probability'].tolist())

    return p_list

def read_eqt_picks(nets, stations, chs, p_times, p_probs, s_times, s_probs,
                   p_snrs, s_snrs,ev_probs):
    """Iterate over picks to generate list of Picks objects

    Parameters
    ----------
    nets : list
        Networks list
    stations : list
        Stations list
    chs : list
        Channels list
    p_times : list
        P pick times
    p_probs : list
        P picks probabilities
    s_times : list
        S picks times
    s_probs : list
        S picks probabilities
    """
    picks_list = []
    for i in range(len(s_times)):
        net, station, ch = nets[i], stations[i], chs[i]
        s_pick = None
        loc = '00'
        if ch == 'EH':
            loc = '20'
        elif ch == 'HN':
            loc = '10'

        ch += 'Z'

        ev_prob = ev_probs[i]
        p_t, p_prob, p_snr  = p_times[i], p_probs[i], p_snrs[i]

        # filtering bad picks according to statistical analysis doing in
        # the notebook EDA_EQTransformer_probs.ipynb
        if p_snr == 'no pick':
            p_snr = 0
        p_snr = float(p_snr)

        # if p_prob >= 0.01 and p_snr >= -1.7:
        p_pick = eqt_pick_constructor(p_t, p_prob, net,
                                        station, loc, ch, 'P',p_snr,ev_prob)
        picks_list.append(p_pick)

        s_t, s_prob, s_snr = s_times[i], s_probs[i], s_snrs[i]
        if s_t != 'no pick':
            if s_snr == 'no pick':
                s_snr = 0
            s_snr = float(s_snr)

            # if float(s_prob) >= 0.01 and s_snr >= 0:
            s_pick = eqt_pick_constructor(s_t, s_prob, net,
                                            station, loc, ch, 'S',s_snr,ev_prob)
            picks_list.append(s_pick)

        # print(f'{i}. {station}, p_t:{p_t}, p_prob:{p_prob}, s_t:{s_t}, s_prob:{s_prob}, p_snr:{p_snr}, s_snr:{s_snr}')

    return picks_list

def picks2xml(pick_list):
    """Transform pick list into an Seiscomp3 XML file

    Parameters
    ----------
    pick_list : list
        List with Pick objects
    
    Returns
    -------
    str
        Returns the string of the XML file
    """
    xml_top = '''<?xml version="1.0" encoding="UTF-8"?>
<seiscomp xmlns="http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/0.10" version="0.10">
  <EventParameters>'''

    xml_bottom = '''
  </EventParameters>
</seiscomp>'''

    xml_S_bottom = '''<comment>
        <text>{}</text>
        <id>RefPickID</id>
      </comment>
    </pick>'''

    xml_file = xml_top

    pick_dic = {}
    for pick in pick_list:
        if pick.phaseHint == 'S':
            try:
                P_pick = pick_dic[pick.station+'_'+'P']
            except KeyError:
                print('No se encontro fase P de referencia para:')
                print(f'\t{pick.publicID}')
                continue
            xml_file += pick.toxml()+xml_S_bottom.format(P_pick.publicID)
        elif pick.phaseHint == 'P':
            # se crea un diccionario para almacenar las fases P a las que
            # luego se relacionaran las fases S. (Seiscomp lo exige)
            pick_dic[pick.station+'_'+pick.phaseHint] = pick 
            xml_file += pick.toxml()
    
    xml_file += xml_bottom

    return xml_file

def id_maker(pick_time, net, station, loc, ch, phaseHint, ai_type):
    """Creates the seiscomp Pick PublicID
    
    Parameters
    ----------
    pick_time : Obspy UTCDateTime object
        Time for phase pick
    
    Returns
    ------
    str
       Seiscomp pick PublicID
    """
    dateID = pick_time.strftime('%Y%m%d.%H%M%S.%f')[:-4]
    if phaseHint == 'P':
        publicID = dateID+f'-{ai_type}-{net}.{station}.{loc}.{ch}'
    elif phaseHint == 'S':
        publicID = dateID+f'-{ai_type}-{net}.{station}.{loc}.{ch}'
    return publicID

def eqt_pick_constructor(time, prob, net, station, loc, ch, ph,snr,ev_prob):
    time = UTCDateTime(time)
    id_ = id_maker(time, net, station, loc, ch, ph, 'EQTransformer')
    creation_t = UTCDateTime()

    evaluation = 'automatic'
    if prob >= 0.95:
        evaluation = 'manual'

    pick = Pick(id_, time, net, station, loc, ch, prob,
                ph, creation_t, evaluation, 'EQTransformer',snr,ev_prob)

    return pick


class Pick:
    '''Class to represent a phase pick

    Atributes
    ---------
    xml_P_block : str
        Template of a pick XML block

    Methods
    -------
    toxml()
        Create seiscomp3 xml block
    '''
    xml_P_block = '''
    <pick publicID="{publicID}">
      <time>
        <value>{pick_time}</value>
      </time>
      <waveformID networkCode="{net}" stationCode="{station}" locationCode="{loc}" channelCode="{ch}"/>
      <filterID>Probability_{prob}+snr_{snr}+eventprob_{ev_prob}</filterID>
      <methodID>{author}</methodID>
      <phaseHint>{phaseHint}</phaseHint>
      <evaluationMode>{evaluation}</evaluationMode>
      <creationInfo>
        <agencyID>SGC</agencyID>
        <author>{author}</author>
        <creationTime>{creation_time}</creationTime>
      </creationInfo>
    </pick>'''

    xml_S_block = '''
    <pick publicID="{publicID}">
      <time>
        <value>{pick_time}</value>
      </time>
      <waveformID networkCode="{net}" stationCode="{station}" locationCode="{loc}" channelCode="{ch}"/>
      <filterID>Probability_{prob}+snr_{snr}+eventprob_{ev_prob}</filterID>
      <methodID>{author}</methodID>
      <phaseHint>{phaseHint}</phaseHint>
      <evaluationMode>{evaluation}</evaluationMode>
      <creationInfo>
        <agencyID>SGC</agencyID>
        <author>{author}</author>
        <creationTime>{creation_time}</creationTime>
      </creationInfo>
      '''

    def __init__(self, publicID, pick_time,
                 net, station, loc, ch, prob,
                 phaseHint, creation_time,
                 evaluation, author='PhaseNet',snr=0,
                 ev_prob=0):
        """
        Parameters
        ----------
        publicID : str
            Seiscomp3 pick pubicID
        pick_time : Obspy UTCDateTime object
            Time for phase pick.
        net : str
            Agency network.
        station : str
            Station name.
        loc : str
            Location code.
        ch : str
            Channel
        phaseHint : str
            Phase type, can be P or S.
        creation_time : Obspy UTCDateTime object
            Time for creation time.
        """
        self.publicID = publicID
        self.pick_time = pick_time
        self.net = net
        self.station = station
        self.loc = loc
        self.ch = ch
        self.prob = prob
        self.ev_prob = ev_prob
        self.snr = snr
        self.phaseHint = phaseHint
        self.creation_time = creation_time
        self.evaluation = evaluation
        self.author = author
    
    def toxml(self):
        """Create a seiscomp xml block
        """

        if self.phaseHint == 'P':
            return self.xml_P_block.format(
                publicID=self.publicID,
                pick_time=self.pick_time,
                net=self.net,
                station=self.station,
                loc=self.loc,
                ch=self.ch,
                prob=self.prob,
                snr=self.snr,
                ev_prob=self.ev_prob,
                phaseHint=self.phaseHint,
                creation_time=self.creation_time,
                evaluation=self.evaluation,
                author=self.author
                )
        elif self.phaseHint == 'S':
            return self.xml_S_block.format(
                publicID = self.publicID,
                pick_time = self.pick_time,
                net = self.net,
                station = self.station,
                loc = self.loc,
                ch = self.ch,
                prob = self.prob,
                snr=self.snr,
                ev_prob=self.ev_prob,
                phaseHint = self.phaseHint,
                creation_time = self.creation_time,
                evaluation = self.evaluation,
                author = self.author
                )

if __name__=='__main__':
    csv_file = "/home/ecastillo/repositories/aipicker_stats/juldaypicks/CM_2020-01/001/eqt_picks.csv"
    df = pd.read_csv(csv_file)
    list_picks = prepare_eqt(df)

    print(list_picks[0].toxml() )
    
    # import sys

    # if len(sys.argv) == 2:
    #     main_picks(sys.argv[1])
    # elif len(sys.argv) == 3:
    #     main_picks(sys.argv[1], sys.argv[2])
    # elif len(sys.argv) == 4:
    #     main_picks(sys.argv[1], sys.argv[2], ai='eqt')
    # else:
    #     main_picks()

