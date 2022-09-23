import shutil
import numpy as np
import pandas as pd
import os
import json 
import numpy as np
import matplotlib.pyplot as plt

def get_d_az(p1,p2):
    """
    Parameters:
    -----------
    p1: tuple
        lon, lat
    p2: tuple
        lon,lat
    
    Returns:
    --------
    Get distance and azimuth 
    """

    d2 = (p2[1]-p1[1])**2 + (p2[0]-p1[0])**2
    d = np.sqrt(d2)

    if (p2[0]-p1[0])==0:
        theta = 90
    else:
        theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) * 180/np.pi

    return d,theta

def get_t_points(p1,p2,d,d2km=114):
    """
    Parameters:
    -----------
    p1 : tuple
        lon,lat 
    p2 : tuple
        lon, lat
    d : float
        distance in km of the transect
    d2km: float
        conversor degrees to kilometers
    """

    _,theta = get_d_az(p1,p2)
    alpha = 90 - theta

    x = d*np.cos(alpha * np.pi/180)
    y = d*np.sin(alpha * np.pi/180)
    x= x/d2km
    y= y/d2km

    pm = ((p2[0]+p1[0])/2 , (p2[1]+p1[1])/2)

    tp1 = (pm[0]-x ,pm[1]+y)
    tp2 = (pm[0]+x,pm[1]-y)

    return tp1, tp2

def get_centers(N,p1,p2):
    
    """
    Paramters:
    ----------
    N: int
        Number of divisions (transects)
    p1: tuple
        (lon, lat) order.
    p2: tuple
        (lon, lat) order.
    Return:
    -------
    center: np.array
        arra of points that indicate the center
    azi: float 
        azimuth between p1 and p2
    """
    d2 = (p2[1]-p1[1])**2 + (p2[0]-p1[0])**2
    d = np.sqrt(d2)
    if (p2[0]-p1[0])==0:
        theta = 90
    else:
        theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) * 180/np.pi

    d = np.linspace(0,d,num=N,endpoint=True)
    x = p1[0] + d*np.cos(theta*np.pi/180)
    y = p1[1] + d*np.sin(theta*np.pi/180)
    center = np.array(list(zip(x,y)))

    azi = 90 - theta

    return center, azi

def get_line_in_map(center,distance,azimuth,d2k=114,save=None):
    """
    Parameters:
    -----------
    center: tuple
        (lon, lat) order.  
    distance: tuple
        (left_distance,rigth_distance) order in km.
    azimuth: float
        degrees
    d2k: float
        factor of conversion kilometers to degree
    save: str
        Path
    """
    cx,cy = center 
    dl,dr = distance
    azi = azimuth * np.pi/180

    xr = cx + np.sin(azi)*dr/d2k
    yr = cy + np.cos(azi)*dr/d2k
    xl = cx - np.sin(azi)*dl/d2k
    yl = cy - np.cos(azi)*dl/d2k

    alpha =  np.arctan((yr-yl)/(xr-xl) )* 180/np.pi
    # alpha = 90 -alpha
    # print(xr,yr,xl,yl,alpha)

    if save != None:
        df = pd.DataFrame({'lon': [xl, xr], 'lat': [yl, yr]})
        df.to_csv(save,sep = " ", header=False,index=False)
    
    return (xr,yr,xl,yl,alpha)

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def get_df(csv, between=None,
            columns=["longitude","latitude","depth"]):

    events = pd.read_csv(csv)
    events["time"] = events["time_event"]
    events["time"] = pd.to_datetime(events["time"])
    if between != None:
        events = events[(events["time"]>between[0]) & (events["time"]<between[1]) ]

    events = events[columns]
    return events

if __name__ == "__main__":
    x = get_line_in_map((-76.45,1.542),(5,5),35,d2k=114,save=None)
    print(x)