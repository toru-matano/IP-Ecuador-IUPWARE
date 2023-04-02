# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 21:45:46 2023

@author: toru1
"""
from pyETo import NetEvapoTrans
import numpy as np
import pandas as pd

METHODS = ['M1', 'M2', 'M3', 'M4']

#   End of DOY: Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec
NOY = np.array([31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 335, 365])
lat = -3.62   #degree N:+, S:-
method = METHODS[0]

df_pcp = pd.read_csv('pcp.csv', index_col=0)
df_tmax = pd.read_csv('tmp_max.csv', index_col=0)
df_tmin = pd.read_csv('tmp_min.csv', index_col=0)
df_tmean = (df_tmax + df_tmin)/2

df_pcp.index = pd.to_datetime(df_pcp.index, format='%d/%m/%Y')
df_tmax.index = pd.to_datetime(df_tmax.index, format='%d/%m/%Y')
df_tmin.index = pd.to_datetime(df_tmin.index, format='%d/%m/%Y')
df_tmean.index = pd.to_datetime(df_tmean.index, format='%d/%m/%Y')

NumYear = df_tmean.index[-1].year - df_tmean.index[0].year +1

ETO_class = NetEvapoTrans()
Ra = ETO_class.extraterrestrial_radiation(lat, NOY)
Ra = np.repeat(Ra, NumYear)
Stations = df_tmean.columns

Pcp = df_pcp.resample('M').mean()
Tmean = df_tmean.resample('M').mean()
Tmax = df_tmax.resample('M').mean()
Tmin = df_tmin.resample('M').mean()


ETo_all = {}
for st in Stations:
    try:
        ETo = pd.DataFrame()
        for method in METHODS:
            eto = ETO_class.GetETo(method, Tmean[st], Tmax[st], Tmin[st], Ra, Pcp=Pcp[st])
            ETo[method] = eto
        #    eto = [ETO_class.GetETo(method, Tmean, Tmax, Tmin, Ra) for method in METHODS]
        #    eto = pd.concat(eto, axis=1)
        ETo.index = ETo.index.strftime('%m/%Y')
        ETo_all[st] = ETo
    except:
        print('Missing data at St. ', st)

