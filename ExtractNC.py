# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 20:02:03 2023

@author: toru1
"""
import netCDF4
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import os
import pandas as pd


ncfile = r'RAIN4PE_daily_0.1d_1981_2015_v1.0.nc'
dnc= xr.open_dataset(ncfile)
keys = list(dnc.variables)
LAT = keys[0]
LON = keys[1]
TIME = keys[2]
PCP = keys[3]

lat = dnc.variables[LAT][:]
lon = dnc.variables[LON][:]
time = dnc[TIME].to_index()
time.name = 'date'
da = dnc[PCP]

df = pd.DataFrame(index=time)
Coord_stations=pd.read_csv('csv/WGEN_stat.csv')

for _, row in Coord_stations.iterrows():
    da2=da.sel({
        LAT:row['lon'],
        LON:row['lat']},
        method='nearest')
    df[row['name']] = da2.to_dataframe()[keys[3]]

df.index = df.index.strftime('%d/%m/%Y')
df.to_csv('pcp.csv')