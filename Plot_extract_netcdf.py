# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 23:26:30 2023

@author: ALEXIS
"""

###OPEN NETCDF FILES FOR THE REGIONAL CLIMATIC MODELS####
### EXTRACT INFORMATION FROM RCM SAM Rotated coordinates###

import netCDF4
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import os
import pandas as pd
from mpl_toolkits.basemap import Basemap

directory=os.getcwd()
print(directory)
## open netcdf file, define directorys
path=(r"E:\KULEUVEN\Tesis\RCMs\Merged_SA\RCMs_Bolivia")
pathcsv=(r"D:\KU LEUVEN\Courses\Thesis\08_Planillas")
Pathout=(r"D:\KU LEUVEN\Courses\Thesis\10_RCMs\0_Extracted files")
##Here change the file name
os.chdir(path)
RCM_name='15_pr_historical_RCA4_15'
data= xr.open_dataset('15_pr_historical_RCA4_15_Bolivia_mmd.nc')

#print (data)
print(data.data_vars)
#print(data.variables.keys())
# Storing the lat and lon data into the variables 
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]

### Storing the lat and lon of station into variables 
os.chdir(pathcsv)
Coord_stations=pd.read_csv('08_Rotated_Coords.csv')
cds=Coord_stations

##Rotated coordinates
lat_station = -0.76
lon_station =  172.99

##Variables
da=data['pr']
pt=data.isel(time=1000)
#pt.plot()
#In case the precipitation is not in mm/d
#pr_RCM=data['pr']*86400
date_range=data['time'].to_index()
df=pd.DataFrame(index=date_range)
os.chdir(Pathout)
for i in range(len(cds)):
    lt_i=cds.loc[i,'rlat']
    ln_i=cds.loc[i,'rlon']
    name=cds.loc[i,'Name']
    da2=da.sel(rlat=lt_i,rlon=ln_i,method='nearest') ###Extractor of information
    da2=da2.to_dataframe()
    df[name]=da2['pr']
    #df.to_csv('0_'+str(name)+'RCM_0'+'.csv')

df.to_csv(str(RCM_name)+'_All_stations'+'.csv',encoding = 'utf-8-sig')    
#df[name]=da2['pr']
#print(df)

