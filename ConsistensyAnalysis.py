# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 21:08:31 2023

@author: toru1
"""

import pandas as pd
import pymannkendall as pmk 
# install package: conda install -c conda-forge pymannkendall
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import os
import glob

def homogeneity_MannKendall(rainseries, alpha):
    # https://www.statology.org/mann-kendall-test-python/
    # significance = [0.01,0.05,0.1] possible alpha values
    test = pmk.original_test(rainseries, alpha)
    
    return test
    
def homogeneity_Buishand(rainseries, title):
    # Buishand, 1982
    mean=np.nanmean(rainseries)
    std=np.nanstd(rainseries)
    totcount = rainseries.size
    cum_dev=np.zeros(rainseries.size)
    
    for i in range(len(rainseries.index)):
        cum_dev[i]=np.nansum(rainseries[0:i+1]-mean)
            
    normal_Q = cum_dev.max()/(std*np.sqrt(totcount))
    normal_R = (cum_dev.max()-cum_dev.min())/(std*np.sqrt(totcount))
    # print('Q/n^.5: {}\nR/n^.5: {}'.format(normal_Q, normal_R))
    
    return cum_dev, totcount, mean, std
    
# %% main

def main(path):
    if path[-3:] == 'pcp':     
        ds = np.loadtxt(path, dtype=float, skiprows=3).T
    elif path[-3:] == 'csv':
        ds = np.loadtxt(path, delimiter=',')
    title = path.split('\\')[-1][:-4]
    
    # convert YEAR & DAY OF YEAR to datetime
    index = pd.to_datetime(range(ds.shape[1]), unit='D', origin=str(int(ds[0, 0])))
    rainseries = pd.DataFrame(ds[2], index=index)
    rainseries.replace(np.nan, 0, inplace=True)

    style.use('ggplot')
    
    # # rainfall time series
    # plt.figure()
    # plt.scatter(rainseries.index, rainseries, s=2, marker='o')
    # plt.ylim(bottom=0)
    # plt.xlabel("Year")
    # plt.title(title)
    # plt.show()
    
    # homogeneity
    test = homogeneity_MannKendall(rainseries, alpha=0.01)
    cum_dev, totcount, mean, std = homogeneity_Buishand(rainseries, title)
    print('{}\t{}'.format(title, test.p))


    reject_ratio = {
        #when sample size n >> 20
        '90%': 1.22,
        '95%': 1.36,
        '99%': 1.63,            
        }
    
    fig, ax = plt.subplots(1, 1)
    ax.plot(rainseries.index, cum_dev/std, linestyle='--')
    ax.axhline(y=0, color='k', linewidth=0.5)
    for key, ratio in reject_ratio.items():
        ax.axhline(y=ratio*np.sqrt(totcount), color='k', linewidth=0.5)
        ax.axhline(y=-ratio*np.sqrt(totcount), color='k', linewidth=0.5)
        ax.text(x=rainseries.index[10], y=ratio*np.sqrt(totcount), s=key)
    
    #plt.text(x=rainseries.index[10], y=rainseries.max(), s=test.trend)

    ax.set_xlabel("Year")
    ax.set_ylabel('Rescaled cumulative deviation')
    yabs_max = abs(max(ax.get_ylim(), key=abs))
    ax.set_ylim(ymin=-yabs_max, ymax=yabs_max)
    ax.set_title(title)


if __name__ == '__main__':
    folder = r'C:\Users\toru1\OneDrive - KU Leuven\Y2-S2\IP Ecuador\datasets_SWAT+\weather_swat_plus\SWAT_generator'
#    folder = r'C:/Users/toru1/OneDrive - KU Leuven/Y2-S2/IP Ecuador/SWAT inputs IP/example/example/example/Scenarios/Default/TxtInOut'
    PATH = glob.glob(folder+'\\*.pcp')
    for path in PATH:
        try:
            main(path)
        except:
            print('\nskipped: ', path)
            
