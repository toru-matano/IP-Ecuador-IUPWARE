# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 15:53:05 2023

@author: toru1
"""
import numpy as np

class NetEvapoTrans:
    def __init__(self):
        pass
    
    def extraterrestrial_radiation(self, lat, noy):
        '''

        Parameters
        ----------
        noy : int
            number of year (0 to 365 or 366)
        lat : float
            latitude
        
        omega_s: sunset hour angle
        dr: inverse relative distance Earth-Sun
        sigma: solar decimation

        '''
        lat = lat*np.pi/180
        solar_const = 0.0820
        dr = 1 + 0.033*np.cos(2*np.pi/365*noy)
        sigma = 0.409*np.sin(2*np.pi/365*noy-1.39)
        omega_s = np.arccos(-np.tan(lat)*np.tan(sigma))
        gamma = (omega_s*np.sin(sigma)*np.sin(lat) + np.cos(sigma)*np.cos(lat)*np.sin(omega_s))
        
        Ra = 24*60*solar_const*dr*gamma/np.pi
        
        return Ra
    
    def _Hargreaves_M1(self, Tmean, Tmax, Tmin, Ra):
        ETo = 0.408*0.0030*(Tmean+20)*(Tmax-Tmin)**0.4*Ra
        return ETo

    def _Hargreaves_M2(self, Tmean, Tmax, Tmin, Ra):
        ETo = 0.408*0.0025*(Tmean+16.8)*(Tmax-Tmin)**0.5*Ra
        return ETo

    def _Hargreaves_M3(self, Tmean, Tmax, Tmin, Ra, Pcp):
        ETo = 0.408*0.0013*(Tmean+17)*(Tmax-Tmin-0.0123*Pcp)**0.76*Ra
        return ETo

    def _Hargreaves_M4(self, Tmean, Tmax, Tmin, Ra):
        ETo = 0.408*0.0023*(Tmean+17.8)*(Tmax-Tmin)**0.424*Ra
        return ETo
        
    def _JensenHaise(self, Tmean, Rs, LH):
        Tx = -3
        ETo = 0.025*(Tmean-Tx)*Rs/LH
        return ETo

    def _McGuinnessBordne(self, Tmean, Rs):
        ETo = (0.0082*Tmean - 0.19)*Rs/1500*2.54
        return ETo
    
    def _PenmanMonteith(self, z, Tmean, Tmax, Tmin, Rn, U):
        CPs = 0.14
        P = 101.3*(293-0.0065*z/293)**5.26
        gamma = 0.665*10**-3*P
        
        def saturation_vapour_pressure(T):
            return 0.6108*np.exp(17.27*T/(T+237.3))
        es = (saturation_vapour_pressure(Tmax) + saturation_vapour_pressure(Tmin))/2
        ea = saturation_vapour_pressure(Tmin)
        
        delta = 4098*saturation_vapour_pressure(Tmean)/(Tmean+273.3)**2
        
        G = CPs*(Tmean)

        ETo = (0.408*delta*(Rn-G)+gamma*900/(Tmean+273)*U*(es-ea))/(delta+gamma*(1+0.34*U))
        
        return ETo
    
    def GetETo(self, method, Tmean, Tmax, Tmin, Ra, **kwarg):
        try:
            Pcp = kwarg['Pcp']
            LH = kwarg['LH']
            Rs = kwarg['Rs']
        except:
            pass
        
        if method=='M1':
            return self._Hargreaves_M1(Tmean, Tmax, Tmin, Ra)
        elif method == 'M2':
            return  self._Hargreaves_M2(Tmean, Tmax, Tmin, Ra)
        elif method == 'M3':
            return  self._Hargreaves_M3(Tmean, Tmax, Tmin, Ra, Pcp)
        elif method == 'M4':
            return  self._Hargreaves_M4(Tmean, Tmax, Tmin, Ra)
#        elif method == 'Jensen':
#            return self._JensenHaise(Tmean, Rs, LH)
#        elif method == 'McGuinness':
#            return self._McGuinnessBordne(Tmean, Rs)
