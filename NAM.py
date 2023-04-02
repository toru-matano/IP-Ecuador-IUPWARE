# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 07:54:15 2023

@author: toru1
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats

class NAM:
    def __init__(self, Start='1900-01-01', End='1900-01-01', EndWU='1900-01-01', time_step='H', area=1):
        
        self.Start = pd.Timestamp(Start)
        self.End = pd.Timestamp(End)
        self.EndWU = pd.Timestamp(EndWU)
        self.time_step = time_step
        
        self.sim_time = pd.date_range(self.Start, self.End, freq=self.time_step)

        self.QOF = np.zeros(self.sim_time.size)
        self.OF_r = np.zeros(self.sim_time.size)        # Overlandflow
        self.IF_r = np.zeros(self.sim_time.size)        # Interflow
        self.G1 = np.zeros(self.sim_time.size)
        self.BF_r = np.zeros(self.sim_time.size)        # Baseflow
        self.QIF = np.zeros(self.sim_time.size)
        self.Qt = np.zeros(self.sim_time.size)
        self.TSK1 = np.zeros(self.sim_time.size)
        self.TSK2 = np.zeros(self.sim_time.size)
        self.S = np.zeros(self.sim_time.size)           # Snow pack
        self.ET = np.zeros(self.sim_time.size)          # Evapotranspiration

        self.sub_area = area                   #km2
        
        self.flow = np.zeros(self.sim_time.size)
        self.pet = np.zeros(self.sim_time.size)
        self.precip = np.zeros(self.sim_time.size)
        self.temp = np.zeros(self.sim_time.size)
        
        self.relU = 0
        self.relL = 0
        self.OFmin = 0
        self.Beta = 0
        self.IF = 0
        self.IFp = 0
        self.BF = 0
        self.Ss = 0
        self.Sw = 0
        self.OFp = 0
        self.OF = 0
        self.U = 0
        self.L = 0
        
        if time_step == 'H':
            self.fA = 3.6                           # L/h (mm/h) --> m3/s  (3600/1000)
        elif time_step == 'D':
            self.fA = 3.6*24
        elif time_step == 'M' or 'MS':
            self.fA = 3.6*24*30
        
        self.snow = False
        
        #               Umax, Lmax, CQOF, CKIF, CK1, CK2, TOF, TIF, TG, CKBF, Csnow, T0
        self.para_name = ['Umax', 'Lmax', 'CQOF', 'CKIF', 'CK1', 'CK2', 'TOF', 'TIF', 'TG', 'CKBF', 'Csnow', 'T0']
        self.paramMin = [ 5,   50, 0.0,  150,  3,  3,   0,   0,    0,  500, 2, -3]
        self.paramMax = [35, 1000, 1.0, 2000, 48, 48, 0.9, 0.9, 0.99, 7500, 5,  3]


    def input_timeseries(self, flow_df, pet_df, precip_df, temp_df=None):

        start_flow = pd.Timestamp(flow_df['Year'][0],flow_df['Month'][0],flow_df['Day'][0], flow_df['Hour'][0])
        end_flow = pd.Timestamp(flow_df['Year'].iloc[-1],flow_df['Month'].iloc[-1],flow_df['Day'].iloc[-1], flow_df['Hour'].iloc[-1])
        start_pet = pd.Timestamp(pet_df['Year'][0],pet_df['Month'][0],pet_df['Day'][0], pet_df['Hour'][0])
        end_pet = pd.Timestamp(pet_df['Year'].iloc[-1],pet_df['Month'].iloc[-1],pet_df['Day'].iloc[-1], pet_df['Hour'].iloc[-1])
        start_precip = pd.Timestamp(precip_df['Year'][0],precip_df['Month'][0],precip_df['Day'][0], precip_df['Hour'][0])
        end_precip = pd.Timestamp(precip_df['Year'].iloc[-1],precip_df['Month'].iloc[-1],precip_df['Day'].iloc[-1], precip_df['Hour'].iloc[-1])
        
        time_flow = pd.date_range(start_flow, end_flow, freq=self.time_step)
        time_pet = pd.date_range(start_pet, end_pet, freq=self.time_step)
        time_precip = pd.date_range(start_precip, end_precip, freq=self.time_step)
        
        flow_df = pd.DataFrame(flow_df.values[:, 4:], index=time_flow)
        pet_df = pd.DataFrame(pet_df.values[:, 4:], index=time_pet)
        precip_df = pd.DataFrame(precip_df.values[:, 4:], index=time_precip)

        if self.Start < max(start_flow, start_pet, start_precip):
            self.Start = max(start_flow, start_pet, start_precip)
            warnings.warn('The selected start date is before the first input record.' /
                  'Setting the start date to the firt date where all data is available: {}\n'.format(
                      self.Start), Warning)
            
        if self.End > min(end_flow, end_pet, end_precip):
            self.End = min(end_flow, end_pet, end_precip)
            warnings.warn('The selected end date is after the last input record' /
                  'Setting the end date to the last date where all data is available: {}\n'.format(
                      self.End), Warning)    

        flow_df = flow_df.loc[self.sim_time]
        flow_df = flow_df.where(flow_df>=0)
        self.flow = flow_df.values.flatten()
        pet_df = pet_df.loc[self.sim_time]
        pet_df = pet_df.where(pet_df>=0)
        self.pet = pet_df.values.flatten()
        precip_df = precip_df.loc[self.sim_time]
        precip_df = precip_df.where(precip_df>=0)
        self.precip = precip_df.values.flatten()

        if self.snow:
            try:
                start_temp = pd.Timestamp(temp_df['Year'][0],temp_df['Month'][0],temp_df['Day'][0], temp_df['Hour'][0])
                end_temp = pd.Timestamp(temp_df['Year'].iloc[-1],temp_df['Month'].iloc[-1],temp_df['Day'].iloc[-1], temp_df['Hour'].iloc[-1])
                time_temp = pd.date_range(start_temp,end_temp, freq = self.time_step)
                temp_df = pd.DataFrame(data=temp_df.iloc[:,4:].values, index=time_temp)
                temp_df = temp_df.loc[self.sim_time]
                self.temp = temp_df.values.flatten()
            except:
                warnings.warn("Temperature file not provided. Snow melt cannot be simulated", Warning)
                self.snow = False

    
    def set_param(self, Fpar, index_Fpar):
        # calibration parameters
        self.index_Fpar = index_Fpar
        
        # check if parameters are within range
        for elem in list(zip(self.paramMin, Fpar, self.paramMax, self.para_name)):
            if elem[0]<=elem[1]<=elem[2]:
                self.invalid=False
            else:
                print('The parameter %s is out of range' %(elem[3]), '\n')
                self.invalid=True
                    
        self.parameters = {}
        for i in range(len(self.para_name)):
            self.parameters[self.para_name[i]]=Fpar[i]
        if self.snow:
            self.parameters['CsnowM']=self.parameters['Csnow']/20
            self.parameters['Cswr'] = 0.00001
           
        # setting initial values
        self.TSK1[0]=self.parameters['CK1']
        self.TSK2[0]=self.parameters['CK2']
        self.U = self.relU*self.parameters['Umax']
        self.L = self.relL*self.parameters['Lmax']


    def write_outputs(self, OUTPUT_dir):
        i = str('{:02d}'.format(self.index_Fpar))
        print('Writing results to file {}.csv ...\n'.format(i))
        data_df = np.column_stack((self.ET,self.OF_r+self.IF_r,self.BF_r,self.S, self.IF_r, self.OF_r))
        df = pd.DataFrame(data = data_df, index = self.sim_time)
        df = df[self.EndWU:]
        df.columns= ['Evapotranspiration', 'Quick flow', 'Slow flow/Baseflow', 'Snow pack', 'Interflow', 'Overlandflow']
        df.index.name = 'Time'
        fname = os.path.join(OUTPUT_dir,'WB_{}.csv'.format(i))
        df.to_csv(fname)
        
        # Qt_df = pd.DataFrame(data = self.Qt, index=self.sim_time)
        # Qt_df.columns= ['Total flow']
        # Qt_df.index.name = 'Time'
        # fname1 = os.path.join(OUTPUT_dir,'TF_{}.csv'.format(i))
        # Qt_df[(self.Start+self.sim_time.freq):self.End].to_csv(fname1, sep=',')
        
    def plot_outputs(self, title):
        # one year at the time
        Sim_all = pd.DataFrame(data = self.Qt, index = self.sim_time)
        Sim_all.columns = ['SimFlow']
        Obs_all = pd.DataFrame(data = self.flow, index = self.sim_time)
        Obs_all.columns=['ObsFlow']
        #Start = np.argwhere(self.sim_time==np.datetime64('2005-03-02 00'))[0][0]
        #End = np.argwhere(self.sim_time==np.datetime64('2005-12-27 23'))[0][0]
        Sim = Sim_all[self.Start:self.End]
        Obs = Obs_all[self.Start:self.End]
        ax = Obs.plot(y='ObsFlow')
        Sim.plot(ax=ax, y = 'SimFlow', title=title)
        
        
    def _onestep(self, i):
        #print('Current time step: ',self.sim_time[i])
        if self.snow: 
            QS = 0.
            if self.temp[i]<=self.parameters['T0']: # snow
                if self.Sw==0.:
                    self.Ss+=self.precip[i]
                else:
                    QS = self.parameters['Csnow']*(self.parameters['T0']-self.temp[i])
                    if QS<=self.Sw:
                        self.Ss+=self.precip[i]+QS
                        self.Sw -= QS
                    else:
                        self.Ss+=self.precip[i]+self.Sw
                        self.Sw = 0.
            else: # snowmelt/rainfall
                if self.Ss==0.:
                    self.U+=self.precip[i]
                else:
                    if self.Ss>=0.9999:
                        QS = self.parameters['Csnow']*(self.temp[i]-self.parameters['T0'])
                    else:
                        QS = self.parameters['CsnowM']*(self.temp[i]-self.parameters['T0'])
                    if QS<=self.Ss:
                        self.Sw +=self.precip[i]+QS
                        self.Ss -=QS
                    else:
                        self.Sw+=self.precip[i]+self.Ss
                        self.Ss = 0.
                    if self.Sw>=(self.parameters['Cswr']*(self.Ss+self.Sw)):
                        self.U+=self.Sw
                        self.Sw=0.
            self.S[i] = self.Ss
            if self.Ss>0.999:
                self.pet[i] = 0.
                
        else:
            self.U+=self.precip[i] # no snow

        # Evapotranspiration
        if self.U>=self.pet[i]:
            evap = self.pet[i]
            self.U-=evap
        else:
            Ea = (self.pet[i]-self.U)*(self.L/self.parameters['Lmax'])
            evap = Ea + self.U
            self.L-=Ea
            self.U=0.
        self.ET[i]=evap
        
        # Interflow
        if (self.L/self.parameters['Lmax'])>self.parameters['TIF']:
            self.QIF[i]=(1/self.parameters['CKIF'])*(((self.L/self.parameters['Lmax'])-self.parameters['TIF'])/(1-self.parameters['TIF']))*self.U
            self.U-=self.QIF[i]
        else:
            self.QIF[i]=0.
        
        # Surface storage
        if self.U>self.parameters['Umax']:
            Pn = self.U-self.parameters['Umax']
            self.U = self.parameters['Umax']
            
            # overland flow
            if (self.L/self.parameters['Lmax'])>self.parameters['TOF']:
                self.QOF[i] = self.parameters['CQOF']*((self.L/self.parameters['Lmax'])-self.parameters['TOF'])/(1-self.parameters['TOF'])*Pn
            else:
                self.QOF[i] = 0.
                
        # groundwater recharge
            if (self.L/self.parameters['Lmax'])>self.parameters['TG']:
                G = (Pn-float(self.QOF[i]))*((self.L/self.parameters['Lmax'])-self.parameters['TG'])/(1-self.parameters['TG'])
            else:
                G = 0.
            self.G1[i]=G
            if self.G1[i]<0.:
                self.G1[i]=0.
        else:
            self.QOF[i] = 0.
            G = 0.
            self.G1[i]=0.
            Pn=0
        
        self.L+=(Pn-G-float(self.QOF[i]))
        
        # Flow routing
        # Linear routing baseflow
        bfp1 = (self.G1[i]*self.sub_area/self.fA)+((self.BF-(self.G1[i]*self.sub_area/self.fA)))*np.exp(-1/self.parameters['CKBF'])
        self.BF = bfp1
        self.BF_r[i] = self.BF
        
        # Two linear routing overland flow
        ifp1 = (self.QIF[i]*self.sub_area/self.fA)+(self.IFp-((self.QIF[i]*self.sub_area/self.fA)))*np.exp(-1/self.parameters['CK1'])
        self.IFp = ifp1
        ifp21 = self.IFp + (self.IF-self.IFp)*np.exp(-1/self.parameters['CK2'])
        self.IF = ifp21
        self.IF_r[i]=self.IF
        
        # Two linear routing overlabnd flow
        if (self.OFp*self.fA/self.sub_area)<self.OFmin:
            self.TSK1[i]=self.parameters['CK1']
        else:
            self.TSK1[i]=self.parameters['CK1']*(self.OF*self.fA/self.sub_area/self.OFmin)**(-self.Beta)
            
        if (self.OF*self.fA/self.sub_area)<self.OFmin:
            self.TSK2[i]=self.parameters['CK2']
        else:
            self.TSK2[i]=self.parameters['CK2']*(self.OF*self.fA/self.sub_area/self.OFmin)**(-self.Beta)
            
        if i == 0:
            self.OFp = ((1-np.exp(-1/self.TSK1[i]))*np.exp(-1/self.TSK1[i])*self.OFp)/(1-np.exp(-1/self.TSK1[i])) + (-1/self.TSK1[i])*self.QOF[i]*self.sub_area/self.fA
            self.OF = ((1-np.exp(-1/self.TSK2[i]))*np.exp(-1/self.TSK2[i])*self.OF)/(1-np.exp(-1/self.TSK2[i])) + (-1/self.TSK2[i])*self.OFp
            self.OF_r[i]=self.OF
        else:
            self.OFp = ((1-np.exp(-1/self.TSK1[i]))*np.exp(-1/self.TSK1[i-1])*self.OFp)/(1-np.exp(-1/self.TSK1[i-1])) + (-1/self.TSK1[i])*self.QOF[i]*self.sub_area/self.fA
            self.OF = ((1-np.exp(-1/self.TSK2[i]))*np.exp(-1/self.TSK2[i-1])*self.OF)/(1-np.exp(-1/self.TSK2[i-1])) + (-1/self.TSK2[i])*self.OFp
            self.OF_r[i]=self.OF

        # Total flow
        self.Qt[i] = self.OF+self.IF+self.BF
        
        
    def run_NAMModel(self, list_Fpar, OUTPUT_dir='.', isoutput=False, isplot=False):
        self.Total_flow = {}
        self.Evapotranspiration = {}
        self.Quick_flow = {}
        self.Base_flow = {}
        self.Snow_pack = {}
        self.Inter_flow = {}
        self.Overland_flow = {}
        
        #calculation of each parameter
        print('running NAM models.......')
        for j, Fpar in enumerate(list_Fpar):
            # Initialize the model by calling the constructor of the NAM class (in NAM.py)
            self.set_param(Fpar, j)
            # Running simulation
            if self.invalid:
                print('Simulation aborted because one of your parameters is out of range (see message above).')
            else:
                print('**********************************************************************\n' \
                      'Summary of your inputs:\nSTART DATE: \t{}END DATE: \t{}TIME STEP: \t{}' \
                      '\nCalibration parameters:\n'.format(
                          self.Start, self.End, self.time_step))
                print(*self.parameters)
                print(*self.parameters.values())
                print('**********************************************************************\n' \
                      'Running simulation...')
                #for date in self.sim_time:
                for i in range(len(self.Qt)):
                    self._onestep(i)
                self.Total_flow[j] = np.array(self.Qt.copy())
                self.Evapotranspiration[j] = np.array(self.ET.copy())
                self.Quick_flow[j] = np.array(self.OF_r + self.IF_r)
                self.Base_flow[j] = np.array(self.BF_r.copy())
                self.Snow_pack[j] = np.array(self.S.copy())
                self.Inter_flow[j] = np.array(self.IF_r.copy())
                self.Overland_flow[j] = np.array(self.OF_r.copy())
    
                #writing results
                if isoutput:
                    self.write_outputs(OUTPUT_dir)
            
                # check outflow
                if isplot:
                    self.plot_outputs(title=str(j))
    
        self.Total_flow = pd.DataFrame(self.Total_flow, index=self.sim_time)
        self.Evapotranspiration = pd.DataFrame(self.Evapotranspiration, index=self.sim_time)
        self.Quick_flow = pd.DataFrame(self.Quick_flow, index=self.sim_time)
        self.Base_flow = pd.DataFrame(self.Base_flow, index=self.sim_time)
        self.Snow_pack = pd.DataFrame(self.Snow_pack, index=self.sim_time)
        self.Inter_flow = pd.DataFrame(self.Inter_flow, index=self.sim_time)
        self.Overland_flow = pd.DataFrame(self.Overland_flow, index=self.sim_time)
        
        pd.DataFrame(list_Fpar).to_csv(OUTPUT_dir+'/Flow_Parameters.csv', header=self.para_name)
        if isoutput:
            fh_TF = os.path.join(OUTPUT_dir,'TF.csv')
            self.Total_flow.to_csv(fh_TF)

class AnalizeNAM():
    def __init__(self, list_Fpar, dataO, Total_flow, folderName='.'):
# read observed data
        self.dataO = dataO
        self.Total_flow = Total_flow
        self.aveO = np.nanmean(self.dataO)
# get file paths
        self.data = np.zeros((Total_flow.shape[1], 15))
        self.params = list_Fpar
        self.folderName = folderName
  
# calculation
    def calc(self):
        for i, dataS in self.Total_flow.items():
            self.dataS = dataS.values
            self.data[i] = np.array([self.NSE(), self.LNSE(), self.PBIAS(), *self.params[i]])
                        
            # print('NSE: {:5.3f}\nRMSE: {:5.3f}\nPBIAS: {:5.3f}'.format(*self.data[i]))

            #self.R_squared()
            #print(self.data[i])
        #print(self.data)
        self.write_errors()
        print('finished analyzing')
        
    def NSE(self):
        A = np.nansum(np.square(self.dataO-self.dataS))
        B = np.nansum(np.square(self.dataO-self.aveO))
        return 1-A/B

    def LNSE(self):
        etha = 10**(-4)
        A = np.nansum(np.square(np.log(self.dataO+etha)-np.log(self.dataS+etha)))
        B = np.nansum(np.square(np.log(self.dataO+etha)-np.log(self.aveO+etha)))
        return 1 - A/B    
    
    def PBIAS(self):
        n = np.count_nonzero(~np.isnan(self.dataO))
        A = np.nansum(self.dataO-self.dataS)
        return 100/n*A
    
    def MAE(self):
        n = np.count_nonzero(~np.isnan(self.dataO))
        A = np.nansum(np.abs(self.dataO-self.dataS))
        return A/n
    
    def RMSE(self):
        n = np.count_nonzero(~np.isnan(self.dataO))
        A = np.nansum(np.square(self.dataO-self.dataS))
        return (A/n)**.5

    def R_squared(self):
        res = stats.linregress(self.dataS, self.dataO)
        print(res)

# output
    def write_errors(self):
        header = ['NSE', 'LNSE', 'PBIAS', 'Umax', 'Lmax', 'CQOF', 'CKIF', 'CK1', 'CK2', 'TOF', 'TIF', 'TG', 'CKBF', 'Csnow', 'T0']
        self.data = pd.DataFrame(self.data, columns=header)
        self.data.to_csv("./"+self.folderName+"/analysis.csv")
        
