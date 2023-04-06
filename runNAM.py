"""
Editted on Tue Dec 14 2021 @author: toru
"""

import NAM
import pandas as pd
import os

# %% Setup values ##################################################################################

Menu = {
    1: 'One-time',
    2: 'Calibration'
    }
selectMode = 2                                                            #select 1, 2

##Mode=='Calibration'
caliblen = 10                                                            #separation number

Parameters = {
    # Min = [5, 50, 0.0, 150, 3, 3, 0, 0, 0, 500, 2, -3]
    # Max = [35, 1000, 1.0, 2000, 48, 48, 0.9, 0.9, 0.99, 7500, 5, 3]
    
    # name, Min, Max, Default
    'Umax': {'Min':   5,  'Max':   35,  'Default':   15},
    'Lmax': {'Min':  50,  'Max': 1000,  'Default':  500}, 
    'CQOF': {'Min':  .2,  'Max':  1.0,  'Default':   .5},
    'CKIF': {'Min': 500,  'Max': 2000,  'Default': 1000},
    'CK1':  {'Min':  15,  'Max':   48,  'Default':   20},
    'CK2':  {'Min':   3,  'Max':   48,  'Default':   20},
    'TOF':  {'Min':   0,  'Max':   .9,  'Default':   .5},
    'TIF':  {'Min':   0,  'Max':   .9,  'Default':   .5},
    'TG':   {'Min':   0,  'Max':  .99,  'Default':   .5},
    'CKBF': {'Min': 500,  'Max': 7500,  'Default': 3000},
    'Csnow':{'Min':   2,  'Max':    5,  'Default':    3},
    'T0' :  {'Min':   0,  'Max':    3,  'Default':    2},
}

# Define simulation warm-up period (YYYY-MM-DD HH)

StartDate = '2000-01-01'                                                # Defining starting date of the simulation
EndDate = '2006-12-31'                                                  # Defining ending date of the simulation
EndWUDate = '2001-12-31'                                                # Defining ending date of warming-up period
    
time_step = 'D'                                                         # D: daily, H: hourly

OUTPUT_dir = 'test'                                                     # output folder
DATA_dir=r'Data'                                                        # input folder

sub_area = 'SubCatchmArea.csv'                                          # input data files
flow = 'Flow.xlsx'
pet = 'PET.xlsx'
precip = 'Precip.xlsx'

snow = False                                                            # will you use the snow module? (you need temperature data!) 
Temp = 'Temp.xlsx'

relU = 0.08                                                             # Initial values for states
relL = 0.05
OF = 0.0
OFp = 0.0
IF = 0.0
IFp = 0.0
BF = 10
Beta = 0.40
OFmin = 0.3631   
Ss = 0                                                                  # extra states
Sw = 0

max_iter = 20
DEBUG = False

#####################################################################################
# %% functions

def List_FlowPar(Fpar, Mode):            
    if Mode == 'Calibration':
        calibnum = [val['Default'] for val in Parameters.values()]
        calibMin = [val['Min'] for val in Parameters.values()]
        calibMax = [val['Max'] for val in Parameters.values()]
        list_Fpar=[]
        len_Fpar = caliblen*(len(calibnum))
        for j, name in enumerate(Parameters.keys()):
            for i in range(caliblen):
                Min = Parameters[name]['Min']
                Max = Parameters[name]['Max']
                Fpar2 = Fpar.copy()
                Fpar2[j] = ((caliblen-i-1)*Min+i*Max)/(caliblen-1)
                list_Fpar.append(Fpar2)
    
    else:             #one-time calculation
        list_Fpar=[Fpar]
    
    return list_Fpar

def Main():
    # %% Setup
    Mode = Menu[selectMode]

    fh_flow = os.path.join(DATA_dir,flow)
    fh_pet = os.path.join(DATA_dir,pet)
    fh_precip = os.path.join(DATA_dir,precip)
    fh_sub_area = os.path.join(DATA_dir,sub_area)

    # Obtain time series
    flow_df = pd.read_excel(fh_flow, sheet_name ='Flow', usecols=[0,1,2,3,4])
    pet_df = pd.read_excel(fh_pet, sheet_name ='PET', usecols=[0,1,2,3,4])
    precip_df = pd.read_excel(fh_precip, sheet_name='Precip', usecols=[0,1,2,3,4])
    if snow:
        temp_df = pd.read_excel(fh_temp, sheet_name='Temp', usecols = [0,1,2,3,4])
    
    with open(fh_sub_area) as f:
        lines = f.readlines()
    area = float(lines[0])
    
    # setup NAM
    MyModel=NAM.NAM(StartDate, EndDate, EndWUDate, time_step, area)
    # update time series
    MyModel.input_timeseries(flow_df, pet_df, precip_df)
    MyModel.snow = snow
    
    # %% ========================= optional ===================================
    MyModel.OFmin = OFmin
    MyModel.Beta = Beta
    MyModel.IF = IF
    MyModel.IFp = IFp
    MyModel.BF = BF
    MyModel.Ss = Ss
    MyModel.Sw = Sw
    MyModel.OFp = OFp
    MyModel.OF = OF
    MyModel.relU = relU
    MyModel.relL = relL
    MyModel.DEBUG = DEBUG
    #======================================================================
    
    # %% Running model
    Optimal = [val['Default'] for val in Parameters.values()]
    i = 0

    if Mode == 'Calibration':
        while True:
            # to prevent large numbers of iteration
            if i > max_iter:
                break
            
            # create matrix of parameters
            list_Fpar = List_FlowPar(Optimal, Mode=Mode)
            
            # run NAM model
            MyModel.run_NAMModel(list_Fpar, OUTPUT_dir)
            
            # check errors
            efficiency = NAM.AnalizeNAM(list_Fpar, MyModel.flow_df, MyModel.Total_flow, OUTPUT_dir)
            efficiency.calc()
            
            # find the best set of parameters
            data = efficiency.data.sort_values(by=['NSE'], ascending=False)
            New_Optimal = data.iloc[0].values[3:].tolist()
            
            # break if the best parameters was the same as previous ones
            if New_Optimal == Optimal:
                print('Iteration: {}\nOptimal values\n{}'.format(i+1, Optimal))
                break
            # update parameters
            Optimal = New_Optimal
            
            i += 1
        
        #========================= optional ===================================
        # Run the optimal model with plot and csv output
        list_Fpar = List_FlowPar(Optimal, Mode='One-time')
        #running with output
        MyModel.run_NAMModel(list_Fpar, OUTPUT_dir, isplot=True, isoutput=True)

    else:
        isplot=True
        isoutput=True
        
        list_Fpar = List_FlowPar(Optimal, Mode)
        
        #running NAM with one time
        MyModel.run_NAMModel(list_Fpar, OUTPUT_dir, isplot=isplot, isoutput=isoutput)
        #analyze output data
        efficiency = NAM.AnalizeNAM(list_Fpar, MyModel.flow_df, MyModel.Total_flow, OUTPUT_dir)
        efficiency.calc()
        
if __name__ == '__main__':
    Main()
