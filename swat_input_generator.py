import os
import datetime
import pandas as pd

'''
Station_Datasets:
    {
        'name'      : The weather station name
        'lat'       : Latitude of weather station used to create statistical parameters (degrees).
        'lon'       : Longitude of weather station (degrees)
        'elve'      : Elevation of weather station (m)
        'rain_yrs'  : The number of years of maximum monthly 0.5 h rainfall data 
                        If no value is input for RAIN_YRS, SWAT will set RAIN_YRS = 10
        #*above data are obtained from WGEN_stat.csv (read_station_data())
        
        'daily_datasets': pandas DataFrame:
            year, day,  hmd,  pcp, ...,
            1990,   1, 0.30, 0.50, ...,
             :      :    :    :    ...,
    }
'''

class SWATplusInputs:
    '''
    Generate inputs files for SWAT+ editor.
    Parameters
    -------------------------------------------------
        file_types: dict
            {extention1: description1, ...},
        Station_Datasets: dict
            see above
    '''
    def __init__(self, file_types, Station_Datasets):
        self.path = '.'
        self.editor = 'SWAT+ editor'        #file generation editor
        self.ct = datetime.datetime.now()   #file generation time

        self.file_types = file_types
        self.Datasets = Station_Datasets

        
    def weather_inputs_out(self, ext):
        '''
        main function to generate inputs inside folder 'weather_swat_plus' 
        create one <ext>.cli file and <station name>.<ext> files for each station

        Parameters
        ----------
        ext : str
            extention of selected weather input file type (see File_TYPES).

        '''
        os.makedirs(self.path, exist_ok=True)

        with open(os.path.join(self.path, ext+'.cli'), 'w') as f:
            #create one <ext>.cli file            
            f.write('{}: {} file names - file written by {} {}\n'.format(
                ext+'.cli', self.file_types[ext], self.editor, self.ct))
            f.write('filename\n')
            
            for ID, data in self.Datasets.items():
                try:
                    station_name = data['name']
                    datafilename = station_name+'.'+ext   
                    #create <station name>.<ext> files
                    self._weather_data_files_out(ID, ext, datafilename)

                    f.write('{}\n'.format(datafilename))

                except:
                    pass

        print('saved {}'.format(ext+'.cli'))


    def _weather_data_files_out(self, ID, ext, datafilename):
        datasets = self.Datasets[ID]
        df_daily = datasets[ext]
        
        with open(os.path.join(self.path, datafilename), 'w') as f:
            f.write('{}: {} data - file written by {} {}\n'.format(
                datafilename, self.file_types[ext], self.editor, self.ct))
                    
            number_of_year = df_daily.index[-1].year - df_daily.index[0].year + 1
            time_step = 0
            
            f.write('nbyr     tstep       lat       lon      elev\n')
            f.write('{:>4d} {:>9d} {:>9.3f} {:>9.3f} {:>9.3f}\n' .format(
                number_of_year, time_step, datasets['lat'], datasets['lon'], datasets['elev']))
            
                #       year, day of year, data
            if ext=='tmp':
                for index, row in df_daily.iterrows():
                    f.write('{:>4d} {:>4d} {:>10.5f} {:>10.5f}\n'.format(
                        index.year, index.timetuple().tm_yday, row['tmp_max'], row['tmp_min']))
            else:
                for index, row in df_daily.iteritems():
                    f.write('{:>4d} {:>4d} {:>10.5f}  \n'.format(
                        index.year, index.timetuple().tm_yday, row))
        
        print('saved '+datafilename)


    def Update_Station_Data(self, kwarg):
        self.__dict__.update(**kwarg)
        

def read_station_data(path):
    '''
    Read station info from WGEN_stat.csv

    Parameters
    ----------
    path : str
        relative path to WGEN_stat.csv

    Returns
    -------
    Station_Datasets : dict
        see above.
    '''
    Station_Datasets = {}
    DataFrame = pd.read_csv(path, header=0)
    for index, row in DataFrame.iterrows():
        station_data = {
            'name': row['name'],
            'lat' : row['lat'],
            'lon' : row['lon'],
            'elev': row['elev'],
            'rain_yrs': row['rain_yrs'],
            }
        Station_Datasets[index] = station_data
    
    return Station_Datasets

# %%
File_TYPES = {
#    'hmd' : 'Relative humidity',
    'pcp' : 'Precipitation',
    'tmp' : 'Temperature',
#    'wnd' : 'Wind speed',
#    'slr' : 'Solar radiation',
}

# %%
if __name__ == '__main__':
    CSV_folder = r'./csv/'
    Station_Datasets = read_station_data(os.path.join(CSV_folder, 'WGEN_stat.csv'))
    nan = -99
    
    for ext in File_TYPES.keys():
        if ext == 'tmp':
            datafile_max = 'tmp_max.csv'
            datafile_min = 'tmp_min.csv'
            print('Reading {}, {}...'.format(datafile_max, datafile_min))
            daily_datasets_min = pd.read_csv(os.path.join(CSV_folder, datafile_min))
            daily_datasets_max = pd.read_csv(os.path.join(CSV_folder, datafile_max))
            
            for key, value in Station_Datasets.items():
                try:
                    daily_datasets = pd.concat([daily_datasets_max[value['name']], daily_datasets_min[value['name']]], axis=1)
                    daily_datasets.index = pd.to_datetime(daily_datasets_min['date'], format='%d/%m/%Y')
                    daily_datasets.columns = ['tmp_max', 'tmp_min']
                    daily_datasets = daily_datasets.fillna(nan)

                    Station_Datasets[key][ext] = daily_datasets
                except:
                    print('missing data at ', value['name'])
        else:
            datafile = ext+'.csv'
            print('Reading {}...'.format(datafile))
            daily_datasets = pd.read_csv(os.path.join(CSV_folder, datafile))
            daily_datasets.index = pd.to_datetime(daily_datasets['date'], format='%d/%m/%Y')
            daily_datasets = daily_datasets.fillna(nan)
    
            for key, value in Station_Datasets.items():
                try:
                    Station_Datasets[key][ext] = daily_datasets[value['name']]
                except:
                    print('missing data at ', value['name'])

        SWAT = SWATplusInputs(File_TYPES, Station_Datasets)
        SWAT.path = './weather_swat_plus'
        SWAT.weather_inputs_out(ext)
