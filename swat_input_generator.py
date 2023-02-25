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
        with open(os.path.join(self.path, ext+'.cli'), 'w') as f:
            #create one <ext>.cli file            
            f.write('{}: {} file names - file written by {} {}\n'.format(
                ext+'.cli', self.file_types[ext], self.editor, self.ct))
            f.write('filename\n')
            
            for ID, data in self.Datasets.items():
                if ext in self.Datasets[ID]['daily_datasets'].keys():
                    station_name = data['name']
                    datafilename = station_name+'.'+ext
                    f.write('{}\n'.format(datafilename))
                    
                    #create <station name>.<ext> files
                    self._weather_data_files_out(ID, ext, datafilename)

        print('saved {}'.format(ext+'.cli'))


    def _weather_data_files_out(self, ID, ext, datafilename):
        datasets = self.Datasets[ID]
        df_daily = datasets['daily_datasets']
        
        with open(os.path.join(self.path, datafilename), 'w') as f:
            f.write('{}: {} data - file written by {} {}\n'.format(
                datafilename, self.file_types[ext], self.editor, self.ct))
                    
            number_of_year = df_daily.shape[0]//365+1
            time_step = 0
            
            f.write('nbyr     tstep       lat       lon      elev\n')
            f.write('{:>4d} {:>9d} {:>9.3f} {:>9.3f} {:>9.3f}\n' .format(
                number_of_year, time_step, datasets['lat'], datasets['lon'], datasets['elev']))
            
            for index, row in df_daily.iterrows():
                #       year, day of year, data
                if ext=='tmp':
                    f.write('{:>4d} {:>4d} {:>10.5f} {:>10.5f}\n'.format(
                        int(row['year']), int(row['day']), row['tmp_max'], row['tmp_min']))
                else:                        
                    f.write('{:>4d} {:>4d} {:>10.5f}  \n'.format(
                        int(row['year']), int(row['day']), row[ext]))
        
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
    'hmd' : 'Relative humidity',
    'pcp' : 'Precipitation',
    'tmp' : 'Temperature',
    'wnd' : 'Wind speed',
    'slr' : 'Solar radiation',
}

if __name__ == '__main__':
    CSV_folder = r'./csv/'
    Station_Datasets = read_station_data(os.path.join(CSV_folder, 'WGEN_stat.csv'))
    
    for ID, station in Station_Datasets.items():
        datafile = station['name']+'.csv'
        daily_datasets = pd.read_csv(os.path.join(CSV_folder, datafile))
        
        print('Reading {}...'.format(datafile))
        Station_Datasets[ID]['daily_datasets'] = daily_datasets
    
    SWAT = SWATplusInputs(File_TYPES, Station_Datasets)
    SWAT.path = './weather_swat_plus'
    os.makedirs(SWAT.path, exist_ok=True)

    for ext in File_TYPES.keys():
        SWAT.weather_inputs_out(ext)
