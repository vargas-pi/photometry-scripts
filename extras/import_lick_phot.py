import tdt
import datetime
import numpy as np
import pandas as pd
import sys
sys.path.append('../')
from analysis_pipeline import *
from sys_parameters import DEFAULT_SAVE_CHANNELS 


def open_save_lick(fpath, mouse1 = None, mouse2 = None, channels = DEFAULT_SAVE_CHANNELS ):
    
    """
    a modified version of the open_save function in the import_plot.py file.
    This version is meant specifically for experiments using our custom built lickometer system
    along with fiber photometry. The goal is to automate the task of creating mouse_data objects
    with the appropriate stim times given lickometer data

    Parameters
    -----------
    fpath: str
        path to the folder of lickometer data
    channels: dict
        a dictionary of lists of photometry channel labels. the keys should be wavelengths (490 and 405)
        and the lists should be ordered meaningfully and consistently for a given system. For example,

    """

    fpath = Path(fpath).resolve()
    l_files = list(fpath.iterdir())
    csvs = list(filter(lambda x: x != "link.txt", l_files)) # get the paths to all csv files
    licks = {} # initialize a dict to store lick data
    for i in csvs: # loop over csv files
        df = pd.read_csv(i, index_col = 0) # read the file
        df['datetime'] = pd.to_datetime(df.datetime) # cnvert the datetime column to datetimes
        if df.loc[0,['left lick count',  'right lick count']].sum()>0: # check if this file was initialized with a blank data point
            #if not calculate the start time and append a blank value to the beginning of the dataframe
            new_st = df.datetime[0] - datetime.timedelta(milliseconds = df['time offset (ms)'].astype('float64')[0])
            df = pd.DataFrame( dict(zip(df.columns, [new_st, 0, 0, 0])), index = [0]).append(df, ignore_index = True)
        
        licks.update({i.name.split('_')[1]:df})

    with open(fpath/"link.txt", 'w') as f: link = Path(f.readline()).resolve() # get the address to the photometry data block
    block = list( list( link.iterdir() )[0].iterdir() )[0] # we need to go 2 levels deep for the tank path
    phot_data = tdt.read_block(block) # read in the photometry data
    subj = phot_data.info.subject.split("_") # pull the subjects
    tstart = phot_data.info.start_date # also the datetime of the beginning
    data = []  # initialize list of mouse data objects
    
    if mouse1:
        m = mouse_data( mouse1, phot_data.streams[channels[1][490]].data, 
                        phot_data.streams[channels[1][405]].data,
                        phot_data.streams[channels[1][490]].fs,
                        t_start = tstart,
                        t_stim = (licks[mouse1][0]['datetime'] - tstart).total_seconds())
        m.left_licks = licks[mouse1].loc[licks[s]['left lick count'],'time offset (ms)'].values
        m.right_licks = licks[mouse1].loc[licks[s]['right lick count'],'time offset (ms)'].values
        data.append(m)
    
    if mouse2:
        m = mouse_data( mouse2, phot_data.streams[channels[2][490]].data, 
                        phot_data.streams[channels[2][405]].data,
                        phot_data.streams[channels[2][490]].fs,
                        t_start = tstart,
                        t_stim = (licks[mouse2][0]['datetime'] - tstart).total_seconds())
        m.left_licks = licks[mouse2].loc[licks[s]['left lick count'],'time offset (ms)'].values
        m.right_licks = licks[mouse2].loc[licks[s]['right lick count'],'time offset (ms)'].values
        data.append(m)

    np.save(block/'exported_data.npy', data)


if __name__ == '__main__':
    pass
