import tdt
import datetime
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as py
sys.path.append('../')
from analysis_pipeline import *
from sys_parameters import *


def open_save_lick(fpath, mouse1 = None, mouse2 = None, cond1 = 0, cond2 = 0, channels = DEFAULT_SAVE_CHANNELS ):
    
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
    csvs = list(filter(lambda x: x.name != "link.txt", l_files)) # get the paths to all csv files
    licks = {} # initialize a dict to store lick data
    for i in csvs: # loop over csv files
        df = pd.read_csv(i, index_col = 0) # read the file
        df['datetime'] = pd.to_datetime(df.datetime) # cnvert the datetime column to datetimes
        if df.loc[0,['left lick count',  'right lick count']].sum()>0: # check if this file was initialized with a blank data point
            #if not calculate the start time and append a blank value to the beginning of the dataframe
            new_st = df.datetime[0] - datetime.timedelta(milliseconds = df['time offset (ms)'].astype('float64')[0])
            df = pd.DataFrame( dict(zip(df.columns, [new_st, 0, 0, 0])), index = [0]).append(df, ignore_index = True)
        
        licks.update({i.name.split('_')[1]:df})

    with open(fpath/"link.txt", 'r') as f: link = Path(f.readline()).resolve() # get the address to the photometry data block
    block = list( list( link.iterdir() )[0].iterdir() )[0] # we need to go 2 levels deep for the tank path
    phot_data = tdt.read_block(block) # read in the photometry data
    tstart = phot_data.info.start_date # also the datetime of the beginning
    data = []  # initialize list of mouse data objects
    
    if mouse1:
        m = mouse_data( mouse1, 
                        phot_data.streams[channels[1][490]].data, 
                        phot_data.streams[channels[1][405]].data,
                        phot_data.streams[channels[1][490]].fs,
                        t_start = tstart,
                        t_stim = (licks[mouse1].loc[0,'datetime'] - tstart).total_seconds(),
                        cond = cond1)
        m.left_licks = licks[mouse1].loc[licks[mouse1]['left lick count'],'time offset (ms)'].values
        m.right_licks = licks[mouse1].loc[licks[mouse1]['right lick count'],'time offset (ms)'].values
        data.append(m)
    
    if mouse2:
        m = mouse_data( mouse2, 
                        phot_data.streams[channels[2][490]].data, 
                        phot_data.streams[channels[2][405]].data,
                        phot_data.streams[channels[2][490]].fs,
                        t_start = tstart,
                        t_stim = (licks[mouse2].loc[0,'datetime'] - tstart).total_seconds(),
                        cond = cond2)
        m.left_licks = licks[mouse2].loc[licks[mouse2]['left lick count'],'time offset (ms)'].values
        m.right_licks = licks[mouse2].loc[licks[mouse2]['right lick count'],'time offset (ms)'].values
        data.append(m)

    np.save(block/'exported_data.npy', data)
    return data


if __name__ == '__main__':

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-path',help='the path to the output directory if you are viewing data retroactively')
    parser.add_argument('-mouse1',help='the identifier of the mouse on the first line')
    parser.add_argument('-mouse2',help='the identifier of the mouse on the second line')
    parser.add_argument('-m1_490',help='the name of the 490 channel for line 1', default=None)
    parser.add_argument('-m1_405',help='the name of the 405 channel for line 1', default=None)
    parser.add_argument('-m2_490',help='the name of the 490 channel for line 2', default=None)
    parser.add_argument('-m2_405',help='the name of the 405 channel for line 2', default=None)
    parser.add_argument('-cond1', help='condition for mouse 1')
    parser.add_argument('-cond2', help='condition for mouse 2')
    args=parser.parse_args()

    if not (args.mouse1 or args.mouse2):
        parser.error('must enter the identifier for at least one mouse')
    
    if args.path:
        # update save channels dict
        channels = DEFAULT_SAVE_CHANNELS
        channels[1][490] = args.m1_490 if args.m1_490 else channels[1][490]
        channels[1][405] = args.m1_405 if args.m1_405 else channels[1][405]
        channels[2][490] = args.m2_490 if args.m2_490 else channels[2][490]
        channels[2][405] = args.m2_405 if args.m2_405 else channels[1][405]

        # read block and save
        d=open_save_lick(args.path, 
                        mouse1=args.mouse1, 
                        mouse2=args.mouse2,
                        cond1=args.cond1,
                        cond2=args.cond2,
                        channels = channels)
        #plot
        fig , ax = py.subplots( *([2]*len(d)) )
        if len(d)==2:
            ax[0,0].set_title(d[0].mouse_id)
            ax[0,0].plot(d[0].t,d[0].F490,'g', linewidth=0.5)
            ax[0,1].set_title(d[1].mouse_id)
            ax[0,1].plot(d[1].t,d[1].F490,'g', linewidth=0.5)
            ax[1,0].plot(d[0].t,8*np.array(d[0].F405),'r', linewidth=0.5)
            ax[1,1].plot(d[1].t,8*np.array(d[1].F405),'r', linewidth=0.5)
        elif len(d)==1:
            ax[0].set_title(d[0].mouse_id)
            ax[0].plot(d[0].t,d[0].F490,'g', linewidth=0.5)
            ax[1].plot(d[0].t,8*d[0].F405,'r', linewidth=0.5)
        py.show()
