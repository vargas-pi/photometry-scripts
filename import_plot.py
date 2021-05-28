import tdt
import matplotlib.pyplot as py
from matplotlib.animation import FuncAnimation
import numpy as np
import json
import os
from datetime import datetime as dt
from threading import Thread, Event
import time
from collections import Callable

output_filename='exported_data'

class mouse_data:
     #class to store data we want to keep
    def __init__(self,mouse_id,F490,F405,fs,n=None,t=np.array([]),t_stim=None):
        self.mouse_id=mouse_id
        self.F490=F490
        self.F405=F405
        self.fs=fs
        self.n=n if n else len(F490)
        self.t=t if t.any() else np.arange(0,(self.n)/fs,1/fs)
        self.t_stim=t_stim if t_stim else None


class Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return {'ndarray': obj.tolist()}
        if isinstance(obj,Callable):
            return {'function':obj.__name__}
        if isinstance(obj,mouse_data):
            return {'mouse_data':obj.__dict__}
        return json.JSONEncoder.default(self, obj)

def hook(obj):
    if 'ndarray' in obj:
        return np.array(obj['ndarray'])
    if 'mouse_data' in obj:
        return mouse_data(**obj['mouse_data'])
    return obj
    
def stream_plot(ylim:list,mouse1=None,mouse2=None):
    try:
        syn=tdt.SynapseAPI('localhost')
    except ConnectionRefusedError:
        print('There was an error connecting to Synapse. Make sure a Synapse is recording and the Synapse API is turned on before running this script')
    mice=[(mouse1,'A'),(mouse2,'B')] #create a list of mouse ids with the channel character attached
    mice=[m for m in mice if m[0]] #only include mice that have been specified

    mice_data=[mouse_data(m[0],[],[],1) for m in mice] #an array of mouse_data objects for each mouse. they're initialized as empty
    for i in range(len(mice_data)):
        mice_data[i].lim490=ylim.copy()
        mice_data[i].lim405=ylim.copy()
        mice_data[i].limt=[0,60]
        mice_data[i].updated=False
        
    def stream(start,green_light):
        while green_light.is_set():
            
            for i in range(len(mice_data)):
                #NOTE: At some point this should be updated so the gizmo/parameter name is a variable. 
                # Right now this function would probably only work on our setup
                new_490=syn.getParameterValue('FibPho1','Response-1'+mice[i][1])
                new_405=syn.getParameterValue('FibPho1','Response-2'+mice[i][1])
                new_t=(dt.now()-start).total_seconds()
                mice_data[i].F490.append(new_490)
                mice_data[i].F405.append(new_405)
                mice_data[i].t=np.append(mice_data[i].t,new_t)

                if new_t>=mice_data[i].limt[1]:
                    mice_data[i].limt[1]+=60
                    mice_data[i].updated=True
                    
                if new_490>=mice_data[i].lim490[1]:
                    mice_data[i].lim490[1]+=1
                    mice_data[i].updated=True
                elif new_490<=mice_data[i].lim490[0]:
                    mice_data[i].lim490[0]-=1
                    mice_data[i].updated=True
                    
                if new_405*8>=mice_data[i].lim405[1]:
                    mice_data[i].lim405[1]+=1
                    mice_data[i].updated=True
                elif new_405*8<=mice_data[i].lim405[0]:
                    mice_data[i].lim405[0]-=1
                    mice_data[i].updated=True
                    
            time.sleep(.1)

    green_light=Event()
    green_light.set()
    start=dt.now()
    stream_thread=Thread(target=stream,args=(start,green_light))
    stream_thread.start()

    try:
        fig , lines, ax= initialize_plots(mice_data,ylim)
        def animate(i):
            for m in range(len(mice_data)):
                if mice_data[m].updated:
                    mice_data[m].updated=False
                    if len(mice_data)==2:
                        ax[0,m].set_ylim(mice_data[m].lim490)
                        ax[1,m].set_ylim(mice_data[m].lim405)
                        ax[0,m].set_xlim(mice_data[m].limt)
                        ax[1,m].set_xlim(mice_data[m].limt)
                    else:
                        ax[0].set_ylim(mice_data[m].lim490)
                        ax[1].set_ylim(mice_data[m].lim405)
                        ax[0].set_xlim(mice_data[m].limt)
                        ax[1].set_xlim(mice_data[m].limt)
                        
                    fig.canvas.draw()
            update_plots(mice_data,lines)
            return lines
        anim=FuncAnimation(fig,animate,interval=1,blit=True)
        py.show()
        while True:
            time.sleep(.1)
    except KeyboardInterrupt:
        print('stopping the stream...')
        green_light.clear()
        stream_thread.join()
        print('stream successfully stopped')

        

def open_save(path,format='npy',mouse1=None,mouse2=None,m1_490='x19a',m1_405='x15a',m2_490='x29B',m2_405='x25B'):
    raw_data=tdt.read_block(path)
    print('data pulled!')

    data=[]
    if mouse1:
        m1=mouse_data(mouse1, raw_data.streams[m1_490].data, raw_data.streams[m1_405].data, raw_data.streams[m1_490].fs) #create an instance of the data class
        data.append(m1)
    if mouse2:
        m2=mouse_data(mouse2, raw_data.streams[m2_490].data, raw_data.streams[m2_405].data, raw_data.streams[m2_490].fs)
        data.append(m2)

    print('saving data...')

    if format=='json':
        with open(os.path.join(path,output_filename+'.json'), 'w') as f:
            json.dump(data,f,cls=Encoder)
    elif format=='npy':
        np.save(os.path.join(path,output_filename),data)
    else:
        raise Exception('Unrecognized File Format!')

    print('save successful!')
    
    return data

def update_plots(d,lines):
    u00=min(len(d[0].t),len(d[0].F490))
    u10=min(len(d[0].t),len(d[0].F405))
    if len(d)==2:
        u01=min(len(d[1].t),len(d[1].F490))
        u11=min(len(d[1].t),len(d[1].F405))
        
        lines[0].set_data(d[0].t[:u00],d[0].F490[:u00])
        lines[1].set_data(d[1].t[:u01],d[1].F490[:u01])
        lines[2].set_data(d[0].t[:u10],8*np.array(d[0].F405[:u10]))
        lines[3].set_data(d[1].t[:u11],8*np.array(d[1].F405[:u11]))
    elif len(d)==1:
        lines[0].set_data(d[0].t[:u00],d[0].F490[:u00])
        lines[1].set_data(d[0].t[:u10],8*np.array(d[0].F405[:u10]))

def initialize_plots(d,ylim:list):
    fig , ax = py.subplots( *([2]*len(d)) )

    if len(d)==2:
        ax[0,0].set_title(d[0].mouse_id)
        l0,=ax[0,0].plot([],[],'g', linewidth=0.5)
        ax[0,0].set_ylim(*ylim)
        ax[0,0].set_xlim(0,60)

        ax[0,1].set_title(d[1].mouse_id)
        l1,=ax[0,1].plot([],[],'g', linewidth=0.5)
        ax[0,1].set_ylim(*ylim)
        ax[0,1].set_xlim(0,60)
        
        l2,=ax[1,0].plot([],[],'r', linewidth=0.5)
        ax[1,0].set_ylim(*ylim)
        ax[1,0].set_xlim(0,60)
        
        l3,=ax[1,1].plot([],[],'r', linewidth=0.5)
        ax[1,1].set_ylim(*ylim)
        ax[1,1].set_xlim(0,60)
        
        lines=[l0,l1,l2,l3]
    elif len(d)==1:
        ax[0].set_title(d[0].mouse_id)
        l0,=ax[0].plot([],[],'g', linewidth=0.5)
        ax[0].set_ylim(*ylim)
        ax[0].set_xlim(0,60)
        
        l1,=ax[1].plot([],[],'r', linewidth=0.5)
        ax[1].set_ylim(*ylim)
        ax[1].set_xlim(0,60)
        
        lines=[l0,l1]
        
    return [fig,lines,ax]

if __name__=='__main__':
    
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-path',help='the path to the output directory if you are viewing data retroactively')
    parser.add_argument('-mouse1',help='the identifier of the mouse on the first line')
    parser.add_argument('-mouse2',help='the identifier of the mouse on the second line')
    parser.add_argument('-m1_490',help='the name of the 490 channel for line 1', default='x19A')
    parser.add_argument('-m1_405',help='the name of the 405 channel for line 1', default='x15A')
    parser.add_argument('-m2_490',help='the name of the 490 channel for line 2', default='x29B')
    parser.add_argument('-m2_405',help='the name of the 405 channel for line 2', default='x25B')
    parser.add_argument('-ylim',help='initial y axis limits for streaming and plotting', nargs=2, default=[0.,3.], type=float)
    parser.add_argument('-format',help='file format to save the data to (either npy or json)', default='npy')
    args=parser.parse_args()

    if not (args.mouse1 or args.mouse2):
        parser.error('must enter the identifier for at least one mouse')
    
    if args.path:
        d=open_save(args.path, 
                    format=args.format,
                    mouse1=args.mouse1, 
                    mouse2=args.mouse2, 
                    m1_490=args.m1_490, 
                    m1_405=args.m1_405, 
                    m2_490=args.m2_490, 
                    m2_405=args.m2_405)
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
    else:
        stream_plot(args.ylim,mouse1=args.mouse1,mouse2=args.mouse2)
        

 
