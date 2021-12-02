#TODO: i should write something to automate adding data from a bunch of files
#TODO: need to make sure t_stim can be saved in import plot


import matplotlib.pyplot as py
from utilities import *
import os
from scipy import stats as st
import pandas as pd
from typing import Dict
from pathlib import Path
from copy import deepcopy
from scipy.io import savemat
from pandas import IndexSlice as idx



class analysis:
    """
    Class for performing post-hoc analyses of fiber photometry data
    ...

    Attributes
    ----------
    norm_method : func
        a callable function for cropping and normalizing all data in 
        this analysis.this function must take the following as arguments:                 

        rec: mouse_data
            an instance of a mouse_data object storing the raw data for a
            given mouse
        t_endrec: float
            the amount of time in seconds from stim to the end of the recording 
            to keep in the analysis
            
        in turn, it should return the following:

        normed_490: numpy.ndarray
            1D array of normalized 490 data 
        normed_405: numpy.ndarray
            1D array of normalized 405 data 
        start_ind: int
            index of the beginning of the cropped region of the recording relative
            to the raw data
        stim_ind: int
            index of the stim time relative to the raw data
        end_ind: int
            index of the end of the cropped region of the recording relative
            to the raw data

    t_endrec : float
        the amount of time in seconds, from the stim to the end of the recording, 
        to keep in the analysis

    ex : float, optional
        the number of standard deviations above and below the mean of a given mouse's
        data beyond which to exclude (default=4)
        
    file_format: str, optional
        the file type to save the analysis to. This can be either 'npy' or 'json' (default npy). 
        NOTE: json is mostly still an option due to older versions of the code that used json 
        for cross platform compatibility. It is strongly recommended to use npy files as they are
        much faster to read and write to. If you would like to export normalized data from this analysis
        to further analyze with another language, the export_to_mat method may be useful
        
    file_loc: str, optional
            the path to a directory to save the analysis to

    Methods
    -------
    load_append_save(file:Path,t_stims:Dict[str,float])
        load, normalize and add data to the analysis. when finished, save the analysis

    normalize_downsample(rec:mouse_data,plot=True)
        call the speicified normalization method and downsample to ~1 Hz. when finished,
        plot the normalized data.

    compute_means()
        compute the average and standard error accross recordings

    recompute()
        re-crop/normalize all data and recompute the mean signal after updating parameters

    save(file_format=False)
        export the analysis object to a file in the output folder if the output folder 
        hasn't been created, this will create it

    export_to_mat()
        forward fill all nan values in the normalized 490 data and export the resulting
        array to a .mat file

    remove_mouse(mouse:str)
        remove a mouse from the analysis

    retrieve_excluded(mouse:str)
        retrieve an excluded mouse from the analysis

    plot_both()
        plot the average normalized 490 and 405 signal with error

    plot_ind_trace(mouse:str)
        plot the individual trace for a given mouse

    plot_490()
        plot the average normalized 490 signal with error

    bin_data(binsize:int,save:bool)
        bin the data average the data in each bin for each mouse

    bin_plot(binsize:int,save:bool)
        bin the data and plot it

    bin_avg(start:int,end:int)
        average over a specified section of data for each mouse

    bin_auc(start:int,end:int)
        area under the curve for a specified section of data for each mouse

    ind_peak_df_f(extrema:str)
        determine either the min or max ∆f/f and time for each mouse separately

    mean_peak_df_f(extrema:str)
        determine either the min or max ∆f/f in the mean signal and identify the values of 
        the normed 490 at this time point in each individual mouse

    time_to_half_pk(extrema:str)
        ...



    """

    def __init__(self, norm_method,t_endrec,t_prestim=300,ex=4,file_format='npy',file_loc='analyses'):
        """
        Parameters
        ----------
        norm_method : func
            see Attributes
        t_endrec : float
            see Attributes
        ex : float, optional
            see Attributes
        file_format: str, optional
            see Attributes
        file_loc: str, optional
            see Attributes

        Returns
        --------
        self: analysis

        """

        self.conditions=None
        self.norm_method=norm_method
        self.t_endrec=t_endrec
        self.raw_data= []
        self.excluded_raw = []
        self.loaded=False
        self.ex=ex
        self.file_loc=file_loc
        self.file_format=file_format
        self.t_prestim=t_prestim

        # NOTE: the following properties are no longer used but keeping here
        # in case of compatibility issues for earlier versions of the code
        self.stim_ind=None 
        self.normed_data = []
        self.excluded_normed = []
    """
    define a few getter and setter functions for relevant parameters in the 
    analysis such that when they are updated, everything is recomputed automatically
    """

    @property
    def t_endrec(self): return self._t_endrec
    @t_endrec.setter
    def t_endrec(self,value):
        self._t_endrec=value
        self.compute()
    
    @property
    def t_prestim(self): return self._t_prestim
    @t_prestim.setter
    def t_prestim(self,value):
        self._t_prestim=value
        self.compute()

    @property
    def ex(self): return self._ex
    @ex.setter
    def ex(self,value):
        self._ex=value
        self.compute()
    
    @property
    def norm_method(self): return self._norm_method
    @norm_method.setter
    def norm_method(self,value):
        self._norm_method=value
        self.compute()

    @property
    def file_format(self): return self._file_format
    @file_format.setter
    def file_format(self,value):
        if value not in ['json','npy']:
            raise Exception('Unrecognized File Format')
        else:
            self._file_format=value

    
    def load_append_save(self,file:Path,t_stims:Dict[str,float]=None,conds:Dict[str,float]=None):
        """
        load, normalize and add data to the analysis. when finished, save the analysis
        
        Parameters
        ----------
        file: pathlib.Path
            the path to the file to load into the analysis
        t_stims: Dict[str,float]
            a dictionary matching each mouse in the specified file to its respective stimulus time
        """

        if file.suffix=='.json':
            with open(file,'r') as f:
                print('loading data from file...')
                d=json.load(f,object_hook=hook)
        elif file.suffix=='.npy':
            print('loading data from file...')
            d=np.load(file,allow_pickle=True).tolist()
        else:
            raise Exception('Unrecognized File Format!')

        for i in range(len(d)): 
            if t_stims is not None:
                #should i add an option to deal with multiple stim times?
                d[i].t_stim=t_stims[d[i].mouse_id]
            elif not hasattr(d[i],'t_stim'):
                raise NoStimTime
            if conds is not None:
                #should i add an option to deal with multiple stim times?
                d[i].cond=conds[d[i].mouse_id]
            elif not hasattr(d[i],'cond'):
                d[i].cond=0
            #increment the trial number for the given mouse
            d[i].trial=len(list(filter( lambda x: x==d[i].mouse_id, self.raw_data)))
        
        self.raw_data.extend(d)
        self.recompute()
        self.loaded=True
        self.save()

    def normalize_downsample(self,rec:mouse_data,plot=True):
        """
        call the speicified normalization method and downsample to ~1 Hz. when finished, plot the normalized data.
        NOTE: the sampling rate from synapse isn't a whole number so we can't necessarilly get exaclt 1Hz but it's close

        Parameters
        ----------
        rec: mouse_data
            the mouse_data object with the raw data to be normalized and downsampled
        plot: bool,optional
            whether or not to plot the normalized data at the end

        Returns
        -------
        m: mouse_data
            the normalized downsampled data for the given mouse
        """

        normed_490, normed_405, t = self.norm_method(rec,self.t_endrec,self.t_prestim)
        #resample the data to 1Hz via linear interpolation
        ds_t,ds_490=resample(t,normed_490)
        _,ds_405=resample(t,normed_405)

        self.t=ds_t.copy()

        #exclude outliers
        excl=lambda x: (x>np.mean(x)+self.ex*np.std(x))|(x<np.mean(x)-self.ex*np.std(x))
        mask=excl(ds_405)
        
        ds_t[mask]=np.nan
        ds_490[mask]=np.nan
        ds_405[mask]=np.nan

        m=deepcopy(rec)
        m.F490, m.F405, m.t, m.t_stim= ds_490, ds_405, ds_t, 0

        if hasattr(rec,'t_start'):
            #if we saved the timestamps of this recording, downsample the time stamps as well
            t_start=rec.t_start+rec.t_stim*timedelta(seconds=1)
            #create a new mouse_data object for the cleaned data
            m.t_start=t_start

        #plot the data
        if plot:
            py.plot(m.t,100*m.F490,'g',linewidth=0.5)
            py.plot(m.t,100*m.F405,'r',linewidth=0.5)
            py.axvline(x=0, c='k',ls='--', alpha=0.5)
            py.xlabel('Time Relative to Stimulus (s)')
            py.ylabel(r'$\frac{\Delta F}{F}$ (%)')
            py.show()

        return m

    def get_mouse_raw(self,mouse):
        return list(filter( lambda x: x.mouse_id==mouse, self.raw_data))
    
    def compute(self):
        """
        crop/normalize all included raw data and compute the mean signal
        """
        if hasattr(self,'excluded_normed') and hasattr(self,'normed_data'): #make sure the fields for normed data exist
            if (len(self.raw_data)+len(self.excluded_raw))>0: #make sure there is data loaded in the analysis (excluded or not)
                print('recomputing...')

                #clear the dataframes
                cols=pd.MultiIndex(levels=[[]]*3,codes=[[]]*3,names=('cond','mouse','trial'))
                self.all_490=pd.DataFrame([],columns=cols)
                self.all_405=pd.DataFrame([],columns=cols)

                for i,r in enumerate(self.raw_data): #loop through the raw data and redo the normalization/downsampling

                    if not hasattr(r,'cond'): self.raw_data[i].cond,r.cond=0,0
                    if not hasattr(r,'trial'):
                        get_trial=lambda x: x.trial if hasattr(x,'trial') else -1
                        mi=list(map( get_trial, self.get_mouse_raw(r.mouse_id) ))
                        self.raw_data[i].trial,r.trial=[max(mi)+1]*2

                    m=self.normalize_downsample(r,plot=False)
                    self.all_490[m.cond,m.mouse_id,m.trial]=pd.Series(m.F490,index=m.t)
                    self.all_405[m.cond,m.mouse_id,m.trial]=pd.Series(m.F405,index=m.t)

                self.conditions=self.all_490.columns.get_level_values('cond')

                #consolidate dataframe by averaging within mice within condition
                self.all_490=self.all_490.groupby('cond',axis=1).apply(lambda x: x.groupby('mouse',axis=1).mean())
                self.all_405=self.all_405.groupby('cond',axis=1).apply(lambda x: x.groupby('mouse',axis=1).mean())

                #calculate mean and standard error for 490
                self.mean_490=self.all_490.groupby('cond',axis=1).mean()
                self.err_490=self.all_490.groupby('cond',axis=1).sem()

                #calculate mean and standard error for 405
                self.mean_405=self.all_405.groupby('cond',axis=1).mean()
                self.err_405=self.all_405.groupby('cond',axis=1).sem()

                print('successful')

    def save(self,file_format=False):
        """
        export the analysis object to a file in the output folder
        if the output folder hasn't been created, this will create it

        Parameters
        ----------
        file_format: str, optional
            see Attributes

        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
            
        if not os.path.exists(self.file_loc): os.mkdir(self.file_loc)
        ftype= self.file_format if not file_format else file_format

        if ftype not in ['npy','json']:
            raise Exception('Unrecognized file format!')

        if ftype=='json':
            with open(os.path.join(self.file_loc,'analysis_'+'_'.join([r.mouse_id for r in self.raw_data])+'.json'),'w') as f:
                print('saving the analysis...')
                json.dump(self.__dict__,f,cls=Encoder)
                print('save successful!')
        elif ftype=='npy':
            np.save(os.path.join(self.file_loc,'analysis_'+'_'.join([r.mouse_id for r in self.raw_data])+'.npy'),self)
            print('save successful!')
    
    def export_to_mat(self):
        """
        forward fill all nan values in the normalized 490 data and export
        the resulting array to a .mat file
        """

        d={}
        for c in self.all_490.columns.get_level_values('cond'):
            data=self.all_490[c].values.T
            if np.isnan(data).any():
                #check if there are nans and forward fill them first
                r,c=np.where(np.isnan(data))
                for i,j in zip(r,c):
                    data[i,j]=data[i,j-1] if j>0 else 0
            d.update({f'_{c}':data})
        
        savemat(os.path.join(self.file_loc,'analysis_'+'_'.join([r.mouse_id for r in self.raw_data])+'.mat'),d)


    def remove_mouse(self,mouse:str):
        """
        remove a mouse from the analysis

        Parameters
        -----------
        mouse: str
            the id of the mouse to remove
        """

        excl_raw=list(filter(lambda x: x.mouse_id==mouse, self.raw_data))

        if not hasattr(self,'excluded_raw'):
            self.excluded_raw=[]

        self.excluded_raw.extend(excl_raw)
        self.raw_data=list(filter(lambda x: x.mouse_id!=mouse, self.raw_data))

        if len(self.raw_data)==0:
            self.loaded=False
            
        self.compute()

    def retrieve_excluded(self,mouse:str):
        """
        retrieve an excluded mouse from the analysis

        Parameters
        -----------
        mouse: str
            the id of the excluded mouse to retrieve
        """
        if not hasattr(self,'excluded_raw') or len(self.excluded_raw)==0:
            print('this analysis has no excluded data')
            return

        ret_raw=list(filter(lambda x: x.mouse_id==mouse, self.excluded_raw))
        self.raw_data.extend(ret_raw)
        self.excluded_raw=list(filter(lambda x: x.mouse_id!=mouse, self.excluded_raw))
        self.loaded=True

        self.compute()


    def plot_both(self,show=True,ax=None):
        """
        plot the average normalized 490 and 405 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        if ax is None: _,ax=py.subplots(1,1)

        for i in self.mean_490:
            ax.fill_between(self.t, 100*self.mean_490[i] + 100*self.err_490[i],
                            100*(self.mean_490[i] - self.err_490[i]), color='g',alpha=0.2 )
            ax.fill_between(self.t, 100*(self.mean_405[i] + self.err_405[i]),
                            100*(self.mean_405[i] - self.err_405[i]), color='r',alpha=0.2  )
            ax.plot(self.t, 100*self.mean_490[i] , 'g', linewidth=0.5)
            ax.plot(self.t, 100*self.mean_405[i], 'r', linewidth=0.5)
        ax.axvline(x=0, c='k',ls='--', alpha=0.5)
        ax.set_xlabel('Time Relative to Stimulus (s)')
        ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')

        if show:
            py.show()

    def plot_ind_trace(self,mouse:str,cond=None,plot_405=True,ax=None,show=True):
        """
        plot the individual trace for a given mouse

        Parameters
        ----------
        mouse:mouse_data
            mouse to plot
        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        if ax is None: _,ax=py.subplots(1,1)
        if cond is not None:
            ax.plot(100*self.all_490.loc[:,cond,mouse] , 'g', linewidth=0.5)
            if plot_405:
                ax.plot(100*self.all_405.loc[:,cond,mouse], 'r', linewidth=0.5)
        else:
            #TODO: make the 490 traces different shades of green for each condition and add a legend
            ax.plot(100*self.all_490.loc[:,idx[:,mouse]] , 'g', linewidth=0.5)
            if plot_405:
                ax.plot(100*self.all_405.loc[:,idx[:,mouse]], 'r', linewidth=0.5)
        ax.axvline(x=0, c='k',ls='--', alpha=0.5)
        ax.set_xlabel('Time Relative to Stimulus (s)')
        ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
        
        if show:
            py.show()

        
    def plot_490(self,show=True,ax=None):
        """
        plot the average normalized 490 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        if ax is None: _,ax=py.subplots(1,1)
        for i in self.mean_490:
            ax.fill_between(self.t, 100*self.mean_490[i] + 100*self.err_490[i],
                            100*(self.mean_490[i] - self.err_490[i]), color='g',alpha=0.2 )
            ax.plot(self.t, 100*self.mean_490[i] , 'g', linewidth=0.5)
        ax.axvline(x=0, c='k',ls='--', alpha=0.5)
        ax.set_xlabel('Time Relative to Stimulus (s)')
        ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')

        if show:
            py.show()


    def bin_plot(self,binsize,save=False,show=True,ax=None):
        """
        run bin_plot and then plot the data

        Parameters
        ----------
        binsize: int
            bin size in seconds
        save: bool
            indicates whether or not to save the dataframe to a csv file for further analysis
        Returns
        -------
        binned: pandas.DataFrame
            binned data
        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        df=self.bin_data(binsize,save=save)
        
        if ax is None: _,ax=py.subplots(1,1)
        for i in df.columns.get_level_values('cond'):
            ax.errorbar(x=df.index,y=100*df['mean'][i],yerr=100*df['sem'][i])
        ax.set_xlabel('Time Relative to Stimulus (s)')
        ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
        
        if show:
            py.show()


    def bin_data(self,binsize,save=False):
        """
        bin the data average the data in each bin for each mouse

        Parameters
        ----------
        binsize: int
            bin size in seconds

        Returns
        -------
        binned: pandas.DataFrame
            binned data
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        #start the bins in both directions from 0
        bins=np.append(-np.arange(0, self.t_prestim+binsize, binsize)[:0:-1],
                        np.arange(0, self.t[-1]+binsize, binsize ))

        mice=[m.mouse_id for m in self.raw_data] # get the names of the kept mice in the analysis

        df=self.all_490.copy()
        for i in df.columns.get_level_values('cond'):
            df[i,'bin start (s)']=bins[np.digitize(self.t,bins)-1]

        def mean_f_trapz(y):
            #calculate the mean within each bin by numerical integration
            y=y.dropna()
            return np.trapz(y,x=y.index)/(y.index[-1]-y.index[0])

        def do_binning(x):
            x=x.droplevel('cond',axis=1)
            return x.pivot_table(values=x.columns[x.columns.isin(mice)], 
                                 index=['bin start (s)'],aggfunc=mean_f_trapz)

        binned=df.groupby('cond',axis=1).apply(do_binning)
        binned_stats=pd.concat([binned.groupby('cond',axis=1).mean(),binned.groupby('cond',axis=1).sem()],keys=['mean','sem'],axis=1)


        if not hasattr(self,'binned'):
            self.binned={}
        self.binned.update({f'{binsize}':binned_stats})

        if save:
            binned_stats.to_csv(os.path.join(self.file_loc,f"binned_{binsize}s {'_'.join([r.mouse_id for r in self.raw_data])}.csv"))
        
            
        return binned_stats


    def bin_avg(self,start:int,end:int,save=False):
        """
        average over a specified section of data for each mouse

        Parameters
        ----------
        start: int
            start time of the bin
        end: int
            end time of the bin

        Returns
        -------
        avg_dict: Dict[str,float]
            dictionary pairing each mouse with the avg normed 490 value
            in the specified bin
            
        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        aucs=self.bin_auc(start,end,save=False,pr=False)
        if (end>self.t_endrec) or (start<-self.t_prestim):
            print('Warning, bin is out of bounds')

        if aucs is not None:
            avgs=aucs/(end-start)

            print(' ')
            print(avgs)

            if not hasattr(self,'avgs'):
                self.avgs={}
            self.avgs.update({f'_{start}_{end}':avgs})
            
            if save:
                avgs.to_csv(os.path.join(self.file_loc,f"avg_{start}s_{end}s_{'_'.join([r.mouse_id for r in self.raw_data])}.csv"))

            return avgs
    
    def bin_auc(self,start:int,end:int,save=False,pr=True):
        """
        area under the curve of a specified section of data for each mouse

        Parameters
        ----------
        start: int
            start time of the bin
        end: int
            end time of the bin

        Returns
        -------
        avg_dict: Dict[str,float]
            dictionary pairing each mouse with the avg normed 490 value
            in the specified bin
            
        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        nearest=lambda arr,val: arr[np.abs(arr-val).argmin()]
        start_t=nearest(self.t,start)
        end_t=nearest(self.t,end)

        if (end>self.t_endrec) or (start<-self.t_prestim):
            print('Warning, bin is out of bounds')

        try: 
            y=self.all_490.loc[start_t:end_t].copy()
            aucs=y.apply(lambda x: np.trapz(x=x.dropna().index,y=x.dropna().values),axis=0)
            aucs=pd.DataFrame(aucs,columns=['Mean ∆F/F'])

            if pr:
                print(' ')
                print(aucs)

            if not hasattr(self,'aucs'):
                self.aucs={}
            self.aucs.update({f'_{start}_{end}':aucs})

            if save:
                aucs.to_csv(os.path.join(self.file_loc,f"auc_{start}s_{end}s_{'_'.join([r.mouse_id for r in self.raw_data])}.csv"))

        except IndexError:
            print('Time stamps are outside the scope of the analysis! Choose a different start and end time or restart the analysis with a different t_endrec.')
            return None
        
        return aucs

    def ind_peak_df_f(self,extrema:str,save=False):
        """
        determine either the min or max ∆f/f and time for each mouse separately
        TODO: we should also be saving the times of the peaks

        Parameters
        ----------
        extrema: str
            extrema (either 'min' or 'max') of the peak to search for
        
        Returns
        -------
        pks_dict: Dict[str,float]
            dictionary pairing each mouse with the peak value
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return


        if extrema=='max':
            pk_inds=self.all_490.loc[0:].apply(lambda x: x.argmax(skipna=True))
        elif extrema=='min':
            pk_inds=self.all_490.loc[0:].apply(lambda x: x.argmin(skipna=True))
        else:
            print('Unrecognized extrema!')
        if ((self.t_endrec-pk_inds)<=20).any(): print(f'Warning! {extrema} is on the edge of the recording')
        pk_inds[self.t_endrec-pk_inds<=20]-=20

        peaks={}
        for i,v in pk_inds.items():
            y=self.all_490[0:].iloc[v-5:v+6][i].dropna()
            peaks.update({i:np.trapz(x=y.index,y=y.values)})

        peaks=pd.DataFrame(pd.Series(peaks), columns=[f'{extrema} ∆F/F'])
        pk_times=pd.DataFrame(pk_inds.apply(lambda x: self.all_490[0:].index[x]), columns=['time (s)'])

        peaks=pd.concat([peaks,pk_times],axis=1)

        print(peaks)
        print('')
        if not hasattr(self,'pks'):
            self.pks={}

        self.pks.update({f'{extrema}':peaks})

        if save:
            df=pd.DataFrame(peaks,columns=[f'{extrema} ∆F/F','time (s)']).T
            df.to_csv(os.path.join(self.file_loc,f"ind_{extrema}_{'_'.join([r.mouse_id for r in self.raw_data])}.csv"))

        return peaks


    def mean_peak_df_f(self,extrema:str,save=False):
        """
        determine either the min or max ∆f/f in the mean signal and identify the values of the normed 490
        at this time point in each individual mouse

        Parameters
        ----------
        extrema: str
            extrema (either 'min' or 'max') of the peak to search for
        
        Returns
        -------
        pks_dict: Dict[str,float]
            dictionary pairing each mouse with the peak value
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        pk_inds=self.mean_490.loc[0:].apply(lambda x: x.argmax(skipna=True))
        if ((self.t_endrec-pk_inds)<=20).any(): print(f'Warning! {extrema} is on the edge of the recording')
        pk_inds[(self.t_endrec-pk_inds)<20]-=20

        peaks={}
        for i,v in pk_inds.items():
            y=self.all_490[0:].iloc[v-10:v+10][i].dropna()
            y=y.apply(lambda x: np.trapz(x=x.dropna().index,y=x.dropna().values))
            y.index=y.index.map(lambda x: (i,x))
            peaks.update(y.to_dict())

        peaks=pd.DataFrame(pd.Series(peaks), columns=[f'{extrema} ∆F/F'])
        ts=pk_inds.apply(lambda x: self.all_490[0:].index[x])
        peaks['ts']=np.nan
        for i in ts.index:
            peaks.loc[i,:]['ts']=ts[i]

        print('')
        print(peaks)
        if not hasattr(self,'mean_pks'):
            self.mean_pks={}

        self.mean_pks.update({f'{extrema}':peaks})

        if save:
            peaks.to_csv(os.path.join(self.file_loc,f"mean_{extrema}_{'_'.join([r.mouse_id for r in self.raw_data])}.csv"))

        return peaks


    def time_to_half_pk(self):
        raise NotImplementedError

    #old versions of functions
    def compute_means(self):
        """
        compute the average and standard error accross recordings
        (Note: this function has been consolidated into compute)
        """
        raise NotImplementedError

    
    def recompute(self):
        """
        re-crop/normalize all data and recompute the mean signal after updating parameters
        (Note: this function has been consolidated into compute)
        """
        raise NotImplementedError


class loaded_analysis(analysis):
    """a class for analyses loaded from json files"""
    def __init__(self,fpath):
        """
        Parameters
        ----------
        fpath: str
            path to the anaalysis file
        """
        with open(fpath,'r') as f: d=json.load(f,object_hook=hook)
        if d['norm_method']=='function':
            print('WARNING: No normalization  method is specified in this loaded analysis! Please update analysis parameters before adding any data.')
        else:
            d['norm_method']=getattr(analysis,d['norm_method']['function'])
            
        for key in d:
            if isinstance(d[key],list):
                if len(d[key])>0:
                    setattr(self,key,d[key])
            else:
                setattr(self,key,d[key])

def load_analysis(fpath):
    """
    load an exported analysis (whether json or npy)

    Parameters
    ----------
    fpath: str
        path to analysis file

    Returns
    -------
    a: analysis
        the analysis object
    """
    fpath=Path(fpath).resolve()
    if fpath.suffix=='.json':
        a=loaded_analysis(fpath)
    elif fpath.suffix=='.npy':
        a=np.load(fpath,allow_pickle=True).item()
    else:
        raise Exception('Unrecognized file format!')
    
    a.file_loc=fpath.parent
    a.compute()
    return a