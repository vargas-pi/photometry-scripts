#TODO: i should write something to automate adding data from a bunch of files
#TODO: raise actual errors instead of printing issues, from the cli we can use try, except clauses to catch the errors


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
import seaborn as sns



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
    multi_index = pd.MultiIndex(levels = [[]]*3, codes = [[]]*3, names = ('cond','mouse','trial'))

    def __init__(self, norm_method, t_endrec, ds_freq = 1, t_prestim=300, ex=5, 
                 detrend = False, detrend_method = detrend_405_constrained, 
                 file_format='npy', file_loc='analyses',fname=None):
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

        self.loaded = False
        self.conditions = None
        self.ds_freq = ds_freq
        self.norm_method = norm_method
        self.t_endrec = t_endrec
        self.detrend = detrend
        self.detrend_method = detrend_method
        self.normed_data = pd.Series([], index = analysis.multi_index, dtype = object)
        self.raw_data = pd.Series([], index = analysis.multi_index, dtype = object)
        self.excluded_raw = pd.Series([], index = analysis.multi_index, dtype = object)
        self.ex = ex
        self.file_loc = file_loc
        self.file_format = file_format
        self.t_prestim = t_prestim
        self.fname = fname

        # NOTE: the following properties are no longer used but keeping here
        # in case of compatibility issues for earlier versions of the code
        self.stim_ind = None 
        self.excluded_normed = []
    """
    define a few getter and setter functions for relevant parameters in the 
    analysis such that when they are updated, everything is recomputed automatically
    """
    @property
    def t(self):
        return np.arange( - self.t_prestim, self.t_endrec + 1/self.ds_freq, step = 1/self.ds_freq)

    @property
    def file_format(self): return self._file_format
    @file_format.setter
    def file_format(self,value):
        if value not in ['json','npy']:
            raise Exception('Unrecognized File Format')
        else:
            self._file_format=value

    
    def load_append_save(self, file:Path,t_stims:Dict[str,float]=None, conds:Dict[str,float]=None):
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
            if len(self.raw_data)>0:
                try: d[i].trial=len(self.raw_data.sort_index().loc[d[i].cond,d[i].mouse_id, :])
                except KeyError: d[i].trial=0
            else: d[i].trial=0
            self.raw_data[d[i].cond, d[i].mouse_id, d[i].trial] = d[i]
        
        self.loaded=True

    def normalize_downsample(self, rec:mouse_data,plot=True, scale=1):
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
        m = deepcopy(rec)
        #resample the data to 1Hz via linear interpolation
        m.center_stim()
        m.resample(self.t)
        if self.detrend: 
            m.F490, m.F405 = self.detrend_method(m)
        m.F490, m.F405 = self.norm_method(m)

        #plot the data
        if plot:
            py.plot(m.t,scale*m.F490,'g',linewidth=0.5)
            py.plot(m.t,scale*m.F405,'r',linewidth=0.5)
            py.axvline(x=0, c='k',ls='--', alpha=0.5)
            py.xlabel('Time Relative to Stimulus (s)')
            if scale==100:
                py.ylabel(r'$\frac{\Delta F}{F}$ (%)')
            else:
                py.ylabel(r'$\frac{\Delta F}{F}$')
            py.show()

        return m
    
    def rename_cond(self, old_name, new_name):
        """
        rename a condition in the 
        """
        def rename(x):
            x.cond = new_name if x.cond==old_name else x.cond 
        self.raw_data.apply(rename) # rename the condition within the mouse_data object
        self.excluded_raw.apply(rename) # same for the excluded
        tmp = self.raw_data.reset_index()
        tmp['cond']= tmp[0].apply(lambda x: x.cond) # set the cond column equal to whatever the cond field is in the updated mouse data object
        self.raw_data = tmp.set_index(['cond','mouse','trial'])[0] # convert back to a series
        tmp = self.excluded_raw.reset_index()
        tmp['cond']= tmp[0].apply(lambda x: x.cond) # same for excluded
        self.excluded_raw = tmp.set_index(['cond','mouse','trial'])[0]
        self.compute()

    def compute(self, log = True):
        """
        crop/normalize all included raw data and compute the mean signal
        """
        if self.loaded: #make sure there is data loaded in the analysis
            if log: print('recomputing...')

            # clear the fields for normed data
            self.all_490=pd.DataFrame([],columns = analysis.multi_index )
            self.all_405=pd.DataFrame([],columns = analysis.multi_index )
            self.normed_data = pd.Series([], index = analysis.multi_index, dtype = object)            

            # loop through the raw data and redo the normalization/downsampling
            for index, r in self.raw_data.items():
                m = self.normalize_downsample(r, plot=False)
                self.normed_data[index] = m
                self.all_490[index] = m.F490
                self.all_405[index] = m.F405

            self.all_490.index = self.t
            self.all_405.index = self.t
            self.all_490.fillna(method='ffill')
            
            self.conditions=self.all_490.columns.get_level_values('cond')

            # consolidate dataframe by averaging within mice within condition
            self.all_490=self.all_490.groupby('cond',axis=1).apply(lambda x: x.groupby('mouse',axis=1).mean())
            self.all_405=self.all_405.groupby('cond',axis=1).apply(lambda x: x.groupby('mouse',axis=1).mean())

            # calculate mean and standard error for 490
            self.mean_490=self.all_490.groupby('cond',axis=1).mean()
            self.err_490=self.all_490.groupby('cond',axis=1).sem()

            # calculate mean and standard error for 405
            self.mean_405=self.all_405.groupby('cond',axis=1).mean()
            self.err_405=self.all_405.groupby('cond',axis=1).sem()

            # update properties of the analysis
            self.mice_by_cond=self.all_490.groupby('cond',axis=1).apply(lambda x:x.columns.get_level_values('mouse').tolist()).to_dict()
            self.n_by_cond=self.all_490.groupby('cond',axis=1).apply(lambda x:x.columns.get_level_values('mouse').size).to_dict()
            self.mice=self.all_490.columns.get_level_values('mouse').unique()
            self.n_tot=self.mice.size
            self.conds=self.all_490.columns.get_level_values('cond').unique()

            if log: print('successful')

    def save(self, file_format=False):
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
        
        if hasattr(self,'fname'):
            fname= 'analysis_'+'_'.join(list(self.mice)) if self.fname==None else self.fname
        else:
            fname= 'analysis_'+'_'.join(list(self.mice))

        if ftype not in ['npy','json']:
            raise Exception('Unrecognized file format!')

        if ftype=='json':
            with open(os.path.join(self.file_loc,fname+'.json'),'w') as f:
                print('saving the analysis...')
                json.dump(self.__dict__,f,cls=Encoder)
                print('save successful!')
        elif ftype=='npy':
            np.save(os.path.join(self.file_loc,fname+'.npy'),self)
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
        
        savemat(os.path.join(self.file_loc,'analysis_'+'_'.join(self.mice)+'.mat'),d)


    def remove_mouse(self, mouse:str, cond=None, recompute=True):
        """
        remove a mouse from the analysis

        Parameters
        -----------
        mouse: str
            the id of the mouse to remove
        """

        if not hasattr(self,'excluded_raw'):
            self.excluded_raw = pd.Series([], index = analysis.multi_index, dtype = object)
        if cond is None:
            self.excluded_raw.append(self.raw_data.xs(pd.IndexSlice[:,mouse], drop_level=False))
            self.raw_data=self.raw_data.drop(index=[mouse],level=1)
        else:
            self.excluded_raw.append(self.raw_data.xs(pd.IndexSlice[cond,mouse], drop_level=False))
            self.raw_data=self.raw_data.drop(index=(cond,mouse))
        if len(self.raw_data)==0:
            self.loaded=False        
        if recompute:
            self.compute()

    def retrieve_excluded(self,mouse:str, cond=None, recompute=True):
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

        if cond is None:
            #self.raw_data.append(self.excluded_raw.xs(pd.IndexSlice[:,mouse], drop_level=False))
            #cond is 0 for now
            self.raw_data[:, mouse, 0]=self.excluded_raw[:, mouse, 0]
            self.excluded_raw = self.excluded_raw.drop(index=[mouse],level=1)
        else:
            #self.raw_data.append(self.excluded_raw.xs(pd.IndexSlice[cond,mouse], drop_level=False))
            tr=self.excluded_raw[cond, mouse].index[0]
            self.raw_data[cond, mouse, tr]=self.excluded_raw[cond, mouse, tr]
            self.excluded_raw = self.excluded_raw.drop(index=(cond,mouse))

        self.loaded=True

        if recompute:
            self.compute()


    def plot_both(self, cond=None, c490='g', c405='r', ylim = None,
                  show=True, ax=None, alpha=.3, figsize=(12,5), scale = 1):
        """
        plot the average normalized 490 and 405 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        cond = self.conds[0] if self.conds.size==1 else cond
        if cond is not None:
            if ax is None: _,ax=py.subplots(1,1)
            ax.fill_between(self.t, scale*(self.mean_405[cond] + self.err_405[cond]),
                            scale*(self.mean_405[cond] - self.err_405[cond]), color=c405,alpha=alpha  )
            ax.plot(self.t, scale*self.mean_405[cond], color=c405, linewidth=.5)
            ax.fill_between(self.t, scale*self.mean_490[cond] + scale*self.err_490[cond],
                            scale*(self.mean_490[cond] - self.err_490[cond]), color=c490,alpha=alpha )
            ax.plot(self.t, scale*self.mean_490[cond] , color=c490, linewidth=.5)
            if scale==100:
                ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
            else:
                ax.set_ylabel(r'$\frac{\Delta F}{F}$')
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.set_xlabel('Time Relative to Stimulus (s)')
        
        else:
            if ax is None:
                _,ax=py.subplots(1,self.conds.size,figsize=figsize)
            elif isinstance(ax, np.ndarray):
                if ax.size<self.conds.size:
                    raise Exception('provided axes have invalid dimensions. creating a new one...')
            bnds=[]
            if self.conds.size>=1:
                for j,i in enumerate(self.conds):
                    ax.flatten()[j].fill_between(self.t, scale * (self.mean_405[i] + self.err_405[i]),
                                                 scale * (self.mean_405[i] - self.err_405[i]), 
                                                 color = c405, alpha = alpha  )
                    ax.flatten()[j].plot(self.t, scale * self.mean_405[i], color = c405, linewidth=.5)
                    bnds.append(ax.flatten()[j].get_ylim())
                
                for j,i in enumerate(self.conds):
                    ax.flatten()[j].fill_between(self.t, scale * self.mean_490[i] + scale*self.err_490[i],
                                                 scale * (self.mean_490[i] - self.err_490[i]), 
                                                 color = c490, alpha = alpha)
                    ax.flatten()[j].plot(self.t, scale * self.mean_490[i] , color=c490, linewidth=.5)
                    bnds.append(ax.flatten()[j].get_ylim())
                    
                    ax.flatten()[j].set_title(i)
                    if scale==100:
                        ax.flatten()[j].set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
                    else:
                        ax.flatten()[j].set_ylabel(r'$\frac{\Delta F}{F}$')
                    ax.flatten()[j].axvline(x=0, c='k',ls='--', alpha=0.5)
                    if ylim is not None:
                        ax.flatten()[j].set_ylim(ylim)
                ax.flatten()[-1].set_xlabel('Time Relative to Stimulus (s)')
                mins,maxes=np.array(bnds).T
                if ylim is None:
                    for i in range(self.conds.size): 
                        ax.flatten()[i].set_ylim(min(mins),max(maxes))
        if show:
            py.show()

        return ax

    def plot_ind_trace(self,mouse:str, cond=None, cm405='Reds', cm490='Greens', c490='g',
                       c405='r', plot_405=True, ax=None, show=True, scale=1, linewidth = 0.5,  ylim = None):
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
        cm490=sns.color_palette(cm490,self.conds.size)
        cm405=sns.color_palette(cm405,self.conds.size)
        
        if cond is not None:
            ax.plot(scale * self.all_490.loc[:,idx[cond,mouse]] , 
                    color = c490, linewidth = linewidth, label= cond)
            if plot_405:
                ax.plot(scale * self.all_405.loc[:,idx[cond,mouse]], 
                        color = c405, linewidth = linewidth)
            ax.legend() #added
        else:
            d = scale * self.all_490.loc[:, idx[:, mouse]]
            d.columns = d.columns.get_level_values('cond')
            d5 = scale * self.all_405.loc[:, idx[:,mouse]]
            d5.columns = d.columns.get_level_values('cond')

            if plot_405:
                ax.set_prop_cycle(color = cm405)
                d5.plot.line(linewidth = linewidth, ax = ax)
            
            ax.set_prop_cycle(color = cm490)
            d.plot.line(linewidth = linewidth, ax = ax)
            if self.conds.size > 1: ax.legend()
       
        ax.axvline(x=0, c='k',ls='--', alpha=0.5)
        ax.set_xlabel('Time Relative to Stimulus (s)')
        if scale==100:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
        else:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$')

        if ylim is not None:
            ax.set_ylim(ylim)
        
        if show:
            py.show()

        return ax

        
    def plot_490(self, cond=None, show=True, ax=None, cm='Set1', c490='g', alpha=.3, scale=1, ylim = None):
        """
        plot the average normalized 490 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        if ax is None: _,ax=py.subplots(1,1)
        cm=dict(zip(self.conds,sns.color_palette(cm,self.conds.size)))
        
        if cond is not None:
            ax.fill_between(self.t, scale*self.mean_490[cond] + scale*self.err_490[cond],
                            scale*(self.mean_490[cond] - self.err_490[cond]), color=c490,alpha=alpha )
            ax.plot(self.t, scale*self.mean_490[cond] , color=c490, linewidth=0.5, label=cond)
            ax.legend()

        else:
            for i in self.conds:
                ax.fill_between(self.t, scale*self.mean_490[i] + scale*self.err_490[i],
                                scale*(self.mean_490[i] - self.err_490[i]), color=cm[i],alpha=alpha )
                ax.plot(self.t, scale*self.mean_490[i] , color=cm[i], linewidth=0.5,label=i)
                if self.conds.size>1: ax.legend()
        ax.axvline(x=0, c='k',ls='--', alpha=0.5)
        ax.set_xlabel('Time Relative to Stimulus (s)')
        
        if scale==100:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
        else:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$')

        if ylim is not None:
            ax.set_ylim(ylim)

        if show:
            py.show()
        
        return ax


    def bin_plot(self, binsize, save=False, cond=None, show=True, ax=None, cm='Set2', color=None, scale=1, ylim = None):
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

        df,_=self.bin_data(binsize,save=save)

        if ax is None: _,ax=py.subplots(1,1)
        cm=dict(zip(self.conds,sns.color_palette(cm,self.conds.size)))

        if cond is not None:
            color=cm[cond] if color is None else color
            ax.errorbar(x=df.index,y=scale*df['mean'][cond],yerr=scale*df['sem'][cond],color=color, capsize=3, label=cond)
            ax.legend()
        else:
            for i in self.conds:
                ax.errorbar(x=df.index,y=scale*df['mean'][i],yerr=scale*df['sem'][i],color=cm[i],label=i)
                if self.conds.size>1: ax.legend()
        
        ax.set_xlabel('Time Relative to Stimulus (s)')
        if scale==100:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$ (%)')
        else:
            ax.set_ylabel(r'$\frac{\Delta F}{F}$')
        
        if ylim is not None:
            ax.set_ylim(ylim)
        
        if show:
            py.show()
        
        return ax


    def bin_data(self, binsize, save=False):
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
        bins[-1]+=1


        df=self.all_490.copy()
        for i in df.columns.get_level_values('cond'):
            df[i,'bin start (s)']=bins[np.digitize(self.t,bins)-1]

        def mean_f_trapz(y):
            #calculate the mean within each bin by numerical integration
            y=y.dropna()
            return np.trapz(y,x=y.index)/binsize

        df = df.stack('cond').reset_index('cond')
        df['bin start (s)'] = bins[np.digitize(df.index, bins)-1]
        binned = df.pivot_table(index = ['cond', 'bin start (s)'], aggfunc = mean_f_trapz)
        binned = binned.unstack('cond').swaplevel(axis=1).sort_index(axis=1)
        binned_stats = pd.concat([binned.groupby('cond',axis=1).mean(),
                                  binned.groupby('cond',axis=1).sem()], 
                                 keys=['mean','sem'], axis=1)

        if not hasattr(self,'binned'):
            self.binned={}
        self.binned.update({f'{binsize}':binned_stats})

        if save:
            binned_stats.to_csv(os.path.join(self.file_loc,f"binned_{binsize}s {'_'.join(self.mice)}.csv"))
            binned.to_csv(os.path.join(self.file_loc,f"binned_{binsize}s_inds {'_'.join(self.mice)}.csv"))
        
            
        return binned_stats,binned


    def bin_avg(self, start:int, end:int, save=False, pr=True):
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
            avgs=avgs.rename(columns={'AUC':'Mean ∆F/F'})

            if pr:
                print(' ')
                print(avgs)

            if not hasattr(self,'avgs'):
                self.avgs={}
            self.avgs.update({f'_{start}_{end}':avgs})
            
            if save:
                avgs.to_csv(os.path.join(self.file_loc,f"avg_{start}s_{end}s_{'_'.join(self.mice)}.csv"))

            return avgs
    
    def bin_auc(self, start:int, end:int, save=False, pr=True):
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
            aucs=pd.DataFrame(aucs,columns=['AUC'])

            if pr:
                print(' ')
                print(aucs)

            if not hasattr(self,'aucs'):
                self.aucs={}
            self.aucs.update({f'_{start}_{end}':aucs})

            if save:
                aucs.to_csv(os.path.join(self.file_loc,f"auc_{start}s_{end}s_{'_'.join(self.mice)}.csv"))

        except IndexError:
            print('Time stamps are outside the scope of the analysis! Choose a different start and end time or restart the analysis with a different t_endrec.')
            return None
        
        return aucs

    def ind_peak_df_f(self, extrema:str, save=False, pr=True, trap=False):
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
            pk_inds=self.all_490.loc[0:].apply(lambda x: x.idxmax(skipna=True))
        elif extrema=='min':
            pk_inds=self.all_490.loc[0:].apply(lambda x: x.idxmin(skipna=True))
        else:
            print('Unrecognized extrema!')
            return

        
        if ((self.t_endrec-pk_inds)<=20).any(): print(f'Warning! {extrema} is on the edge of the recording')
        pk_inds[(self.t_endrec-pk_inds)<=20]-=20

        peaks={}
        for i,v in pk_inds.items():
            if trap==True:
                y=self.all_490.loc[v-10:v+10][i].dropna()
                peaks.update({i:np.trapz(x=y.index,y=y.values)})
            else:
                peaks.update({i:self.all_490.loc[v][i]})

        peaks=pd.DataFrame(pd.Series(peaks), columns=[f'{extrema} ∆F/F'])
        pk_times=pd.DataFrame(pk_inds, columns=['time (s)'])

        peaks=pd.concat([peaks,pk_times],axis=1)

        if pr:
            print(peaks)
            print('')
        if not hasattr(self,'pks'):
            self.pks={}

        self.pks.update({f'{extrema}':peaks})

        if save:
            df=pd.DataFrame(peaks,columns=[f'{extrema} ∆F/F','time (s)']).T
            df.to_csv(os.path.join(self.file_loc,f"ind_{extrema}_{'_'.join(self.mice)}.csv"))

        return peaks


    def mean_peak_df_f(self, extrema:str, save=False, pr=True, trap=False):
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

        if extrema=='max':
            pk_inds=self.mean_490.loc[0:].apply(lambda x: x.idxmax(skipna=True))
        elif extrema=='min':
            pk_inds=self.mean_490.loc[0:].apply(lambda x: x.idxmin(skipna=True))
        else:
            print('Unrecognized extrema!')
            return

        if ((self.t_endrec-pk_inds)<=20).any(): print(f'Warning! {extrema} is on the edge of the recording')
        pk_inds[(self.t_endrec-pk_inds)<20]-=20

        peaks={}
        for i,v in pk_inds.items():
            if trap==True:
                y=self.all_490.loc[v-10:v+10][i].dropna()
                y=y.apply(lambda x: np.trapz(x=x.dropna().index,y=x.dropna().values))
            else:
                y=self.all_490.loc[v][i].dropna()
            y.index=y.index.map(lambda x: (i,x))
            peaks.update(y.to_dict())
        peaks=pd.DataFrame(pd.Series(peaks), columns=[f'{extrema} ∆F/F'])

        ts=pk_inds.copy()
        peaks['ts']=ts.loc[peaks.index.get_level_values(0)].values

        if pr:
            print('')
            print(peaks)

        if not hasattr(self,'mean_pks'):
            self.mean_pks={}

        self.mean_pks.update({f'{extrema}':peaks})

        if save:
            peaks.to_csv(os.path.join(self.file_loc,f"mean_{extrema}_{'_'.join(self.mice)}.csv"))

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
    if not hasattr(a, 'ds_freq'):
        a.ds_freq = 1
    if not hasattr(a, 't_prestim'):
        a.t_prestim = a._t_prestim
        a.t_endrec= a._t_endrec
    if not hasattr(a, 'detrend'):
        a.detrend = False
        a.detrend_method = detrend_405_constrained
    if not hasattr(a, 'norm_method'):
        a.norm_method = a._norm_method

    if not isinstance(a.raw_data, pd.Series):
        tmp =  pd.Series([], index = analysis.multi_index, dtype = object)
        for x in a.raw_data:
            if not hasattr(x, 'cond'): x.cond = 0
            if not hasattr(x, 'trial'): 
                get_trial=lambda x: x.trial if hasattr(x,'trial') else -1
                get_mouse_raw = lambda mouse: list(filter( lambda x: x.mouse_id==mouse, a.raw_data))
                mi=list(map( get_trial, get_mouse_raw(x.mouse_id) ))
                x.trial=max(mi)+1
            tmp[x.cond, x.mouse_id, x.trial] = x
        a.raw_data = tmp 

     # make sure the excluded raw data field is formatted correctly
    if not isinstance(a.excluded_raw, pd.Series):
        tmp =  pd.Series([], index = analysis.multi_index, dtype = object)
        for x in a.excluded_raw:
            if not hasattr(x, 'cond'): x.cond = 0
            if not hasattr(x, 'trial'): 
                get_trial=lambda x: x.trial if hasattr(x,'trial') else -1
                get_mouse_raw = lambda mouse: list(filter( lambda x: x.mouse_id==mouse, a.excluded_raw))
                mi=list(map( get_trial, get_mouse_raw(x.mouse_id) ))
                x.trial=max(mi)+1
            tmp[x.cond, x.mouse_id, x.trial] = x
        a.excluded_raw = tmp

    a.compute()
    return a