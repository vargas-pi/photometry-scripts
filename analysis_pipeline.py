import json
import numpy as np
import matplotlib.pyplot as py
from import_plot import mouse_data, Encoder, hook
import os
from scipy import stats as st
import pandas as pd
import math
from typing import Dict
from pathlib import Path
from scipy.io import savemat

class analysis:

    def __init__(self, norm_method,t_endrec,ex=4,file_format='npy',file_loc=None):
        self.norm_method=norm_method
        self.t_endrec=t_endrec
        self.raw_data= []
        self.normed_data = []
        self.excluded_raw = []
        self.excluded_normed = []
        self.loaded=False
        self.ex=ex

        if not file_loc:
            self.file_loc='analyses'
        else:
            self.file_loc=file_loc

        if file_format in ['json','npy']:
            self.file_format=file_format
        else:
            raise Exception('Unrecognized File Format')
    

    def load_append_save(self,file:Path,t_stims:Dict[str,float]):
        """
        allow users to specify the address of a data file exported from python
        load in this data, normalize it and allow the user to decide if they'd like to keep it in the analysis
        for all data kept in the analysis we then characterize it and save to the instance for later use
        """

        if file.suffix=='.json':
            with open(file,'r') as f:
                print('loading data from file...')
                d=json.load(f,object_hook=hook)
        elif file.suffix=='.npy':
              d=np.load(file,allow_pickle=True).tolist()
        else:
            raise Exception('Unrecognized File Format!')
        for i in range(len(d)):
            d[i].t_stim=t_stims[d[i].mouse_id]
            m=self.normalize_downsample(d[i])
            self.raw_data.append(d[i])
            self.normed_data.append(m)
        self.loaded=True
        self.compute_means()
        self.save()


    def compute_means(self):
        """
        compute the average and standard error accross recordings
        """
        if not self.loaded:
            print('no data has been loaded')
            return
        self.all_490=np.array([rec.F490 for rec in self.normed_data])
        self.all_405=np.array([rec.F405 for rec in self.normed_data])

        self.stim_ind=int(np.where(self.t==0)[0])

        self.mean_490=np.nanmean(self.all_490,axis=0)

        err_490=st.sem(self.all_490,axis=0,nan_policy='omit')
        if isinstance(err_490,np.ma.core.MaskedArray):
            err_490.data[err_490.mask]=np.nan
            self.err_490=err_490.data
        else:
            self.err_490=err_490

        self.mean_405=np.nanmean(self.all_405,axis=0)

        err_405=st.sem(self.all_405,axis=0,nan_policy='omit')
        if isinstance(err_405,np.ma.core.MaskedArray):
            err_405.data[err_405.mask]=np.nan
            self.err_405=err_405.data
        else:
            self.err_405=err_405


    def update_params(self,length=False,norm_method=False,ex_crit=False,file_format=False):
        """
        update parameters set in the beginning of the analysis and recompute
        """

        if not (length or norm_method or ex_crit or file_format):
            return
        
        self.t_endrec = length if length else self.t_endrec
        self.norm_method= norm_method if norm_method else self.norm_method
        self.ex=ex_crit if ex_crit else self.ex

        if file_format in ['json','npy']:
            self.file_format=file_format

        
        for i in range(len(self.raw_data)):
            self.normed_data[i]=self.normalize_downsample(self.raw_data[i])
            
        if hasattr(self,'excluded_raw'):
            if len(self.excluded_raw)>0:
                for i in range(len(self.excluded_raw)):
                    self.excluded_normed[i]=self.normalize_downsample(self.excluded_raw[i])

        self.compute_means()

    def remove_mouse(self,mouse:str):
        """
        remove a mouse from the analysis
        """

        excl_raw=list(filter(lambda x: x.mouse_id==mouse, self.raw_data))
        excl_norm=list(filter(lambda x: x.mouse_id==mouse, self.normed_data))

        if not hasattr(self,'excluded_raw'):
            self.excluded_raw=[]
            self.excluded_normed=[]
        
        self.excluded_raw.extend(excl_raw)
        self.excluded_normed.extend(excl_norm)
        
        self.raw_data=list(filter(lambda x: x.mouse_id!=mouse, self.raw_data))
        self.normed_data=list(filter(lambda x: x.mouse_id!=mouse, self.normed_data))

        if len(self.raw_data)==0:
            self.loaded=False

        self.compute_means()

    def retrieve_excluded(self,mouse:str):
        """
        retrieve an excluded mouse from the analysis
        """
        if not hasattr(self,'excluded_raw') or len(self.excluded_raw)==0:
            print('this analysis has no excluded data')
            return

        ret_raw=list(filter(lambda x: x.mouse_id==mouse, self.excluded_raw))
        ret_norm=list(filter(lambda x: x.mouse_id==mouse, self.excluded_normed))
        self.raw_data.extend(ret_raw)
        self.normed_data.extend(ret_norm)
        self.excluded_raw=list(filter(lambda x: x.mouse_id!=mouse, self.excluded_raw))
        self.excluded_normed=list(filter(lambda x: x.mouse_id!=mouse, self.excluded_normed))
        self.loaded=True
        self.compute_means()

    def save(self,file_format=False):
        """
        export the analysis object to a file in the output folder
        if the output folder hasn't been created, this will create it
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
    
    def export_to_mat(self,):
        data=self.all_490.copy()
        if np.isnan(data).any():
            #check if there are nans and forward fill them first
            r,c=np.where(np.isnan(data))
            for i,j in zip(r,c):
                data[i,j]=data[i,j-1] if j>0 else 0
        
        savemat(os.path.join(self.file_loc,'analysis_'+'_'.join([r.mouse_id for r in self.raw_data])+'.mat'),{'data':data})


    def normalize_downsample(self,rec:mouse_data,ex=4):
        """
        call the speicified normalization method and downsample to ~1 Hz'
        NOTE: the sampling rate from synapse isn't a whole number so we can't necessarilly get exaclt 1Hz but it's close
        """
        normed_490, normed_405, start_ind, stim_ind, end_ind = self.norm_method(rec,self.t_endrec)
        t=rec.t[start_ind:end_ind]-rec.t[stim_ind]

        #resample the data to 1Hz
        #NOTE: when downsampling we want the stimulus to be at t=0
        n_stim_ind=np.where(t==0)[0][0]
        ds=lambda x: np.append( np.flip( x[ n_stim_ind::-math.floor(rec.fs) ][1:] ),
                                x[ n_stim_ind::math.floor(rec.fs) ]     )
        ds_490=ds(normed_490)
        ds_405=ds(normed_405)
        ds_t=ds(t)

        self.t=ds_t.copy()

        #exclude outliers
        excl=lambda x: (x>np.mean(x)+self.ex*np.std(x))|(x<np.mean(x)-self.ex*np.std(x))
        mask=excl(ds_405)
        
        ds_t[mask]=np.nan
        ds_490[mask]=np.nan
        ds_405[mask]=np.nan


        #create a new mouse_data object for the cleaned data
        m=mouse_data(rec.mouse_id,ds_490,ds_405,1)
        m.t=ds_t #overwrite the default time vectors generated by the init fcn
        m.t_stim=rec.t_stim

        #plot the data
        py.plot(m.t,m.F490,'g',linewidth=0.5)
        py.plot(m.t,m.F405,'r',linewidth=0.5)
        py.axvline(x=0, c='k',ls='--', alpha=0.5)
        py.xlabel('Time Relative to Stimulus (s)')
        py.ylabel(r'$\frac{\Delta F}{F}$ (%)')
        py.show()

        return m

    @staticmethod
    def select_data(data:mouse_data,t_endrec):
        """
        select the data that will be used in the analysis'
        """

        
        start_ind=math.ceil((data.t_stim-300)*data.fs)+1
        stim_ind=math.ceil((data.t_stim)*data.fs)+1
        end_ind=math.ceil((data.t_stim+t_endrec)*data.fs)+1

        pre_stim_490=data.F490[start_ind:stim_ind]
        pre_stim_405=data.F405[start_ind:stim_ind]
        sel_490=data.F490[start_ind:end_ind+1]
        sel_405=data.F405[start_ind:end_ind+1]


        return (start_ind, stim_ind, end_ind, pre_stim_490, pre_stim_405, sel_490, sel_405)
    
    @staticmethod
    def norm_to_median_pre_stim(data:mouse_data,t_endrec):
        """
        normalize to the median of the data 5 minutes before the stimulus'
        """

        start_ind, stim_ind, end_ind, pre_stim_490, pre_stim_405, sel_490, sel_405=analysis.select_data(data,t_endrec)

        #compute the baseline by taking the median of the 5 miutes pre-stimu data
        f490_baseline=np.median(pre_stim_490)
        f405_baseline=np.median(pre_stim_405)

        #normalize the 490 and 405 to the respctive baseline
        normed_490=(sel_490-f490_baseline)/f490_baseline
        normed_405=(sel_405-f405_baseline)/f405_baseline
        return (normed_490,normed_405, start_ind, stim_ind, end_ind)

    @staticmethod
    def norm_to_405(data:mouse_data,t_endrec):
        """
        normalize to a linear fit of the 405 data'
        """
        start_ind, stim_ind, end_ind, pre_stim_490, pre_stim_405, sel_490, sel_405=analysis.select_data(data,t_endrec)

        x=np.arange(len(sel_405))
        m, b, _, _, _ = st.linregress(x, sel_405)

        fit_405=m*x+b

        normed_490=100*np.divide(sel_490-fit_405, fit_405)
        normed_490-=np.median(normed_490[0:300])

        normed_405=100*np.divide(sel_405-fit_405, fit_405)
        normed_405-=np.median(normed_405[0:300])
        
        return (normed_490,normed_405, start_ind, stim_ind, end_ind)

    @staticmethod
    def zscore(data:mouse_data):
        pass 

    def plot_both(self):
        """
        plot the average normalized 490 and 405 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        py.fill_between(self.t, self.mean_490 + self.err_490,
                        self.mean_490 - self.err_490, color='g',alpha=0.2 )
        py.fill_between(self.t, self.mean_405 + self.err_405,
                        self.mean_405 - self.err_405, color='r',alpha=0.2  )
        py.plot(self.t, self.mean_490 , 'g', linewidth=0.5)
        py.plot(self.t, self.mean_405, 'r', linewidth=0.5)
        py.axvline(x=0, c='k',ls='--', alpha=0.5)
        py.xlabel('Time Relative to Stimulus (s)')
        py.ylabel(r'$\frac{\Delta F}{F}$ (%)')

        py.show()
        
    def plot_490(self):
        """
        plot the average normalized 490 signal with error
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return
        
        py.fill_between(self.t, self.mean_490 + self.err_490,
                        self.mean_490 - self.err_490, color='g',alpha=0.2 )
        py.plot(self.t, self.mean_490 , 'g', linewidth=0.5)
        py.axvline(x=0, c='k',ls='--', alpha=0.5)
        py.xlabel('Time Relative to Stimulus (s)')
        py.ylabel(r'$\frac{\Delta F}{F}$ (%)')

        py.show()

    def bin_plot(self,binsize:int):
        """
        bin the data average the data in each bin for each mouse
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        #TODO: account for binsizes that aren't evenly divisible

        bins=np.append(-300,self.t[self.stim_ind-binsize::binsize])
        binsize_act=bins[-1]-bins[-2]

        added=False

        if self.t[-1] not in bins:
            bins=np.append(bins,self.t[-1]+1)
            added=True
        else:
            bins[-1]+=1

        all_490=pd.DataFrame(self.all_490.transpose())
        for i in range(len(self.normed_data)): all_490[i]=list(zip(all_490[i],self.t))


        all_490['bins']=np.digitize(self.t,bins)
        all_490['border']=np.append(0,np.array(all_490['bins'][1:])-np.array(all_490['bins'][0:-1])==1)

        borders=all_490[all_490.border==1].copy()
        borders.bins-=1
        borders=borders.set_index(borders.index+0.5) #not sure if this is necessary

        all_490.update(borders)

        def nantrapz(x):
            y,x2=np.array(list(zip(*x)))
            y,x2=y[~np.isnan(y)],x2[~np.isnan(y)],
            return np.trapz( y, x=x2 )/(binsize_act)

        binned=pd.pivot_table(all_490,values=range(len(self.normed_data)),index='bins', aggfunc=nantrapz)
        binned['bins']=bins[0:-1]
        binned['binsizes']=bins[1:]-bins[0:-1]
        binned=binned[binned.bins!=-300] #remove the first bin
        binned=binned[:-1] if added else binned #only include bins of the appropriate size (i.e. get rid the end bin if needed)

        print(list(binned.columns))
        binned['mean']=binned[range(len(self.normed_data))].mean(axis=1,skipna=True) 
        binned['sem']=binned[range(len(self.normed_data))].sem(axis=1,skipna=True)

        binned=binned.rename( columns={ i: self.normed_data[i].mouse_id for i in range(len(self.normed_data)) } )

        print(binned)

        py.errorbar(binned['bins'],binned['mean'],yerr=binned['sem'])
        py.xlabel('Time Relative to Stimulus (s)')
        py.ylabel(r'$\frac{\Delta F}{F}$ (%)')
        py.show()

        if not hasattr(self,'binned'):
            self.binned={}
        self.binned.update({f'{binsize}':binned.to_dict()})

        return binned

    def bin_avg(self,start:int,end:int):
        """
        average over a specified section of data for each mouse
        """
        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        nearest= lambda arr,val: np.abs(arr-val).argmin()

        start_t=nearest(self.t,start)
        end_t=nearest(self.t,end)


        try: 
            y=self.all_490[:,start+self.stim_ind:end+self.stim_ind+1]
            x=self.t[start+self.stim_ind:end+self.stim_ind+1]
            mask= ~(np.isnan(y).max(axis=0))
            y=y[:,mask]
            x=x[mask]

            avg=np.trapz(y, x=x,axis=1)/(end_t-start_t)
            avg_dict={}
            print(' ')
            for r,av in zip(self.normed_data,avg): 
                avg_dict.update({r.mouse_id:av})
                print(f'{r.mouse_id}:{av}')

            if not hasattr(self,'avgs'):
                self.avgs={}
            self.avgs.update({f'_{start}_{end}':avg_dict})
            
            
        except IndexError:
            print('Time stamps are outside the scope of the analysis! Choose a different start and end time or restart the analysis with a different t_endrec.')
        
        return avg_dict
    
    def ind_peak_df_f(self,extrema:str):
        """
        determine either the min or max ∆f/f for each mouse
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return


        if extrema=='max':
            pks_ind=300+np.nanargmax(self.all_490[:,300:],axis=1)
        elif extrema=='min':
            pks_ind=300+np.nanargmin(self.all_490[:,300:],axis=1)
        else:
            print('Unrecognized extrema!')

        for i,v in enumerate(pks_ind):
            #check that the extrema aren't on the edge of the recording
            if self.all_490.shape[-1]-20<v:
                pks_ind[i]-=10
                print(f'Warning! {extrema} is on the edge of the recording')

        surr_inds=np.arange(-5,6)+pks_ind[:,np.newaxis]
        row_inds=np.arange(self.all_490.shape[0]).repeat(11)
        x=self.t[surr_inds.flatten()].reshape((-1,11))
        y= self.all_490[row_inds,surr_inds.flatten()].reshape((-1,11))


        peaks=[np.trapz(y=i[~np.isnan(i)], x=j[~np.isnan(i)])/10 for i,j in zip(y,x)]

        pks_dict={}

        print('')
        for r,p,l in zip(self.normed_data,peaks,x[:,5].flatten()):
            pks_dict.update({r.mouse_id:p})
            print(f'{r.mouse_id}:{p}, at {l}s')

        if not hasattr(self,'pks'):
            self.pks={}

        self.pks.update({f'{extrema}':pks_dict})

        return pks_dict


    def mean_peak_df_f(self,extrema:str):
        """
        determine either the min or max ∆f/f for each mouse
        """

        if not self.loaded:
            print('Must have usable data loaded in the analysis first!')
            return

        pk_ind= 300+np.nanargmax(self.mean_490[300:]) if extrema=='max' else 300+np.nanargmin(self.mean_490[300:])
        if self.all_490.shape[-1]-20<pk_ind:
                pk_ind-=10
                print(f'Warning! {extrema} is on the edge of the recording')

        x=self.t[pk_ind-5:pk_ind+6]
        y= self.all_490[:,pk_ind-5:pk_ind+6]
        mask= ~(np.isnan(y).max(axis=0))
        x,y= x[mask],y[:,mask]

        peak=np.trapz(y,x=x,axis=1)/10
        pks_dict={}

        print('')
        for r,p in zip(self.normed_data,peak):
            pks_dict.update({r.mouse_id:p})
            print(f'{r.mouse_id}:{p}, at {x[5]}s')

        if not hasattr(self,'mean_pks'):
            self.mean_pks={}

        self.mean_pks.update({f'{extrema}':pks_dict})


    def time_to_half_pk(self):
        raise NotImplementedError

class loaded_analysis(analysis):
    def __init__(self,fpath):
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
    fpath=Path(fpath).resolve()
    if fpath.suffix=='.json':
        a=loaded_analysis(fpath)
    elif fpath.suffix=='.npy':
        a=np.load(fpath,allow_pickle=True).item()
    else:
        raise Exception('Unrecognized file format!')
    
    a.file_loc=fpath.parent
    return a