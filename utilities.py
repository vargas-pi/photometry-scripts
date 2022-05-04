import json
import numpy as np
from scipy import stats as st
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from datetime import timedelta


output_filename='exported_data'

class mouse_data:
    def __init__(self,mouse_id,F490,F405,fs,t_start=None,
                 t=np.array([]),t_stim=None, cond = 0):

        self.mouse_id = mouse_id
        self.F490 = F490
        self.F405 = F405
        self.fs = fs
        self.n = len(F490)
        self.t = t if t.any() else np.arange(0,(self.n)/fs,1/fs)
        if t_start: self.t_start = t_start
        if t_stim: self.t_stim = t_stim
        self.cond = cond
        self.centered = False
   
    @property
    def pre_stim_490(self):
        if self.centered: return self.F490[self.t<0].copy()
        else:
            raise NotCentered
    @property
    def pre_stim_405(self):
        if self.centered: return self.F405[self.t<0].copy()
        else:
            raise NotCentered

    def center_stim(self):
        self.t -= self.t_stim
        if hasattr(self, 't_start'):
            self.t_start = self.t_start + self.t_stim*timedelta(seconds=1)
        self.t_stim = 0
        self.centered = True

    def resample(self, tn, kind='linear'):
        f = interp1d(self.t, self.F490, kind, 
                     bounds_error=False,
                     fill_value='extrapolate')
        self.F490 = f(tn)
        f = interp1d(self.t, self.F405, kind,
                     bounds_error=False,
                     fill_value='extrapolate')
        self.F405 = f(tn)
        self.t = tn.copy()
        self.n = self.t.size
        self.fs = 1/(self.t[1:] - self.t[:-1]).mean()
    
    def exclude(self, ex):
        excl = lambda x: (x > np.mean(x) + ex*np.std(x)) | (x < np.mean(x) - ex*np.std(x))
        mask = excl(self.F405) | excl(self.F490)
        self.F490[mask] = np.nan
        self.F405[mask] = np.nan
    


class NoStimTime(Exception):
    """
    no stimulus time has been specified for this mouse
    """
    pass
class NotCentered(Exception):
    """
    tried to access pre-stim data of an uncentered mouse_data object
    """
    pass

class Encoder(json.JSONEncoder):
    """json encoder for all json files generated/used in the pipeline"""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return {'ndarray': obj.tolist()}
        if isinstance(obj,Callable):
            return {'function':obj.__name__}
        if isinstance(obj,mouse_data):
            return {'mouse_data':obj.__dict__}
        return json.JSONEncoder.default(self, obj)

def hook(obj):
    """hook for loading data from json files generated by this pipeline"""
    if 'ndarray' in obj:
        return np.array(obj['ndarray'])
    if 'mouse_data' in obj:
        return mouse_data(**obj['mouse_data'])
    return obj


####################
#DETRENDING FUNCTIONS
#####################

def fit_405(data, constrained):
    """
    """
    
    #downsample temporarily for efficiency and roughly low pass filter for a cleaner fit
    #TODO: change this to a butterworth filter instead
    filt_ds490 = gaussian_filter1d(data.F490[::int(data.fs)], sigma = 10)
    filt_ds405 = gaussian_filter1d(data.F405[::int(data.fs)], sigma = 10)

    #fit the 405 data to the 490
    if constrained:
        (m, b), _ = curve_fit(lambda x, m, b: np.abs(m)*x + b, filt_ds405, filt_ds490) 
        f0_ds = np.abs(m) * filt_ds405 + b
    else:
        m, b = np.polyfit(filt_ds405, filt_ds490, 1)
        f0_ds = m * filt_ds405 + b
    
    # upsample the baseline estimate
    f0_fn = interp1d( data.t[::int(data.fs)], f0_ds, 
                      bounds_error = False, 
                      fill_value = 'extrapolate')
    f0 = f0_fn(data.t)

    return f0, m, b

    
def detrend_405_constrained(data:mouse_data):
    """
    detrend the 490 and 405 data by estimating a linear fit between the signals
    and constraining the slope to be positive

    Parameters
    ---------
    rec: mouse_data
        an instance of a mouse_data object storing the raw data for a
        given mouse
    t_endrec: float
        the amount of time in seconds from stim to the end of the recording 
        to keep in the analysis

    Returns
    -------
    normed_490: numpy.ndarray
        1D array of normalized 490 data 
    normed_405: numpy.ndarray
        1D array of normalized 405 data 
    t: np.ndarray
        updated time array for the normalized data
    """

    f0, m, b= fit_405(data, constrained = True)
    f0 -= np.median(f0[data.t<0])

    normed_490 = data.F490 - f0
    #map 405 to 490 and subtract the baseline
    normed_405 = (data.F405*np.abs(m) + b) - f0
    #map back to 405 space for further processing
    normed_405 = (normed_405 - b)/np.abs(m)
    return normed_490, normed_405


def detrend_405(data:mouse_data):
    """
    detrend the 490 and 405 data by estimating a linear fit between the signals

    Parameters
    ---------
    rec: mouse_data
        an instance of a mouse_data object storing the raw data for a
        given mouse
    t_endrec: float
        the amount of time in seconds from stim to the end of the recording 
        to keep in the analysis

    Returns
    -------
    normed_490: numpy.ndarray
        1D array of normalized 490 data 
    normed_405: numpy.ndarray
        1D array of normalized 405 data 
    t: np.ndarray
        updated time array for the normalized data
    """

    f0, m, b = fit_405(data, constrained = False)
    f0 -= np.median(f0[data.t<0])

    normed_490 = data.F490 - f0
    #map 405 to 490 and subtract the baseline
    normed_405 = (data.F405 * m + b) - f0 
    #map back to 405 space for further processing
    normed_405 = (normed_405 - b)/m

    return normed_490, normed_405


########################
#NORMALIZATION FUNCTIONS
########################

def norm_to_median_pre_stim(data:mouse_data):
    """
    normalize to the median of the data before the stimulus. we assume the data
    has also been cropped 

    Parameters
    ---------
    rec: mouse_data
        an instance of a mouse_data object storing the raw data for a
        given mouse
    t_endrec: float
        the amount of time in seconds from stim to the end of the recording 
        to keep in the analysis

    Returns
    -------
    normed_490: numpy.ndarray
        1D array of normalized 490 data 
    normed_405: numpy.ndarray
        1D array of normalized 405 data 
    t: np.ndarray
        updated time array for the normalized data
    """

    #compute the baseline by taking the median of the 5 miutes pre-stimu data
    f490_baseline = np.median(data.pre_stim_490)
    f405_baseline = np.median(data.pre_stim_405)

    #normalize the 490 and 405 to the respctive baseline
    normed_490 = (data.F490 - f490_baseline)/f490_baseline
    normed_405 = (data.F405 - f405_baseline)/f405_baseline
    return normed_490, normed_405


def norm_to_405(data:mouse_data):

    f0, m, b = fit_405(data, constrained = False)
    normed_490 = (data.F490 - f0)/f0
    normed_405 = (data.F405 * m + b - f0)/f0

    return normed_490, normed_405

def norm_to_405_constrained(data:mouse_data):

    f0, m, b = fit_405(data, constrained = True)
    normed_490 = (data.F490 - f0)/f0
    normed_405 = (data.F405 * m + b - f0)/f0

    return normed_490, normed_405


def zscore(data:mouse_data):
    normed_490 = (data.F490 - data.F490.mean())/data.F490.std()
    normed_405 = (data.F405 - data.F405.mean())/data.F405.std()
    return normed_490, normed_405

def zscore_mod(data:mouse_data):
    normed_490 = (data.F490 - data.pre_stim_490.mean())/data.pre_stim_490.std()
    normed_405= (data.F405 - data.pre_stim_405.mean())/data.pre_stim_405.std()
    return normed_490, normed_405