import numpy as np
import pandas as pd
import scipy
from scipy import stats
import multiprocessing as mp
import time
import warnings
from scipy.optimize import curve_fit

def mean_confidence_interval(data, confidence=95):
    data = np.asarray(data)
    data = data[~np.isnan(data)]
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    low = np.percentile(data,(100-confidence)/2)
    
    high = np.percentile(data,100-(100-confidence)/2)
    return m, low, high

def single_bootstrap(x,y,stims,n_samples,func,sat_changes,remove_downtrend, bounds, force_baseline, do_max_adj,
                    remove_outliers,take_means,with_replacement=False,subt_min=True):
    warnings.filterwarnings(action='ignore')
    np.seterr(all='ignore')
    this_y = []
    this_stims = []
    for st in stims:
        
        #
        if with_replacement:
            this_x = np.where(x==st)[0][np.random.randint(0,len(np.where(x==st)[0]),(n_samples))]
        else:
            this_x = np.random.permutation(np.where(x==st)[0])[:n_samples]
        an_y = y[this_x]
        if remove_outliers:
            z = np.abs(stats.zscore(an_y))
            an_y = an_y[z<2.2]
        if take_means:
            this_y.append(np.nanmean(an_y))
            this_stims = stims
        else:
            this_y = this_y + list(an_y)
            this_stims = this_stims + list(x[this_x])
    this_y = np.asarray(this_y)
    this_stims = np.asarray(this_stims)
    this_y, this_stims = process_data(this_y, this_stims, force_baseline, do_max_adj, remove_downtrend,subt_min)
    
    try:
        if bounds is not None:
            sig2opt = basic_fit(func,  this_stims,  this_y, bounds = bounds)#, p0= estimate_p0_logistic(x_mean,y_mean)[0:2])
        else:
            sig2opt = basic_fit(func,  this_stims,  this_y)
        return find_sat_from_fit(sig2opt,func,sat_changes=sat_changes)
    except:
        #print(time.time()-s)
        return np.nan
    
    

def bootstrap_10(x,y,stims,n_samples,func,sat_changes,remove_downtrend, bounds, force_baseline, do_max_adj,
                    remove_outliers,take_means,with_replacement,subt_min, r_seed):
    warnings.filterwarnings(action='ignore')
    np.seterr(all='ignore')
    np.random.seed(r_seed)
    #s=time.time()
    svs = [single_bootstrap(x,y,stims,n_samples,func,sat_changes,remove_downtrend, bounds, force_baseline, do_max_adj,
                    remove_outliers,subt_min,take_means) for _ in range(10)]
    
    #print(time.time()-s)
    return svs

def process_data(this_y,stims, force_baseline, do_max_adj, remove_downtrend,subt_min):
    this_y = np.asarray(this_y)
    if force_baseline:
        this_y[this_y<this_y[0]] = this_y[0]
    if subt_min:
        this_y = this_y -this_y.min()
    if do_max_adj:
        this_y = this_y/this_y.max()
    this_stims=stims.copy()
    if remove_downtrend:
        for _ in range(2): #allow chopping of up to two end values if they're off
            if this_y[-1]+this_y[-1]*.2 < this_y[-2]+this_y[-2]*.2:
                this_stims = this_stims[0:-1]
                this_y = this_y[0:-1]
                
    return this_y, this_stims
    

def find_sat_from_fit(sig2opt,func,sat_changes=True):
    y_vals = func(np.linspace(0,500,501), *sig2opt)
    if sat_changes:
        sat_approach = np.where(y_vals>.85*sig2opt[0])
    else:
        sat_approach = np.where(y_vals>.85)
    if len(sat_approach[0])>0:
        return sat_approach[0][0]
    else:
        return np.nan
    
def alt_sigmoid(t,a,b,r):
    return a/(1+(t/b)**-r)

def alt_sigmoid3(t,a,b,r,c):
    return c+a/(1+(t/b)**-r)
def alt_sigmoid_f(t,a,b,r):
    t[t==0] = 0+.0000001
    return a/(1+(t/b)**-r)
def get_bootstrapped_fits(delta_dff, cell, func, n_samples, sat_changes=True, bounds=None,remove_downtrend=False,
                         parallel=False, force_baseline=False, do_max_adj=True,remove_outliers=False,
                         take_means = False,with_replacement=True, subt_min=True):
    warnings.filterwarnings(action='once')
    np.seterr(over='ignore')
    x = delta_dff.stim[delta_dff.cell==cell].values
    y = delta_dff.differ[delta_dff.cell==cell].values
    #find the minimum number of trials for a stimulus
    stims, counts = np.unique(x, return_counts=True)
    if n_samples > counts.min():
        print('not enough samples')
        return [np.nan]
    sat_vals=[]
    #sample randomly 100 times, taking n_samples from each value
    if parallel:
        pool = mp.Pool(mp.cpu_count())
        
        temp = [(x,y,stims,n_samples,func,
                                                                   sat_changes,remove_downtrend,bounds,
                                     force_baseline,do_max_adj,remove_outliers,take_means, with_replacement,subt_min,
                 np.random.randint(30000))
                                           for zz in range(100)]
        sat_vals = pool.starmap_async(bootstrap_10, temp).get()

        # result_objects is a list of pool.ApplyResult objects
        pool.close()
    else:
        sat_vals = [single_bootstrap(x,y,stims,n_samples,func,
                                     sat_changes,remove_downtrend,bounds,
                                     force_baseline,do_max_adj,remove_outliers,take_means,
                                     with_replacement,subt_min) for _ in range(1000)]
    
    return sat_vals




def basic_fit(func, x, y,bounds=None,p0=None):
    if bounds is not None:
        if p0 is not None:
            opt_params, _ = curve_fit(func,x,y,bounds=bounds,p0=p0)
        else:
            opt_params, _ = curve_fit(func,x,y,bounds=bounds)
    else:
        opt_params, _ = curve_fit(func,x,y)
    return opt_params