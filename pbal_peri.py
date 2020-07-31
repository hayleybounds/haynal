import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import pandas as pd
import sys

from scipy.stats import f_oneway
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from PSTH_creation_utils import *
from chrome_test_holo_utils_2019 import *
from fit_and_bootstrap import basic_fit, find_sat_from_fit
from gen_utils import mean_of_triwise_mean
import sklearn
import h5py

def get_stim_times_peri(file_location):
    with h5py.File(file_location + 'data.mat', 'r') as f:
        stimTimes=f['firstTimes'][:]
    return stimTimes
        
def get_powers_peri(file_location):
    with h5py.File(file_location + 'data.mat', 'r') as f:
        powers = np.array(f['powers'])
    return powers

def get_file_starts_peri(file_location):
    fs = sio.loadmat(glob(file_location + '*file_ends.mat')[0])
    file_ends = fs['file_ends'][0,:]
    print('File ends are', file_ends)
    file_starts = []
    for i,f in enumerate(file_ends):
        if i==0:
            file_starts.append(1)
        else:
            file_starts.append(file_starts[len(file_starts)-1]+f)
    return file_starts

def fix_stim_times_peri(traces, stimTimes):
    #fix if stimtimes doesn't fill in the end nans
    if traces.shape[0]>stimTimes.shape[0]:
        anStimTimes=np.zeros((1,traces.shape[0]))
        anStimTimes[0,0:stimTimes.shape[0]] = stimTimes
        stimTimes=anStimTimes[0,:]
    return stimTimes

#util functions for transforming loaded data
def get_psf_per_online(traces, file_starts, first_stim_times, pre_time, length, dfof_method,
                   df_range = None, do_zscore=False, do_baseline = False, baseline_n=None, skip_ends = 0,
                       skip_starts = 0):
    
    starts = file_starts
    lengths = file_starts.copy()
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, starts, skip_ends, skip_starts)
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    
    #initialize psf array
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    
    for i in range(len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        for cell in range(traces.shape[0]):
            real_start = start+first_stim_times[cell].astype('int16')-pre_time
            real_end = real_start+length
            cut_traces[i,cell,:] = get_traces_segment(real_start, real_end, traces,cell=[cell])
    
    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces

def traces_to_stim_per_online(traces, powers):
    """from the output of get_psf_sing_holo_seq, create a pandas dataframe w/various
    info.
    traces: output of get_psf_sing_holo_seq
    mat_path: path to stim info, has fields 'stimTags' and 'strengthKey' which is ordered list of
        what strength numbers in 'stimTags' corresponds to.
    do_cell_variables: whether to add various cell variables like redness to df
    """
    #change traces into a pandas array with cell and time as cols
    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))
    
    for trial in all_melted['minor'].unique():
        all_melted.loc[all_melted.minor==trial,'stim'] = powers[trial]
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})    
    return all_melted

def label_stimmed_peri(df, stimTimes):
    #label cells actually stimmed
    df = df[~np.isnan(df.stim)]
    cells_stimmed = np.where(stimTimes>0)[0]
    df.loc[:,'stimmed'] = df.cell.isin(cells_stimmed)

    #make sure the targets are actually properly transferred
    assert len(df.cell[df.stimmed].unique()) == len(cells_stimmed)
    for t in cells_stimmed:
        assert t in df.cell[df.stimmed].unique()
    return df


def alt_sigmoid_f(t,a,b,r):
    t[t==0] = 0+.0000001
    return a/(1+(t/b)**-r)

alt_sigmoid = lambda t,a,b,r: a/(1+(t/b)**-r)
alt_sigmoid2 = lambda t,b,r: 1/(1+(t/b)**-r)



def get_sat_vals_peri(mean_df, func, sat_changes=False, remove_downtrend=False, bounds=None, dont_mean=True):
    sat_vals=[]
    for cell in mean_df.cell.unique():
        these_vals = mean_df[mean_df.cell==cell]
        these_vals = these_vals.sort_values(by=['stim'])
        stimwise_mean = mean_of_triwise_mean(these_vals)
        #still make sure sorted
        stimwise_mean = stimwise_mean.sort_values(by=['stim'])
        x_mean=stimwise_mean.stim.values
        y_mean=stimwise_mean.df.values
        semy=stimwise_mean['sem'].values
        if y_mean[0] > y_mean[-1]:
            #print("suppressed")
            sat_vals.append(np.nan)
            continue        
        if remove_downtrend:
            downtrend_remove = []
            for _ in range(1,3): #allow chopping of up to two end values if they're off
                if y_mean[-i]+semy[-i] < y_mean[-(i+1)]-semy[-(i+1)]:
                    downtrend_remove.append(x_mean[-i])
            these_vals = these_vals[~these_vals.stim.isin(downtrend_remove)]
            stimwise_mean = stimwise_mean[~stimwise_mean.stim.isin(downtrend_remove)]
        
        if dont_mean:
            x_vals = these_vals.stim.values
            y_vals = these_vals.df.values
        else:
            x_vals = stimwise_mean.stim.values
            y_vals = stimwise_mean.df.values
            
        #then normalize the y
        y_vals = y_vals - y_mean.min()
        y_vals = y_vals/y_mean.max()
        try:
            if bounds is not None:
                sig2opt = basic_fit(func,  x_vals,  y_vals, bounds=bounds)#, p0= estimate_p0_logistic(x_mean,y_mean)[0:2])
            else:
                sig2opt = basic_fit(func,  x_vals,  y_vals)
            #print('this one worked')
        except:
            print('failed cell',cell)
            sat_vals.append(np.nan)
            continue
            
        sv = find_sat_from_fit(sig2opt, func, sat_changes=sat_changes)
        sat_vals.append(sv)
    return sat_vals


def plot_power_peri(mean_df,cell,func,dont_mean=True):
    r2=np.nan
    omean = mean_df
    mmean_df = mean_of_triwise_mean(mean_df)
    sv = mmean_df[mmean_df.cell==cell].sat_val1.unique()[0]
    #sv4 = mean_df[mean_df.cell==cell].sat_val4.unique()[0]
    these_vals = mmean_df[mmean_df.cell==cell]
    x_mean,y_mean,semy = get_normed_x_y_err(these_vals)
    plt.errorbar(x_mean,y_mean, yerr = semy, marker='o', label = 'online',zorder=1)
    try:
        if dont_mean:
            these_vals = omean[omean.cell==cell]
            x_mean = these_vals.stim.values
            y_mean = these_vals.df.values
            y_mean = y_mean - y_mean.min()
            y_mean = y_mean/y_mean.max()
        
        sig2opt = basic_fit(func,  x_mean,  y_mean)#, p0= estimate_p0_logistic(x_mean,y_mean)[0:2])
        these_vals = omean[omean.cell==cell]
        
        r2=sklearn.metrics.r2_score(y_mean,func(x_mean, *(sig2opt)))
        plot_curve(sig2opt, func, np.linspace(0,x_mean.max()))
        
    except:
        pass
    plt.scatter(sv,1,c='c',marker='*',s=200,zorder=3)
    #plt.scatter(sv4,1,c='k',marker='*',s=200,zorder=3,alpha=.4)
    return r2

def plot_curve(opts, func, xrange):
    plt.plot(xrange, func(xrange,*opts))