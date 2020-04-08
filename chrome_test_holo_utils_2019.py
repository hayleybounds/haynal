
import numpy as np
import scipy.io as sio
from skimage import io
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.stats import sem
from scipy.stats.mstats import zscore
import seaborn as sns
import pandas as pd
import math
from scipy.stats import f_oneway
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import r2_score
from datetime import datetime
from scipy.optimize import curve_fit
from PSTH_creation_utils import *


def mean_by_trial(df, start, stop, add_columns=[]):
    """given df with cell + ƒltrial columns as well as any cols in add_columns,
    takes only values between start time + stop time, groups them by cell,
    trial, and add_columns values (usually add_columns would just be cell-specific
    variables I want to maintain, like redness or something). Then takes the mean
    of all remaining columns per cell, trial + add_columns group
    """
    df = df[(df.time>=start) & (df.time<=stop)]
    cols = ['cell', 'trial'] +add_columns
    mean_on = df.groupby(by=cols).mean()
    
    
    cols = ['cell', 'trial'] + add_columns
    mean_on = mean_on.reset_index(level=cols)
    
    return mean_on

def mean_by_id_with_sem(df,start,stop,col='id',add_columns=[]):
    df = df[(df.time>=start) & (df.time<=stop)]
    cols = ['cell', col] +add_columns
    cols_first = cols+['trial']
    mean_first = df.groupby(cols_first).mean()
    mean_on = mean_first.groupby(by=cols).mean()
    mean_sem = mean_first.groupby(by=cols).df.sem()
    mean_on = mean_on.reset_index(level=cols.append(add_columns))
    mean_on.loc[:,'sem'] = mean_sem.values
    return mean_on

def mean_by_id(df, start, stop, col='id', sem_col='trial',add_columns=[]):
    """given df with cell + id column, 'id', as well as any cols in add_columns,
    takes only values between start time + stop time, groups them by cell,
    id, and add_columns values (usually add_columns would just be cell-specific
    variables I want to maintain, like redness or something). Then takes the mean
    of all remaining columns per cell, id + add_columns group.
    Usually, id is a stim id, so you're getting a mean by stim type
    """
    df = df[(df.time>=start) & (df.time<=stop)]
    cols = ['cell', col] +add_columns
    mean_on = df.groupby(cols).mean()
   
    mean_on = mean_on.reset_index(level=cols.append(add_columns))
    return mean_on

def mean_difference(df, prestart, prestop, start, stop, add_columns = []):
    """get the mean of during stim (start to stop) minus mean of pre stim (prestart to prestop).
    output has columns cell, trial, add_columns, mean_post(mean post stim), mean_pre,
    and differ(mean_post-mean_pre).
    """
    mean_pre = mean_by_trial(df, prestart, prestop, add_columns=add_columns)
    mean_post = mean_by_trial(df,start, stop, add_columns = add_columns)
    mean_post = mean_post.rename(columns = {'df':'mean_post'})
    mean_post['mean_pre'] = mean_pre.df.values
    mean_post['differ'] = mean_post.mean_post - mean_post.mean_pre
    return mean_post


def get_normed_x_y_err(these_vals):
    x_mean=these_vals.stim.values
    y_mean=these_vals.differ.values
    semy=these_vals['sem'].values

    y_mean = y_mean - y_mean.min()
    semy = semy/y_mean.max()
    y_mean = y_mean/y_mean.max()
    return x_mean,y_mean,semy

def mean_of_mean_diffs(mean_post):
    #gives a single value per stim value, rather than per trial

    mean_diffs = mean_post.groupby(['cell','stim']).mean()
    sem_diffs = mean_post.groupby(['cell','stim']).sem()

    mean_diffs = mean_diffs.reset_index(level=['cell','stim'])
    mean_diffs['sem'] = sem_diffs.values
    mean_diffs_no_long = mean_diffs.loc[(mean_diffs.stim!=1000) & (mean_diffs.stim!=500),:]
    #mean_post_no_long = mean_post.loc[(mean_post.stim!=10000) & (mean_diffs.stim!=500),:]
    return mean_diffs_no_long

def label_good_cells(df, cutoff, print_percents = False):
    """given df, add good_cell column that reflects which cells have a dist to target matched less than cutoff"""
    mat = df.file.unique()[0]
    mat = sio.loadmat(mat)
    best_targets = get_index_of_best_targets(mat, dist_cutoff=cutoff)
    df.loc[:,'good_match'] = df.cell.isin(best_targets)
    #make sure the targets are actually properly transferred
    assert len(df.cell[df.good_match].unique()) == len(best_targets)
    for t in best_targets:
        assert t in df.cell[df.good_match].unique()
    
    if print_percents:
        print(df.mouse.values[0], df.dpi.values[0], round(len(best_targets)/len(df.cell.unique()),2))
        #ask how many as a percent of in range targets were identified
        try:
            print(round(len(best_targets)/(len(df.cell.unique())-len(mat['outOfRange'][0])),2))
        except:
            print('MISSING', df.day.values[0])
        return df #you modded the original object anyways, but just in case you want it, it's returned
    
def get_p_stimmable(mean_df,val_col,stim_id_col = 'stim'):
    """use anova to test for each cell in mean_df, whether the responses, given by val_col, are
    significantly different for different stimuli, given by stim_id_col
    """
    ps=[]
    for cell in mean_df.cell.unique():
        #use list comprehension to get a list of all the values for this cell for each stim point
        vals = mean_df[mean_df.cell==cell] #vals for this cell
        stims = vals[stim_id_col].unique()
        [_, p] = f_oneway(*[vals[vals[stim_id_col]==stim][val_col] for stim in stims])
        ps.append(p)
    return ps



#util functions for transforming loaded data
def get_psth_sing_holo_seq(traces, mat_path, stim_starts, stim_spacing, pre_time, length, dfof_method,
                   df_range = None, do_zscore=False, do_baseline = False, baseline_n=None, skip_ends = 0,
                   skip_starts = 0):
    """given traces, gets the peristimulus fluorescence values as n_trials x n_cells x length.
    traces: a n_cell x frames numpy array w/cells ordered in holo target order
    mat_path: path to mat containing FileStruct which has a value fileStart that delims
    starts of files.
    stim_starts: frame number of the first holo stimulus
    stim_spacing: frames between each successive stimulus
    pre_time: n frames before stim to add to psf
    length: total length of psf
    dfof_method: str, one of 'percentile', 'by trial', or 'at start'. percentile is not trialwise,
        so done using all data at start. rcentile takes the 30th 
        prctile value per cell as the f0. by trial does it so that the f0 is the mean of the values
        in df_range per trial, at_start takes the mean of the frames in df_range at the 
        start of the recording
    df_range: used by by trial + at start dfof methods, mean diff things for each of them
    do_zscore: if used, before traces are cut, they are zscored by cell. If dfof_method is percentile, traces
        are first dff'ed
    do_baseline: option to subtract the mean of baseline_n frames from the psf. This is per trial. Done after dff. 
    baseline_n: how many frames, beginning at first frame of each trial (as defined by psf, not from beginning
        of acq, to use for baseline subtraction
    skip_ends: how many files to skip at the end (they sometimes contain errors)
    """
    
    mat = sio.loadmat(mat_path)
    starts = get_trial_starts(mat)
    lengths = get_trial_lengths(mat)
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, lengths, skip_ends, skip_starts)
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    
    #initialize psf array
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    
    for i in range(len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        for cell in range(traces.shape[0]):
            if stim_spacing is None:
                if np.isnan(stim_starts[cell]):
                    real_start = np.nan
                else:
                    real_start = start+math.floor(stim_starts[cell])-pre_time
            else:
                real_start = start+stim_starts+math.floor(cell*stim_spacing)-pre_time
            real_end = real_start+length
            if np.isnan(real_start):
                cut = np.tile(np.nan,length)
            #if this is the first trial and real start is negative, pad it
            elif real_start < 0 and real_end > 0:
                cut = np.concatenate((np.tile(np.nan, (-1*real_start)), traces[cell,0:real_end]),0)
            elif real_end < 0: #this should never happen
                print("UH OH REAL END IS NOT EVEN IN TIME")
                return
            #if real_end longer than available data, but real start isn't
            elif real_end > traces.shape[1] and real_start < traces.shape[1]:
                cut = np.concatenate((traces[cell,real_start:],np.tile(np.nan, (real_end-traces.shape[1]))),0)
            elif real_end > traces.shape[1] and real_start >= traces.shape[1]: #this also should not happen
                print('UH OH REAL START IS LONGER THAN TRACES')
                return
            else:
                cut = traces[cell, real_start:real_end]
            cut_traces[i,cell,:] = cut 
    
    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces


def get_psf_run(traces, mat_path, stim_starts, stim_spacing, pre_time, length, dfof_method,
                   df_range = None, do_zscore=False, do_baseline = False, baseline_n=None, skip_ends = 0):
    """given traces, gets the peristimulus fluorescence values as n_trials x n_cells x length.
    traces: a n_cell x frames numpy array w/cells ordered in holo target order
    mat_path: path to mat containing FileStruct which has a value fileStart that delims
    starts of files.
    stim_starts: frame number of the first holo stimulus
    stim_spacing: frames between each successive stimulus
    pre_time: n frames before stim to add to psf
    length: total length of psf
    dfof_method: str, one of 'percentile', 'by trial', or 'at start'. percentile is not trialwise,
        so done using all data at start. rcentile takes the 30th 
        prctile value per cell as the f0. by trial does it so that the f0 is the mean of the values
        in df_range per trial, at_start takes the mean of the frames in df_range at the 
        start of the recording
    df_range: used by by trial + at start dfof methods, mean diff things for each of them
    do_zscore: if used, before traces are cut, they are zscored by cell. If dfof_method is percentile, traces
        are first dff'ed
    do_baseline: option to subtract the mean of baseline_n frames from the psf. This is per trial. Done after dff. 
    baseline_n: how many frames, beginning at first frame of each trial (as defined by psf, not from beginning
        of acq, to use for baseline subtraction
    skip_ends: how many files to skip at the end (they sometimes contain errors)
    """
    run_path = mat_path.split('holo_raw')
    run_path = run_path[0] + 'holo_run_data'
    runs = sio.loadmat(run_path)['runs']
    mat = sio.loadmat(mat_path)
    starts = get_trial_starts(mat)
    if skip_ends > 0:
        starts = starts[0:len(starts)-skip_ends]
    #percentile dff method is based on whole trace, so do now

    #initialize psf array
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    
    for i in range(0,len(starts)):
        for cell in range(0, traces.shape[0]):
            real_start = stim_starts+math.floor(cell*stim_spacing)-pre_time
            real_end = real_start+length
            #if this is the first trial and real start is negative, pad it
            if real_start < 0 and real_end > 0:
                cut = np.concatenate((np.tile(np.nan, (-1*real_start)), runs[i,0:real_end]),0)
            elif real_end < 0: #this should never happen
                print("UH OH REAL END IS NOT EVEN IN TIME")
                return
            #if real_end longer than available data, but real start isn't
            elif real_end > runs.shape[1] and real_start < runs.shape[1]:
                cut = np.concatenate((runs[i,real_start:],np.tile(np.nan, (real_end-runs.shape[1]))),0)
            elif real_end > runs.shape[1] and real_start >= runs.shape[1]: #this also should not happen
                print(real_end, stim_starts, cell, stim_spacing)
                cut = np.tile(np.nan,length)
                print('UH OH REAL START IS LONGER THAN TRACES')
            
            else:
                cut = runs[i, real_start:real_end]
            cut_traces[i,cell,:] = cut 
   
    return cut_traces

def get_motion_psf(traces, mat_path, stim_starts, stim_spacing, pre_time, length,do_target_order=True,skip_ends=0):
    
    mat = sio.loadmat(mat_path)
    motion = get_motion_data(mat)
    planes = get_cell_planes(mat,target_order = do_target_order)
    starts = get_trial_starts(mat)
    if skip_ends > 0:
        starts = starts[0:len(starts)-skip_ends]
        
    motion = motion**2
    motion = np.sum(motion, axis=1)
    motion = np.sqrt(motion)
    motion = np.squeeze(motion)
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    for i in range(len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        for cell in range(traces.shape[0]):
            cell_plane = planes[cell]
            if stim_spacing is None:
                if np.isnan(stim_starts[cell]):
                    real_start = np.nan
                else:
                    real_start = start+math.floor(stim_starts[cell])-pre_time
            else:
                real_start = start+stim_starts+math.floor(cell*stim_spacing)-pre_time
            real_end = real_start+length
            if np.isnan(real_start):
                cut = np.tile(np.nan,length)
            #if this is the first trial and real start is negative, pad it
            elif real_start < 0 and real_end > 0:
                cut = np.concatenate((np.tile(np.nan, (-1*real_start)), motion[0:real_end,cell_plane-1]),0)
            elif real_end < 0: #this should never happen
                print("UH OH REAL END IS NOT EVEN IN TIME")
                return
            #if real_end longer than available data, but real start isn't
            elif real_end > traces.shape[1] and real_start < traces.shape[1]:
                cut = np.concatenate((motion[real_start:,cell_plane-1],np.tile(np.nan, (real_end-traces.shape[1]))),0)
            elif real_end > traces.shape[1] and real_start >= traces.shape[1]: #this also should not happen
                print('UH OH REAL START IS LONGER THAN TRACES')
                return
            else:
                cut = motion[real_start:real_end,cell_plane-1]
            cut_traces[i,cell,:] = cut 
    return cut_traces


def get_trialwise_motion(traces, mat_path, skip_ends=0, skip_starts=0, do_target_order=True):
    
    mat = sio.loadmat(mat_path)
    motion = get_motion_data(mat)
    planes = get_cell_planes(mat,target_order = do_target_order)
    
    starts = get_trial_starts(mat)
    lengths = get_trial_lengths(mat).astype('int16')
    if traces is None:
        traces = get_traces(mat)
    
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, lengths, skip_ends, skip_starts)
    
    motion = motion**2
    motion = np.sum(motion, axis=1)
    motion = np.sqrt(motion)
    motion = np.squeeze(motion)
    
    #initialize psf array
    length = min(lengths)
    cut_traces = np.zeros((int(np.shape(starts)[0]), len(np.unique(planes)), min(lengths)))
    for i in range(0,len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        for j,plane in enumerate(np.unique(planes)):
            cut_traces[i,j,:] = motion[start:start+length,j-1]
    return cut_traces
            
            
def find_sat_power(func,opt,sat_changes=False,sat_perc = .85):
    y_vals = func(np.linspace(0,1000,1000),*opt)
    if sat_changes:
        sat_approach = np.where(y_vals>sat_perc*opt[0])
    else:
        sat_approach = np.where(y_vals>sat_perc)
    if len(sat_approach[0])>0:
        return sat_approach[0][0]
    else:
        return np.nan
def traces_to_stim_df(traces, mat_path, do_cell_variables = False, is_baseline=False, cut_runs=None, is_new = False):
    """from the output of get_psf_sing_holo_seq, create a pandas dataframe w/various
    info.
    traces: output of get_psf_sing_holo_seq
    mat_path: path to stim info, has fields 'stimTags' and 'strengthKey' which is ordered list of
        what strength numbers in 'stimTags' corresponds to.
    do_cell_variables: whether to add various cell variables like redness to df
    is_baseline: if true, ignore stim intensity info and set all to 0.
    """
    
    mat = sio.loadmat(mat_path)
    
    #change traces into a pandas array with cell and time as cols
    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))
    
    if cut_runs is not None:
        print('will add run data')
    #define trial stim variables
    tags = get_stim_tags(mat)
    if not is_new:
        strengths = get_strength_key(mat)
    if len(all_melted['minor'].unique()) > len(tags):
        raise ValueError('number of stim tags is less than n trials in psf array');
    for trial in all_melted['minor'].unique():
        if cut_runs is not None:
            for cell in all_melted['major'].unique():
                all_melted.loc[(all_melted.minor==trial) & (all_melted.major==cell), 'mean_run'] = np.nanmean(cut_runs[trial,cell,:])
        if is_baseline:
            all_melted.loc[all_melted.minor==trial,'stim'] = 0
            all_melted.loc[all_melted.minor==trial,'minor'] = trial + np.max(all_melted.minor.unique())+1
        elif is_new:
            all_melted.loc[all_melted.minor==trial,'stim'] = tags[trial]
        else:
            all_melted.loc[all_melted.minor==trial,'stim'] = strengths[tags[trial]]
        

    if do_cell_variables:
        #define cell level variables
        redratio = get_cell_stat(mat, 'redratio', target_order = True)
        redcell = get_cell_stat(mat, 'redcell', target_order = True)
        redprob = get_cell_stat(mat, 'redprob', target_order = True)
        try:
            red_val = get_cell_stat(mat, 'red_interior_uncorrected', target_order = True)
        except:
            print("no red val")
            red_val = np.tile(np.nan, red_cell.shape)
        target_id = get_ROI_order(mat)
        target_distance = get_target_distance(mat)
        for cell in all_melted['major'].unique():
            these_cell = all_melted.major==cell
            all_melted.loc[these_cell,'redratio'] = redratio[cell] 
            all_melted.loc[these_cell,'redcell'] = redcell[cell] 
            all_melted.loc[these_cell,'redprob'] = redprob[cell] 
            all_melted.loc[these_cell,'target_id'] = target_id[cell]
            all_melted.loc[these_cell,'target_distance'] = target_distance[cell]
            all_melted.loc[these_cell, 'red_val'] = red_val[cell]

    #define overall variables
    all_melted.loc[:,'file'] = mat_path
    if not is_new:
        mouse, day = get_mouse_and_day(mat)
        all_melted.loc[:,'mouse'] = mouse
        all_melted.loc[:,'day'] = day
        
        mice_inj_dict = {'I119.3':'180802', 'I119.4':'180802', 'HB17.2':'180808', 'HB17.3':'180809',
                         'HB23.1':'180814'}
        dpi = days_between(mice_inj_dict[mouse], day)
        all_melted.loc[:,'dpi'] = dpi
        
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})
    return all_melted


#util functions for online data

def load_and_format_online(online_path, nplanes=3):
    """load online csv from online_path, reindex it + interpolate missing values, 
    and transpose it to be cell x frame
    """
    online_data = pd.read_csv(online_path)
    online_data = online_data.set_index('frameNumber')
    new_idx = list(range(nplanes,online_data.index.max(),nplanes))
    online_data_reindex = online_data.reindex(new_idx, index = 'frameNumber',method='nearest')
    online_data_reindex = online_data_reindex.fillna(method='ffill')
    #turn it into a numpy array, with cells as cols. skip 0 col as thats index
    online_arr=online_data_reindex.values[:,1:]
    online_arr=online_arr.transpose()
    return online_arr

#util functions for accessing parts of the mat file

def get_ordered_traces(mat_path):
    """loads the mat file at mat_path, gets the traces
    and orders them by ROI order.
    """
    mat = sio.loadmat(mat_path)
    traces = get_traces(mat)
    order = get_ROI_order(mat)
    traces = traces[order,:]
    return traces

def get_traces(traces):
    traces = traces['signals']['F'][0][0]
    return traces

def get_mouse_day_dpi(mat):
    mouse, day = get_mouse_and_day(mat)
    mice_inj_dict = {'I119.3':'180802', 'I119.4':'180802', 'HB17.2':'180808', 'HB17.3':'180809',
                    'I119.2':'180801', 'HB23.1': '180814'}
    dpi = days_between(mice_inj_dict[mat['mouse'][0]], mat['day'][0])
    return mouse, day, dpi

def get_trial_starts(ends):
    starts = ends['fileStruct']['fileStart'][0][0][0]
    return starts


def get_trial_lengths(ends):
    starts = ends['fileStruct']['fileLength'][0][0][0]
    return starts

def get_ROI_order(mat):
    try:
        order = mat['targetROI'][0]
    except:
        order = mat['holoDat']['targetROI'][0][0][:,0]
    return order-1

def get_stim_times(mat):
    return mat['holoDat']['firstTimes'][0,0][:,1]

def get_stim_tags(mat):
    try:
        t=mat['stimTags'][0]-1
    except:
        t=mat['holoDat']['the_tags'][0,0][0]
    return t

def get_strength_key(mat):
    return mat['strengthKey'][0]

def get_mouse_and_day(mat):
    return mat['mouse'][0], mat['day'][0]

def get_target_distance(mat):
    try:
        order = mat['targetDistance'][0]
    except:
        order = mat['holoDat']['targetDistance'][0][0][:,0]
        
    return order

def get_motion_data(mat):
    motion = mat['Motion']
    return motion

def get_cell_planes(mat, target_order=True):
    
    planes = mat['Plane'][:,0]
    if target_order:
        rois = get_ROI_order(mat)
        planes = planes[rois]
    return planes

def get_cell_stat(mat, cell_stat, target_order = False):
    rratios = mat['signals'][cell_stat][0][0]
    rratios = np.asarray([cell[0] for cell in rratios])
    if target_order:
        return rratios[get_ROI_order(mat)]
    else:
        return rratios

def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%y%m%d")
    d2 = datetime.strptime(d2, "%y%m%d")
    return abs((d2 - d1).days)

def get_index_of_best_targets(mat, dist_cutoff = 10):
    targets = get_ROI_order(mat)
    dist = get_target_distance(mat)
    best_targets = []
    for target in np.unique(targets):
        target_idx = np.where(targets==target)[0]
        if len(target_idx)>1:
            single_target = target_idx[dist[target_idx].argmin()]
        else:
            single_target = target_idx[0]
        if dist[single_target] < dist_cutoff:
            best_targets.append(single_target)
    return best_targets

def mean_of_mean_diffs(mean_post):
    #gives a single value per stim value, rather than per trial

    mean_diffs = mean_post.groupby(['cell','stim'])['differ'].mean()
    sem_diffs = mean_post.groupby(['cell','stim'])['differ'].sem()

    mean_diffs = mean_diffs.reset_index(level=['cell','stim'])
    mean_diffs['sem'] = sem_diffs.values
    mean_diffs_no_long = mean_diffs.loc[(mean_diffs.stim!=1000) & (mean_diffs.stim!=500),:]
    #mean_post_no_long = mean_post.loc[(mean_post.stim!=10000) & (mean_diffs.stim!=500),:]
    return mean_diffs_no_long

def get_normed_x_y_err(these_vals):
    x_mean=these_vals.stim.values
    try:
        y_mean=these_vals.differ.values
    except:
        y_mean = these_vals.df.values
    semy=these_vals['sem'].values

    y_mean = y_mean - y_mean.min()
    semy = semy/y_mean.max()
    y_mean = y_mean/y_mean.max()
    return x_mean,y_mean,semy
   
def get_sat_vals(mean_df, func, sat_changes=False, remove_downtrend=False, bounds=None):
    mean_df = mean_of_mean_diffs(mean_df)
    sat_vals=[]
    for cell in mean_df.cell.unique():
        these_vals = mean_df[mean_df.cell==cell]
        x_mean=these_vals.stim.values
        y_mean=these_vals.differ.values
        semy=these_vals['sem'].values
        y_mean = y_mean - y_mean.min()
        semy = semy/y_mean.max()
        y_mean = y_mean/y_mean.max()
        idx = np.argsort(x_mean)
        if y_mean[idx[0]] > y_mean[idx[len(y_mean)-1]]:
            #print("suppressed")
            sat_vals.append(np.nan)
            continue
        if remove_downtrend:
            for _ in range(2): #allow chopping of up to two end values if they're off
                if y_mean[-1]+semy[-1] < y_mean[-2]-semy[-2]:
                    x_mean = x_mean[0:-1]
                    y_mean = y_mean[0:-1]
        try:
            if bounds is not None:
                sig2opt = basic_fit(func,  x_mean,  y_mean, bounds=bounds)#, p0= estimate_p0_logistic(x_mean,y_mean)[0:2])
            else:
                sig2opt = basic_fit(func,  x_mean,  y_mean)
            #print('this one worked')
        except:
            sat_vals.append(np.nan)
            continue
            
        y_vals = func(np.linspace(0,2000,2000), *sig2opt)
        if sat_changes:
            sat_approach = np.where(y_vals>.85*sig2opt[0])
        else:
            sat_approach = np.where(y_vals>.85)
        if len(sat_approach[0])>0:
            sat_vals.append(sat_approach[0][0])
        else:
            sat_vals.append(np.nan)
    return sat_vals

