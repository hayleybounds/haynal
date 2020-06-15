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
import h5py
import os
from PSTH_creation_utils import *
from scipy.stats import f_oneway
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import r2_score
from datetime import datetime
from scipy.optimize import curve_fit
from tiff_utils import *
from gen_utils import *
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift
import time

## These are for .mat files <7.1 (ie. not hdf5)

def get_holo_targs_by_seq(mat, is_pbal = False):
    """used by other functions, goes thru all conditions and gets first roi and finds targs it has
    gives you the targets in each sequence (as in big list of sequences) and the hologram roi name"""
    conds = get_trial_conds(mat)
    first_roi = np.array([s[0][0] if len(s)>0 else 0
                          for s in mat['ExpStruct'][0,0]['outParams'][0,0]['sequence'][0]])
    #hack, for now the ensemble size is always the same so use a numpy array
    rois = mat['ExpStruct'][0,0]['holoRequest'][0,0]['rois'][0]
    print('roi length is ', rois.shape)
    rois = np.squeeze(np.dstack(rois))
    ens_size = rois.shape[0]
    targs_in_seq = np.zeros((len(first_roi),ens_size))
    for iroi,roi in enumerate(first_roi):
        if roi==0:
            targs_in_seq[iroi,:] = np.nan
        else:
            #two different -1 for matlab fixes
            targs_in_seq[iroi,:] = rois[:,roi-1]-1 #MATLAB INDEX FIX
    if is_pbal:
        #print(targs_in_seq)
        max_per_seq = targs_in_seq.max(axis=1)
        #print(max_per_seq)
        is_bal = max_per_seq<=int(get_holo_targs(mat).shape[0]/2)
        #print(is_bal)
        targs_in_seq[targs_in_seq>int(get_holo_targs(mat).shape[0]/2)] -= int(get_holo_targs(mat).shape[0]/2)
    else:
        is_bal = np.zeros((targs_in_seq.shape[0],1))
    return targs_in_seq, first_roi, is_bal

def get_trial_stimmed_targs(mat, is_pbal=False):
    """using get_holo_targs_by_seq gets the mapping of targets to condition. Then converts them to trialwise"""
    conds = get_trial_conds(mat)
    targs_in_seq, first_roi, is_bal = get_holo_targs_by_seq(mat, is_pbal=is_pbal)
    targs_per_trial = targs_in_seq[conds-1,:]
    hologram_per_trial = first_roi[conds-1]
    bal_per_trial = is_bal[conds-1]
    
        
    return targs_per_trial, hologram_per_trial, bal_per_trial #hologram per trial is basically sequence

def get_targs_per_holo(mat):
    """for each hologram trial type, get the holo targs in that hologram"""
    targs_in_seq, first_roi, is_bal = get_holo_targs_by_seq(mat)
    #make a dict matching targs to holograms
    targs_per_holo = {}
    for iroi,roi in enumerate(first_roi):
        targs_per_holo[roi]=targs_in_seq[iroi]
    return targs_per_holo

def get_first_stim_times(mat, get_all = False):
    stim_times = mat['ExpStruct'][0,0]['outParams'][0,0]['firstStimTimes'][:,1]
    
    if get_all:
         stim_times = mat['ExpStruct'][0,0]['outParams'][0,0]['firstStimTimes']
    return stim_times

def get_trial_conds(mat):
    trial_conds = mat['ExpStruct'][0,0]['trialCond'][0]
    return trial_conds

def get_powers(mat):
    pwr_list = mat['ExpStruct'][0,0]['outParams'][0,0]['power'][0]
    return pwr_list

def get_pulses(mat):
    pulse_list = mat['ExpStruct'][0,0]['outParams'][0,0]['pulses'][0]
    return pulse_list

def get_ncells(mat):
    cell_list = mat['ExpStruct'][0,0]['outParams'][0,0]['cellsPerHolo'][0]
    return cell_list

def get_trial_power(mat):
    """Returns a list of the power used (in mW) on a per trial basis.

    Inputs:
        mat (dict): data from the NoDaq .mat loaded by scipy.io.loadmat
    """

    trial_power = []
    conds = get_trial_conds(mat)
    pwr_list = get_powers(mat)

    for i,cond in enumerate(conds):
        trial_power.append(pwr_list[int(cond-1)]) # stupid matlab indexing fix

    return trial_power

def get_trial_pulses(mat):
    """returns a list of the number of pulses on a per trial basis
    Inputs:
        mat (dict): data from the NoDaq .mat loaded by scipy.io.loadmat
    """
    trial_power = []
    conds = get_trial_conds(mat)
    pwr_list = get_pulses(mat)

    for i,cond in enumerate(conds):
        trial_power.append(pwr_list[int(cond-1)]) # stupid matlab indexing fix

    return trial_power

def get_trial_ncells(mat):
    """returns a list of the number of ncells on a per trial basis
    Inputs:
        mat (dict): data from the NoDaq .mat loaded by scipy.io.loadmat
    """
    trial_power = []
    conds = get_trial_conds(mat)
    pwr_list = get_pulses(mat)

    for i,cond in enumerate(conds):
        trial_power.append(pwr_list[int(cond-1)]) # stupid matlab indexing fix

    return trial_power

def get_holo_targs(mat):
    holo_targs = mat['ExpStruct'][0,0]['holoRequest'][0,0]['targets']
    return holo_targs

def get_holo_weights(mat):
    return mat['ExpStruct'][0,0]['holoRequest'][0,0]['roiWeights']


""" Matching Holos to Suite2p Fxns """
def get_masks_shift(base_path, masks_path, ran):
    pass

def get_target_mapping(mat, base_path, planes_key=None, pixel_shift = None):
    """Maps the holo targs to suite2p extracted data. Returns a numpy array.

    Inputs:
        mat (dict): data from the NoDaq .mat loaded by scipy.io.loadmat
        base_path (str): path to suite2p folder where 'combined' folder exists
    Optional:
        planes_key (list or arr): 1d arr of all the planes, for if some aren't shot
        pixel_shift (arr, planesx2): shift btwn masks and suite2p.
    """
    mat = mat.copy() #in place stuff was fucking my shit UP
    holo_targs = get_holo_targs(mat)
    medians = get_cell_field(base_path, 'med')
    planes = get_cell_field(base_path, 'iplane')
    if planes_key is None:
        planes_key = np.unique(holo_targs[:,2])
        planes_key.sort()
    target_mapping = np.empty((holo_targs.shape[0],),dtype='int16')
    target_dists = np.empty((holo_targs.shape[0],),dtype='float64')

    for i,targ in enumerate(holo_targs):
        plane_idx = np.where(planes_key==targ[2])[0]
        cells_this_plane = np.where(planes==plane_idx)[0]
        meds_this_plane = medians[cells_this_plane,:]
        
        if len(cells_this_plane)==0:
            print('WARNING THERES NOTHING ON PLANE ', plane_idx)
            target_mapping[i] = 0
            target_dists[i] = 400
            continue

        if plane_idx == 1:
            meds_this_plane = meds_this_plane - np.array([0, 512])
        elif plane_idx == 2:
            meds_this_plane = meds_this_plane - np.array([512,0])
        elif plane_idx == 3:
            meds_this_plane = meds_this_plane - np.array([512,512])

        targ = targ[:2]
        if pixel_shift is not None:
            targ += np.squeeze(pixel_shift[plane_idx,:])
        distances = np.linalg.norm(meds_this_plane-targ, axis=1)
        
        
        idx_best_match = np.argmin(distances)

        target_mapping[i] = cells_this_plane[idx_best_match]
        target_dists[i] = np.min(distances)

    return target_mapping, target_dists


def get_index_of_best_targets(mat_path, base_path, dist_cutoff = 10, planes_key=None,pixel_shift=None):
    mat = sio.loadmat(mat_path)
    #here, targets is a list of the suite2p rois, where the index corresponds the target
    targets, dist = get_target_mapping(mat, base_path, planes_key=planes_key, pixel_shift=pixel_shift)
    best_targets = []
    #take all uniquely appearing suite2p rois
    for target in np.unique(targets):
        #find the target indexes where this suite2p roi was matched as closest cell
        target_idx = np.where(targets==target)[0]
        if len(target_idx)>1:
            #if there are multiple times where this suite2p roi was matched
            #to a target, take only the target index that the closest to this
            #suite2p roi
            single_target = target_idx[dist[target_idx].argmin()]
        else:
            single_target = target_idx[0]
        #if the distance of the suite2p roi to the target of single_idx is 
        #below threshold, it is good
        if dist[single_target] < dist_cutoff:
            best_targets.append(single_target)
    return best_targets

############util functions for transforming loaded data##############33333
def py_get_psth(base_path, mat_path, dfof_method, tifFolder, pre_time, length, epoch, is_online_data = False,
                              get_spikes = False, traces=None,
                              df_range = None, do_zscore=False, do_baseline = False, baseline_n=None,
                first_stim_times=None, stim_aligned=False, planes_key=None):
    """Computes the peristimulus fluorescence values as n_trials x n_cells x length.
    Loads F traces and get trialwise data from python suite2p outputs.
    Determines length of framedata based on the shortest trial.

    Inputs:
        base_path (str): path to suite2p folder where 'combined' folder exists
        dfof_method (str): one of 'percentile', 'by trial', or 'at start'. percentile takes the 30th
            prctile value per cell as the f0. 'rolling_percentile' is the same but calculates rolling f0 based on 200 datapts.
            by trial does it so that the f0 is the mean of the values in df_range per trial, at_start takes the mean of the
            frames in df_range at the start of the recording
        tifFolder (str): path to folder containing tif files for finding file lengths
        pre_time (int): number of frames before trial start to include
        length (int): number of frames to analyze (includes pre-time)
        epoch (int): idx of epoch to analyze


    Kword Inputs:
        df_range (int): used by by trial + at start dfof methods, mean diff things for each of them
        do_zscore (bool): if used, before traces are cut, they are zscored by cell. If dfof_method is percentile, traces
            are first dff'ed
        do_baseline: option to subtract the mean of baseline_n frames from the psf. This is per trial. Done after dff.
        baseline_n: how many frames, beginning at first frame of each trial (as defined by psf, not from beginning
        of acq, to use for baseline subtraction
        first_stim_times: if don't want to auto retrieve stim times of cond 1 (so order 2) can give them here

    """

    # Load mat file
    mat = sio.loadmat(mat_path)

    # Get stim times and target mapping
    if first_stim_times is None:
        first_stim_times = get_first_stim_times(mat)

    fstarts, full_lengths, stim_starts = get_starts_suite2p(tifFolder, base_path, epoch, get_stim_starts=stim_aligned)
    if stim_aligned:
        fstarts = stim_starts
    
    if not is_online_data and epoch>1:
        fstarts = fstarts + get_prev_folder_length(base_path, epoch)
        print('len fstarts', len(fstarts))
        
    #file_starts, full_lengths, _ = get_starts_multifolder(tifFolder, base_path, get_stim_starts=False)
    FR = get_fr(tifFolder + base_path.split('/')[-3].split('_')[0] +'/') 
    first_stim_times *= FR
    
    if traces is None:
        _,_,_, sp, traces = get_iscell_F_Neu(base_path)
    
    if get_spikes:
        traces = sp

    if not is_online_data:
        holo_target_matches, dists = get_target_mapping(mat, base_path, planes_key=planes_key)
        traces = traces[holo_target_matches,:]
    
    if traces.shape[0]>first_stim_times.shape[0]:
        print('padding first stim times length')
        anStimTimes=np.zeros((1,traces.shape[0]))
        anStimTimes[0,0:first_stim_times.shape[0]] = first_stim_times
        first_stim_times=anStimTimes[0,:]

    #fstarts = file_starts[epoch-1]
    #epoch_lengths = [sum(flen) for flen in full_lengths]

    #if not is_online_data:
    #    fstarts = fstarts + sum(epoch_lengths[0:epoch-1])

    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    #initialize psf array
    cut_traces = np.zeros((int(np.shape(fstarts)[0]), int(np.shape(traces)[0]), length))

    for i in range(len(fstarts)):
        start = int(fstarts[i])-1 #convert for 0 indexing
        for cell in range(traces.shape[0]):
            real_start = start+first_stim_times[cell].astype('int16')-pre_time
            real_end = real_start+length
            cut_traces[i,cell,:] = get_traces_segment(real_start, real_end, traces,cell=[cell])

    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces



def nodaq_stimtest_df(traces, mat_path, do_cell_variables = False, base_path = None, planes_key = None):
    """from the output of get_psf_sing_holo_seq, create a pandas dataframe w/various
    info.
    traces: output of get_psf_sing_holo_seq
    mat_path: path to stim info, has fields 'stimTags' and 'strengthKey' which is ordered list of
        what strength numbers in 'stimTags' corresponds to.
    do_cell_variables: whether to add various cell variables like location
    base_path: the base suite2p path to get all the good variables from
    """

    # change traces into a xarray for lts of arrays
    # xr.DataArray(data, [items, major, minor])
    # array = xr.DataArray(np.transpose(traces)).to_dataset(dim='dim_1')
    # array = array.reset_index(level=['major','minor'])
    # all_melted = pd.melt(array, ('major','minor'))

    # legacy code
    #change traces into a pandas array with cell and time as cols

    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))

    # for (no)daq stim Tests, first, get the power used on a trial-by-trial basis
    mat = sio.loadmat(mat_path)
    trial_stim_cond = get_trial_power(mat)
    trial_pulse_cond = get_trial_pulses(mat)
    trial_cell_cond = get_trial_ncells(mat)

    if len(all_melted['minor'].unique()) > len(trial_stim_cond):
        raise ValueError('Number of trial conditions is less than n trials in array.')
    elif len(all_melted['minor'].unique()) < len(trial_stim_cond):
        print('WARNING TRIALS IN DF SMALlER THAN CONDS, CUTTING CONDS')
        trial_stim_cond = trial_stim_cond[:len(all_melted['minor'].unique())]
        trial_pulse_cond = trial_pulse_cond[:len(all_melted['minor'].unique())]
        trial_cell_cond = trial_cell_cond[:len(all_melted['minor'].unique())]
    trial_df = pd.DataFrame({'minor':list(range(len(trial_stim_cond))),
                                         'power': trial_stim_cond, 'pulses': trial_pulse_cond,
                            'ncells_holo': trial_cell_cond})
    all_melted = all_melted.merge(trial_df, on='minor')

    #for trial in all_melted['minor'].unique():
    #    all_melted.loc[all_melted.minor==trial,'power'] = trial_stim_cond[trial]
    #    all_melted.loc[all_melted.minor==trial,'pulses'] = trial_pulse_cond[trial]
        
    
    
    if do_cell_variables and base_path is not None:
        mapping, dists = get_target_mapping(mat, base_path, planes_key = planes_key)
        medians = get_cell_field(base_path, 'med')
        planes = get_cell_field(base_path, 'iplane')
        for cell in all_melted['major'].unique():
            all_melted.loc[all_melted.major==cell, 'medianx'] = medians[mapping[cell]][0]
            all_melted.loc[all_melted.major==cell, 'mediany'] = medians[mapping[cell]][1]
            all_melted.loc[all_melted.major==cell, 'cell_idx'] = mapping[cell]
            all_melted.loc[all_melted.major==cell, 'dist'] = dists[cell]
            all_melted.loc[all_melted.major==cell, 'plane'] = planes[mapping[cell]]
            

    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})

    return all_melted

def multi_ensemble_stim_df(traces, mat_path, cut_motion=None, planes_key=None,
                          base_path=None, pixel_shift=None,is_pbal=False):
    """from the output of get_psf_sing_holo_seq, create a pandas dataframe w/various
    info.
    traces: output of get_psf_sing_holo_seq
    mat_path: path to stim info, has fields 'stimTags' and 'strengthKey' which is ordered list of
        what strength numbers in 'stimTags' corresponds to.
    do_cell_vars: whether to add various cell variables like redness to df
    """
    s=time.time()
    mat = sio.loadmat(mat_path)
    powers = get_trial_power(mat)
    targs_by_trial, holo_by_trial, bal_per_trial = get_trial_stimmed_targs(mat,is_pbal=is_pbal)
    mapping, dist = get_target_mapping(mat, base_path, planes_key,pixel_shift=pixel_shift)
    
    traces = traces.astype('float16')
    #change traces into a pandas array with cell and time as cols
    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))
    print('first stuff', time.time()-s)
    s=time.time()
    if is_pbal:
        trial_df = pd.DataFrame({'minor':list(range(len(powers))),'power':powers,'holo':holo_by_trial,
                                'is_bal':bal_per_trial})
    else:
        trial_df = pd.DataFrame({'minor':list(range(len(powers))),'power':powers,'holo':holo_by_trial})
    all_melted = all_melted.merge(trial_df, on='minor')
    print('first melt', time.time()-s)
    s=time.time()
    
    #NEW 5/18 - ONLY DISALLOW CELLS MAPPED MULTIPLE TIMES
    single_map_targ_idx = get_index_of_best_targets(mat_path, base_path, np.inf,
                                          planes_key = planes_key, pixel_shift=pixel_shift)
    print('THERE ARE', len(np.where(dist[np.setdiff1d(list(range(len(mapping))),single_map_targ_idx)]<5)[0]), 'HOLOS EXCLUDED CLOSER THAN LIM')
    all_holos_targs = get_targs_per_holo(mat)
    for _,h in all_holos_targs.items():
        targs_single_match_trial = np.intersect1d(single_map_targ_idx, h)
        print('shot trial',len(h),'matched trial',len(targs_single_match_trial))
    
    
    for trial in all_melted['minor'].unique():
        this_set = all_melted.minor==trial
        
        targs_for_this_trial = targs_by_trial[trial]
        targs_single_match_trial = np.intersect1d(single_map_targ_idx, targs_for_this_trial)
        
        cells_this_holo = mapping[targs_single_match_trial.astype('uint16')]
        all_melted.loc[(this_set) & (all_melted.major.isin(cells_this_holo)),'stimmed'] = 1
        #all_melted.loc[(this_set) and not (all_melted.major.isin(cells_this_holo)),'stimmed'] = False
        if cut_motion is not None:
            #for time reasons, try trialwise
            s=time.time()
            for iplane in range(len(np.unique(planes))):
                all_melted.loc[all_melted.major.isin(cells_per_plane[iplane])
                               & (this_set),'max_motion'] = np.max(cut_motion[trial,iplane,:])
            print(time.time()-s)
            #if seq[trial] is not 0:                    
            #    for cell in all_melted[(this_set) & all_melted.stimmed==1].major.unique():
            #        this_set_cell = all_melted.major==cell
            #        all_melted.loc[(this_set) & (this_set_cell), 'max_motion'] = np.max(cut_motion[trial,cell,:])
            #        all_melted.loc[(this_set) & (this_set_cell), 'mean_motion'] = np.nanmean(cut_motion[trial,cell,:])
            #else:
            #    for cell in all_melted.major.unique():
            #        this_set_cell = all_melted.major==cell
            #        all_melted.loc[(this_set) & (this_set_cell), 'max_motion'] = np.max(cut_motion[trial,cell,:])
            #        all_melted.loc[(this_set) & (this_set_cell), 'mean_motion'] = np.nanmean(cut_motion[trial,cell,:])
    print('trial_info', time.time()-s)
    s=time.time()
        
        
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})  
    
    #add the hologram target info
    #cell_df = pd.DataFrame({'cell': mapping, 'target_numb': range(len(mapping)), 'target_dist': _})
    #all_melted = all_melted.merge(cell_df,on='cell',how='outer') #have to use outer to make it not delted cells no tin cell_df
    
    #add the plane and median
    if base_path is not None:
        all_melted = add_cell_info(all_melted, base_path, fix_offset=True)
    
        df = all_melted
        df.loc[:,'x']=df.medianx.values
        df.loc[:,'y']=df.mediany.values
    print(time.time()-s)
    return df
### PANDAS TRANSFORM FXNS ######

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


def get_normed_x_y_err(these_vals,val_col = 'df',stim_col='stim'):
    x_mean=these_vals[stim_col].values
    y_mean=these_vals[val_col].values
    semy=these_vals['sem'].values

    y_mean = y_mean - y_mean.min()
    semy = semy/y_mean.max()
    y_mean = y_mean/y_mean.max()
    return x_mean,y_mean,semy

def mean_of_mean_diffs(mean_post):
    #gives a single value per stim value, rather than per trial

    mean_diffs = mean_post.groupby(['cell','stim']).mean()
    sem_diffs = mean_post.groupby(['cell','stim']).df.sem()

    mean_diffs = mean_diffs.reset_index(level=['cell','stim'])
    mean_diffs['sem'] = sem_diffs.values
    mean_diffs_no_long = mean_diffs.loc[(mean_diffs.stim!=1000) & (mean_diffs.stim!=500),:]
    #mean_post_no_long = mean_post.loc[(mean_post.stim!=10000) & (mean_diffs.stim!=500),:]
    return mean_diffs_no_long

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
