
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

from datetime import datetime



def traces_to_behavior_df(traces, mat_path, get_target_rois=False):
    """from the output of get_psf_sing_holo_seq, create a pandas dataframe w/various
    info. This is for mag behavior instead of for chrome tests. Similar, but a little different
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
    
    #define trial stim variables
    stims = get_stim_data(mat)
    manipulated = get_manipulated(mat)
    beh_outcomes = get_beh_outcomes(mat)
    licked = get_licked(mat)
    if len(all_melted['minor'].unique()) > len(stims):
        raise ValueError('number of stims is less than n trials in psf array');
    for trial in all_melted['minor'].unique():
        all_melted.loc[all_melted.minor==trial,'stim'] = stims[trial]
        all_melted.loc[all_melted.minor==trial,'beh_out'] = beh_outcomes[trial]
        all_melted.loc[all_melted.minor==trial,'licked'] = licked[trial]
        all_melted.loc[all_melted.minor==trial,'holo_stim'] = manipulated[trial]
    
    
    if get_target_rois:
        targeted_targets = get_targeted_targets(mat)
        targeted_dist = get_targeted_dist(mat)
        all_melted.loc[all_melted.major.isin(targeted_targets), 'holo_target'] = True
        for i, cell in enumerate(targeted_targets):
            all_melted.loc[all_melted.major==cell,'dist'] = targeted_dist[i]
    #define overall variables
    #mouse, day = get_mouse_and_day(mat)
    #all_melted.loc[:,'mouse'] = mouse
    #all_melted.loc[:,'day'] = day
    #all_melted.loc[:,'file'] = mat_path
    #mice_inj_dict = {'I119.3':'180802', 'I119.4':'180802', 'HB17.2':'180808', 'HB17.3':'180809',
    #                 'HB23.1':'180814'}
    #dpi = days_between(mice_inj_dict[mouse], day)
    #all_melted.loc[:,'dpi'] = dpi
        
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})
    return all_melted

def get_targeted_targets(mat):
    return mat['targetedTargetROIs'][0]-1      

def get_targeted_dist(mat):
    return mat['targetedTargetDistance'][0]       
     
def get_stim_data(mat):
    return mat['stimulusData'][0]

def get_beh_outcomes(mat):
    return mat['behaviorOutcomes'][0]

def get_licked(mat):
    beh_out = mat['behaviorOutcomes'][0]
    beh_out[beh_out==4] = 1
    beh_out[beh_out==3] = 0
    beh_out[beh_out==2] = 0
    return beh_out

def get_manipulated(mat):
    return mat['manipulationLog'][0]

