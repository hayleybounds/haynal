import pandas as pd
import scipy.io as sio
import numpy as np
from glob import glob
import pickle
from scipy.stats.mstats import zscore
from ScanImageTiffReader import ScanImageTiffReader
from PSTH_creation_utils import *
import os
from scipy.stats import f_oneway

from gen_utils import mean_by_id

def get_orientations_py(mat):
    all_grating_info = np.vstack(mat['result']['gratingInfo'][0,0][0,0].tolist()[0:5])
    return all_grating_info[0,:]

def get_orientations(ends):
    starts = ends['visData']['condInfo'][0][0][0,:]
    return starts

def get_locinds(mat):
    locinds = mat['visData']['locinds'][0][0]
    locinds = np.transpose(locinds)
    return locinds


def get_median_loc(mat):

    return mat['allCoM']

def traces_to_vis_df(traces, mat_path, exp_type,do_cell_variables = False, is_baseline=False):
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
    
    #define trial stim variables
    if exp_type == 'ori':
        tags = get_orientations(mat)
        if len(all_melted['minor'].unique()) > len(tags):
            raise ValueError('number of stim tags is less than n trials in psf array');
        for trial in all_melted['minor'].unique():
            all_melted.loc[all_melted.minor==trial,'orientation'] = tags[trial]
    elif exp_type == 'retino':
        tags = get_locinds(mat)
        x=tags[1,:]
        y=tags[0,:]
        for trial in all_melted['minor'].unique():
            all_melted.loc[all_melted.minor==trial,'x'] = x[trial]
            all_melted.loc[all_melted.minor==trial,'y'] = y[trial]
    if do_cell_variables:
        #define cell level variables
        medians = get_median_loc(mat)
           
        for cell in all_melted['major'].unique():
            all_melted.loc[all_melted.major==cell,'medx'] = medians[cell,0]
            all_melted.loc[all_melted.major==cell,'medy'] = medians[cell,1]
            
   
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})
    return all_melted

def get_trial_starts(ends):
    starts = ends['fileStruct']['fileStart'][0][0][0]
    return starts
def get_trial_lengths(ends):
    starts = ends['fileStruct']['fileLength'][0][0][0]
    return starts
def get_traces(traces):
    traces = traces['signals']['F'][0][0]
    return traces
def get_ROI_order(mat):
    roi_order = mat['holoDat']['targetROI'][0][0][:,0]
    return roi_order-1
def get_target_distance(mat):
    roi_order = mat['holoDat']['targetDistance'][0][0][:,0]
    return roi_order-1


def get_index_of_best_targets(mat_path, dist_cutoff = 10):
    mat=sio.loadmat(mat_path)
    targets = get_ROI_order(mat)
    dist = get_target_distance(mat)
    best_targets = []
    for target in np.unique(targets):
        target_idx = np.where(targets==target)[0]
        if len(target_idx)>1:
            single_target = target_idx[dist[target_idx].argmin()]
        else:
            single_target = target_idx[0]
        if dist[single_target] < 10:
            best_targets.append(single_target)
    return best_targets

#util functions for transforming loaded data
def get_trialwise_data(mat_path, dfof_method, traces=None,
                   df_range = None, do_zscore=False, do_baseline = False, baseline_n=None, skip_ends = 0,
                      skip_starts = 0, baseline_subt = False):
    """given a mat path that has the traces and filedata, chop into trialwise traces
    mat_path: path to mat containing FileStruct which has a value fileStart that delims
    starts of files.
    dfof_method: str, one of 'percentile', 'by trial', or 'at start'. prcentile takes the 30th 
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
    lengths = get_trial_lengths(mat).astype('int16')
    if traces is None:
        traces = get_traces(mat)
    
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, lengths, skip_ends, skip_starts)
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    
    #initialize psf array
    length = min(lengths)
    if length<3:
        print('WARNING UR SHORTEST IS TOOO SHORT')
        length = 15
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    for i in range(0,len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        cut_traces[i,:,:] = traces[:,start:start+length]
    
    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces


def get_ordered_traces(mat):
    """loads the mat file at mat_path, gets the traces
    and orders them by ROI order.
    """
    traces = get_traces(mat)
    order = get_ROI_order(mat)
    traces = traces[order[:,0],:]
    return traces

def get_holo_starts(mat):
    starts = mat['holoDat']['firstTimes'][0][0][:,0]
    starts =starts*float(mat['fileStruct']['FR'][0][0][0])
    return starts.astype('int16')

#util functions for transforming loaded data
def get_psf_sing_holo_seq(mat_path, pre_time, length, dfof_method,
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
    lengths = get_trial_lengths(mat).astype('int16')
    traces = get_ordered_traces(mat)
    holo_start_times = get_holo_starts(mat)
    
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, lengths, skip_ends, skip_starts)
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    
    #now order into ordered traces?
    
    #initialize psf array
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    
    for i in range(0,len(starts)):
        start = int(starts[i])-1 #convert for 0 indexing
        for cell in range(0, traces.shape[0]):
            real_start = start+holo_start_times[cell]-pre_time
            real_end = real_start+length
            cut = get_traces_segment(real_end, real_start, traces)
    
    #then dfof/baseline subtract
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces


        

def get_stim_psfs(mat_path, dfof_method, tifFolder, picklefile, pre_time, length, traces = None,
                   df_range = None, do_zscore=False, do_baseline = False, baseline_n=None, skip_ends = 0,
                      skip_starts = 0, baseline_subt = False):
    """given a mat path that has the traces, chop into post stim psfs. Uses the new scanimage stuff to actually
    get the exact stim frames and use that.
    mat_path: path to mat containing FileStruct which has a value fileStart that delims
    starts of files.
    dfof_method: str, one of 'percentile', 'by trial', or 'at start'. percentile takes the 30th 
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
    if traces is None:
        traces = get_traces(mat)
    starts, lengths = get_stim_starts(tifFolder, picklefile)
    print(starts, lengths)
    
    #if skipping anything, get rid of them
    starts, lengths = chop_file_vars(starts, lengths, skip_ends, skip_starts)
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    print(traces.shape)
    #initialize psf array
    cut_traces = np.zeros((int(np.shape(starts)[0]), int(np.shape(traces)[0]), length))
    for i in range(0,len(starts)):
        start = int(starts[i])-pre_time #its already conv for 0 indexing!
        cut_traces[i,:,:] = traces[:,start:start+length]
    
    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces

def get_continuous_stim_starts(tifFolder, picklefile, matfile=None):
    files=glob(picklefile)
    if len(files) > 0:
        print('found pickle file')
        with open(files[0],'rb') as f:
            [stim_starts] = pickle.load(f)
            return stim_starts
    else:
        files = glob(tifFolder+'*.tif')
        if len(files)==0:
            print('no tifs found in ', tifFolder)
            return
        else:
            print(len(files), ' tifs found in ', tifFolder)
        files.sort()
        stim_on_frame = [];
        FR=None
        
        counter = 0
        for file in files:
            with ScanImageTiffReader(file) as reader:
                
                stim_on = []

                for framei in range(0,2,reader.shape()[0]):
                    string = reader.description(framei).split('auxTrigger')[1].split('[')[1]
                    if string[0]!=']' and string[0]!='-':
                        stim_on.append(framei)
                if len(stim_on) > 0:
                    [stim_on_frame.append(so+counter) for so in stim_on]
                
                FR = float(reader.metadata().split('scanVolumeRate = ')[1].split('\n')[0])
                counter = counter+reader.shape()[0]
                
                #int(reader.metadata().split('hBeams.powers = ')[1].split('\n')[0])
      
        stim_starts = np.asarray(stim_on_frame)
        stim_starts = (stim_starts/(2*3)).astype('uint16')
     
        dat = [ stim_starts]
        with open(picklefile, 'wb') as f:
            pickle.dump(dat, f)
        print('saved data!')
        if matfile is not None:
            saveme = {'fileStart':file_starts+1, 'fileLength':file_lengths,'stimStart':stim_starts, 'FR':FR}
            sio.savemat(matfile,saveme)
        return stim_starts

    
"""NEW PYTHON SUITE2P FUNCS"""
def get_trialwise_data_pydat(base_path, dfof_method, tifFolder, pre_time, epoch, get_spikes = False,
                                     stim_aligned=False, 
                                     traces=None,
                                     df_range = None, do_zscore=False, do_baseline = False, baseline_n=None):
    """load F traces and get trialwise data from python suite2p outputs
    Determines length of framedata based on the shortest trial
    
    Inputs:
        base_path (string): root folder of saved suite2p data
        dfof_method (string): one of 'percentile', 'by trial', or 'at start'. percentile takes the 30th 
            prctile value per cell as the f0. 'rolling_percentile' is the same but calculates rolling f0 based on 200 datapts.
            by trial does it so that the f0 is the mean of the values in df_range per trial, at_start takes the mean of the
            frames in df_range at the start of the recording
        tifFolder (string): path to folder containing tif files for finding file lengths
        pre_time (int): number of frames before trial start to include
        epoch (int): idx of epoch to analyze
        
            
    Kword Inputs:
        df_range (int): used by by trial + at start dfof methods, mean diff things for each of them
        do_zscore (bool): if used, before traces are cut, they are zscored by cell. If dfof_method is percentile, traces
            are first dff'ed
        do_baseline: option to subtract the mean of baseline_n frames from the psf. This is per trial. Done after dff. 
        baseline_n: how many frames, beginning at first frame of each trial (as defined by psf, not from beginning
        of acq, to use for baseline subtraction
        
    """
    
    if traces is None:
        _,_,_, sp, traces = get_iscell_F_Neu(base_path)
       
    if get_spikes:
        traces = sp
        
    file_starts, full_lengths, stim_starts = get_starts_multifolder(tifFolder, base_path, get_stim_starts=stim_aligned)
    if stim_aligned:
        file_starts = stim_starts
        
    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    length = min(full_lengths[epoch-1])
    
    #initialize psf array - it has to only be for one epoch bc number of trials varies
    cut_traces = np.zeros((len(file_starts[epoch-1]), np.shape(traces)[0], length))
    
    fstarts = file_starts[epoch-1]
    epoch_lengths = [sum(flen) for flen in full_lengths]
    fstarts = fstarts + sum(epoch_lengths[0:epoch-1])
    
    #for each trial get the file start, and the relative stim start, and add it to traces
    for istim, stim_start in enumerate(fstarts):
        real_start = stim_start-pre_time
        real_end = real_start+length
        cut = get_traces_segment(real_start, real_end, traces)
        cut_traces[istim,:,:] = cut
    
    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces




def traces_to_vis_df_py(traces, mat_path, exp_type, contrasts=None, base_path=None, drop_locs=[]):
    """from the output of trialwise data, create a pandas dataframe w/various
    info.
    traces: output of trialwise data
    mat_path: path to result mat
    base_path: if provided, will use to extract cell info
    drop_locs: for recordings with corrupted tifs, drop the trials corresponding to those tifs
    """
    
    #change traces into a pandas array with cell and time as cols
    print('shape of traces', traces.shape)
    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))
    
    #define trial stim variables
    if exp_type == 'ori':
        mat = sio.loadmat(mat_path)
        tags = get_orientations_py(mat)
        if len(all_melted['minor'].unique()) > len(tags):
            raise ValueError('number of stim tags is less than n trials in psf array');
        for trial in all_melted['minor'].unique():
            all_melted.loc[all_melted.minor==trial,'orientation'] = tags[trial]
    if exp_type == 'noise':
        tags = contrasts
        if len(all_melted['minor'].unique()) > len(tags):
            print(len(all_melted['minor'].unique()), len(tags))
            raise ValueError('number of stim tags is less than n trials in psf array');
        for trial in all_melted['minor'].unique():
            all_melted.loc[all_melted.minor==trial,'contrast'] = tags[trial]
    elif exp_type == 'retino':
        mat = sio.loadmat(mat_path)
        tags = mat['result']['locinds'][0][0].T
        x=tags[1,:]
        y=tags[0,:]
        if len(drop_locs)>0:
            x= np.delete(x,np.array(drop_locs))
            y = np.delete(y, np.array(drop_locs))
        for trial in all_melted['minor'].unique():
            all_melted.loc[all_melted.minor==trial,'x'] = x[trial]
            all_melted.loc[all_melted.minor==trial,'y'] = y[trial]
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})
    if base_path is not None:
        all_melted = add_cell_info(all_melted, base_path)
  
   
    
    return all_melted

        
def get_mean_maps(ret_df, start, stop):
    ret_df.loc[:,'abs_id'] = ret_df.x+ret_df.y*30
    mean_by_loc = mean_by_id(ret_df,start,stop,col='abs_id', add_columns=['x','y'])
    mean_maps = np.empty((len(mean_by_loc.x.unique()), len(mean_by_loc.y.unique()), len(mean_by_loc.cell.unique())))
    mean_by_loc.x  = mean_by_loc.x.astype('int32')
    mean_by_loc.y = mean_by_loc.y.astype('int32')
    for cell in mean_by_loc.cell.unique():
        rel_vals = mean_by_loc[mean_by_loc.cell==cell]
        for x in mean_by_loc.x.unique():
            for y in mean_by_loc.y.unique():
                mean_maps[x-1,y-1,cell] = rel_vals[(rel_vals.x==x)&(rel_vals.y==y)].df
                
    return mean_maps

def get_p_retino(mean_df, col):
    ps=[]
    for cell in mean_df.cell.unique():
        #use list comprehension to get a list of all the values for this cell for each stim point
        vals = mean_df[mean_df.cell==cell] #vals for this cell
        stims = vals[col].unique()
        [_, p] = f_oneway(*[vals.df[vals[col]==idu] for idu in stims])
        ps.append(p)
    return np.asarray(ps)


def get_online_mean_maps(oret_path):
    oret_mat=sio.loadmat(oret_path)
    PSTHs = oret_mat['PSTHs']
    locs = oret_mat['locs']
    PSTHs = PSTHs - np.nanmean(PSTHs[:,:,0:5],2).reshape(PSTHs.shape[0],PSTHs.shape[1],1)
    online_mean_maps = np.zeros((len(np.unique(locs[:,0])), len(np.unique(locs[:,1])), PSTHs.shape[1]))
    for x in np.unique(locs[:,0]):
        for y in np.unique(locs[:,1]):
            these_trials = np.where((locs[:,0]==x) & (locs[:,1]==y))[0]
            online_mean_maps[x-1,y-1,:] = np.nanmean(np.nanmean(PSTHs[these_trials,:,5:],0),1)
    return online_mean_maps
 