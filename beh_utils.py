import pandas as pd
import numpy as np
import h5py
import scipy.io as sio
import matplotlib.pyplot as plt
from gen_utils import *
from scratch_vis_utils import traces_to_vis_df_py
import seaborn as sns
import time

def load_trial_df(tri_df_path, daq_path = None):
    #get pandas df of behavior stims and responses
    beh_trials_df = pd.read_csv(tri_df_path,
                                 sep='\t')
    #psychopy prefills a bunch of trials, but only these are actually run
    beh_trials_df = beh_trials_df[beh_trials_df.ran==1]
    #force handling of nan values
    beh_trials_df.Response = np.genfromtxt(beh_trials_df.Response)
    beh_trials_df.RespTime = np.genfromtxt(beh_trials_df.RespTime)
    if daq_path is not None:
        dataset = get_daq_inputs(daq_path)
        print('dataset shape ', dataset.shape, 'trials shape', len(beh_trials_df))
        if dataset.shape[0] >= len(beh_trials_df):
            dataset=dataset[:len(beh_trials_df),:,:]

            beh_trials_df.loc[:,'post_stim_lick'] = get_first_post_stim_lick(dataset)/20000
            beh_trials_df.loc[:,'pre_stim_lick'] = get_pre_stim_lick(dataset)/20000
            beh_trials_df.loc[:,'reward_time'] = get_reward_time(dataset)/20000
            starts, ends = get_stim_starts_ends(dataset)
            beh_trials_df.loc[:,'start'] = starts/20000
            beh_trials_df.loc[:,'end'] = ends/20000
    return beh_trials_df


def traces_to_beh_df(traces, trials_df, base_path=None):
    """from the output of get_trialwise_data or get_trialwise_data_pydat, and the behavioral output file,
    create pandas dataframe with trial info included.
    """
    all_melted = cut_to_melted(traces)    

    #define trial stim variables
    intensity = trials_df.intensity
    response = trials_df.Response
    response_time = trials_df.RespTime
    trial_vars = np.vstack([intensity.values, response.values, response_time.values]).T
    vectors = np.zeros((len(all_melted), trial_vars.shape[1]))
    for trial in all_melted['trial'].unique():
        vectors[all_melted.trial==trial,:] = trial_vars[trial, :]

    all_melted.loc[:,'intensity'] = vectors[:,0]
    all_melted.loc[:,'response'] = vectors[:,1]
    all_melted.loc[:,'response_time'] = vectors[:,2]
    
    if base_path is not None:
        all_melted = add_cell_info(all_melted, base_path)
    return all_melted


def plot_performance(tri_df_path):
    beh_trials_df = load_trial_df(tri_df_path, daq_path = None)
    performance = beh_trials_df.groupby(['intensity']).Response.value_counts(normalize=True)
    performance = performance.xs(1,level='Response')
    performance = performance.reset_index(level=['intensity'])
    plt.plot(performance.intensity+.5,performance.Response, 'ro-', label='lick rate')


"""%%%%%%%%%%%%%%%%%daq utils%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
def print_daq_notes(daq_path):
    try:
        print(sio.loadmat(daq_path)['ExpStruct']['notes'][0][0][0])
    except:
        print('you didnt make this method')

def get_daq_inputs(daq_path):
    dataset = []
    try:
        with h5py.File(daq_path,'r') as f:
            mat = f['ExpStruct']['inputs']
            for i in range(mat.shape[0]):
                dataset.append(np.array(f[mat[i][0]]))
    except:
        dataset = sio.loadmat(daq_path)['ExpStruct']['inputs'][0][0]
        ds = np.zeros((dataset.shape[1],dataset[0,0].shape[1], dataset[0,0].shape[0]))
        for i in range(dataset.shape[1]):
            ds[i,:,:] = dataset[0,i].T

        dataset = ds
    dataset = np.array(dataset)
    return dataset

def get_stim_starts_ends(dataset, stim_line=1):
    stim_info = np.squeeze(dataset[:,stim_line,:])

    start_stims = np.zeros((dataset.shape[0]))
    end_stims = np.zeros((dataset.shape[0]))
    start_stims[:] = np.nan
    end_stims[:] = np.nan
    start_stims[np.where(np.diff(stim_info,axis=1)>0)[0]] = np.where(np.diff(stim_info,axis=1)>0)[1]
    end_stims[np.where(np.diff(stim_info,axis=1)>0)[0]] = np.where(np.diff(stim_info,axis=1)<0)[1]
    return start_stims, end_stims


def get_mean_stim_start_end(dataset, stim_line=1):
    start_stims, end_stims = get_stim_starts_ends(dataset, stim_line=1)

    return np.nanmean(start_stims), np.nanmean(end_stims)

def get_first_post_stim_lick(dataset, stim_line=1, lick_line=4):
    mean_stim_start, mean_stim_end = get_mean_stim_start_end(dataset, stim_line=stim_line)
    lick_info = np.squeeze(dataset.copy()[:,lick_line,:])
    # set all pts before stim starts to 0 to get post stim_start licks
    lick_info[:,:int(mean_stim_start)] = 0
    first_post_stim_lick = np.zeros((dataset.shape[0]))
    first_post_stim_lick[:] = np.nan
    for tri_n in range(dataset.shape[0]):
        where_lick = np.where(np.diff(lick_info[tri_n,:])>0)[0]
        if len(where_lick)>0:
            first_post_stim_lick[tri_n] = where_lick[0]
    return first_post_stim_lick

def get_pre_stim_lick(dataset, stim_line=1, lick_line=4):
    mean_stim_start, mean_stim_end = get_mean_stim_start_end(dataset, stim_line=stim_line)
    lick_info = np.squeeze(dataset.copy()[:,lick_line,:])
    # set all pts before stim starts to 0 to get post stim_start licks
    lick_info[:,int(mean_stim_start):] = 0
    pre_stim_lick = np.zeros((dataset.shape[0]))
    pre_stim_lick[:] = np.nan
    for tri_n in range(dataset.shape[0]):
        where_lick = np.where(np.diff(lick_info[tri_n,:])>0)[0]
        if len(where_lick)>0:
            pre_stim_lick[tri_n] = where_lick[0]
    return pre_stim_lick

def get_reward_time(dataset, water_line=5):
    water_info = np.squeeze(dataset[:,water_line,:])

    water_time = np.zeros((dataset.shape[0]))
    water_time[:] = np.nan
    water_time[np.where(np.diff(water_info,axis=1)>0)[0]] = np.where(np.diff(water_info,axis=1)>0)[1]
    return water_time







"""%%%%%%%%%%%%%%%%%%%%%%%%%%%% gather data utils %%%%%%%%%%%%%%%%%%%%% """

def make_ret_df(base_path, tif_folder, ret_epoch, ret_path, drop_loc_only, ret_trial_drop):
    s=time.time()
    trial_traces_ret = get_trialwise_data_pydat(base_path, 'rolling_percentile', tif_folder, 5, ret_epoch,
                                   do_zscore=True, do_baseline=True, baseline_n = 4, stim_aligned=True)
    trial_traces_ret_sp = get_trialwise_data_pydat(base_path, '', tif_folder, 5, ret_epoch,
                                       do_zscore=True, get_spikes=True, stim_aligned=True)
    print('traces got in', time.time()-s)
    s=time.time()
    ret_df = traces_to_vis_df_py(trial_traces_ret, ret_path, 'retino', drop_locs = drop_loc_only)
    ret_df_sp = traces_to_vis_df_py(trial_traces_ret_sp, ret_path, 'retino', drop_locs = drop_loc_only)
    ret_df.loc[:,'sp'] = ret_df_sp['df'].values
    ret_df.loc[:,'abs_id'] = ret_df.x+ret_df.y*30
    ret_df = ret_df.drop(index=ret_df.index[ret_df.trial.isin(ret_trial_drop)])
    print('ret df got in', time.time()-s)
    return ret_df

def make_beh_df(base_path, tif_folder, beh_trials_df, beh_epoch, beh_trial_drop):
    traces_beh = get_trialwise_data_pydat(base_path, 'rolling_percentile', tif_folder, 5, beh_epoch, stim_aligned=True,
                                     do_zscore=True, do_baseline=True, baseline_n=4)
    traces_beh_nob = get_trialwise_data_pydat(base_path, 'rolling_percentile', tif_folder, 5, beh_epoch, stim_aligned=True,
                                     do_zscore=True, do_baseline=False, baseline_n=4)

    traces_beh_sp = get_trialwise_data_pydat(base_path, '', tif_folder, 5, beh_epoch, stim_aligned=True,
                                        get_spikes=True, do_zscore=True)
    
    print('cutting added trials off of beh_trials_df. ', len(beh_trials_df), ' trials exist, ', traces_beh.shape[0], ' tifs exist')
    beh_trials_df = beh_trials_df[beh_trials_df.index<traces_beh.shape[0]]

    df_sp = traces_to_beh_df(traces_beh_sp, beh_trials_df)
    df = traces_to_beh_df(traces_beh, beh_trials_df)
    df_nob = traces_to_beh_df(traces_beh_nob, beh_trials_df)
    df['sp'] = df_sp.df.values
    df['nob'] = df_nob.df.values

    print('droppings from both', beh_trial_drop)
    df = df[~df.trial.isin(beh_trial_drop)]
    beh_trials_df = beh_trials_df[~beh_trials_df.index.isin(beh_trial_drop)]

    return df, beh_trials_df


def display_beh_info(beh_trials_df, beh_daq_path, base_path):
    print_daq_notes(beh_daq_path)
    performance = beh_trials_df.groupby(['intensity']).Response.value_counts(normalize=True)
    performance = performance.xs(1,level='Response')
    performance = performance.reset_index(level=['intensity'])
    plt.plot(performance.intensity+.5,performance.Response, 'ro-', label='lick rate')
    plt.xscale('log')
    plt.title(str(len(beh_trials_df)) + ' Trials')
    plt.ylabel('% Lick')
    plt.xlabel('% Contrast')
    ax2= plt.gca().twinx()
    try:
        sns.lineplot(data=beh_trials_df, x='intensity', y='post_stim_lick', err_style='bars')
    except:
        print('NO DAQ DATA')
        sns.lineplot(data=beh_trials_df, x='intensity', y='RespTime', err_style='bars')
    ax2.set_xscale('log')
    plt.savefig(base_path + 'behavior_summary.png')