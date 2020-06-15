import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tiff_utils import *
from PSTH_creation_utils import *
import warnings
import time

def get_trialwise_data_pydat(base_path, dfof_method, tifFolder, pre_time, epoch, is_online_data=False, get_spikes = False,
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

    fstarts, full_lengths, stim_starts = get_starts_suite2p(tifFolder, base_path, epoch, get_stim_starts=stim_aligned)
    
    if stim_aligned:
        fstarts = stim_starts
    if full_lengths[0]==1:
        print('dropping first file with one frame')
        fstarts = fstarts[1:]
        full_lengths = full_lengths[1:]
    
    if not is_online_data and epoch>1:
        fstarts = fstarts + get_prev_folder_length(base_path, epoch)

    #percentile dff method is based on whole trace, so do now
    traces = do_pre_dfof(traces, dfof_method, do_zscore)
    length = min(full_lengths)

    #initialize psf array - it has to only be for one epoch bc number of trials varies
    cut_traces = np.zeros((len(fstarts), np.shape(traces)[0], length))

    
    #for each trial get the file start, and the relative stim start, and add it to traces
    for istim, stim_start in enumerate(fstarts):
        real_start = stim_start-pre_time
        real_end = real_start+length
        cut = get_traces_segment(real_start, real_end, traces)
        cut_traces[istim,:,:] = cut

    #then dfof and or baseline
    cut_traces = do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n)
    return cut_traces

"""%%%%%%% common pandas funcs%%%%%%%%%"""

def add_cell_info(df, base_path,fix_offset=False):
    """adds median and plane info from suite2p to a df"""
    planes = get_cell_field(base_path, 'iplane')
    medians = get_cell_field(base_path, 'med')
    if fix_offset:
        for dpth in range(0,4):
            if dpth==1:
                xoff,yoff = 512,0
            elif dpth==2:
                xoff,yoff=0,512
            elif dpth==3:
                xoff,yoff=512,512
            else:
                xoff,yoff=0,0
            medians[planes==dpth,0]-=yoff
            medians[planes==dpth,1]-=xoff
    s=time.time()
    cell_df = pd.DataFrame({'cell':list(range(len(planes))),'plane':planes,'mediany':medians[:,0],'medianx':medians[:,1]})
    df = df.merge(cell_df, on='cell',how='left')
    print('celldat merge time', time.time()-s)
    #cell_fields = np.hstack([planes.reshape(len(planes),1), medians])
    #vectors = np.zeros((len(df),cell_fields.shape[1]))
    #for cell in df.cell.unique():
    #    vectors[df.cell==cell,:] = cell_fields[cell,:]

    #df.loc[:,'plane'] = vectors[:,0]
    #df.loc[:,'mediany'] = vectors[:,1]
    #df.loc[:,'medianx'] = vectors[:,2]
    return df

def cut_to_melted(traces):    
    #change traces into a pandas array with cell and time as cols
    warnings.simplefilter(action='ignore', category=FutureWarning)
    frame = pd.Panel(np.transpose(traces)).to_frame(filter_observations = False)
    frame = frame.reset_index(level=['major','minor'])
    all_melted = pd.melt(frame, ('major','minor'))
    all_melted = all_melted.rename(columns = {'major':'cell', 'minor':'trial', 'variable':'time',
                       'value': 'df'})
    return all_melted

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555suite2p funcs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55"""


def get_iscell_F_Neu(base_path, nplanes=3, combined=True):
    """extract useful suite2p info for all the planes. Create F_subt, neuropil subtracted traces for id'ed cells
    returns iscell, F, Neu, Sp, F_subt
    """
    if not combined:
        print("HEY SOMETHING IS SUPER WRONG WITH TIMING HERE")
        iscell = [np.load(base_path + 'plane' + str(iplane)+'/iscell.npy') for iplane in range(nplanes)]
        iscell = np.concatenate(iscell)


        F = [np.load(base_path + 'plane' + str(iplane)+'/F.npy') for iplane in range(0,nplanes)]
        if min([f.shape[1] for f in F]) != max([f.shape[1] for f in F]):
            min_len = min([f.shape[1] for f in F])
            F = [f[:,:min_len] for f in F]


        Sp = [np.load(base_path + 'plane' + str(iplane)+'/spks.npy') for iplane in range(0,nplanes)]
        if min([f.shape[1] for f in Sp]) != max([f.shape[1] for f in Sp]):
            min_len = min([f.shape[1] for f in Sp])
            Sp = [f[:,:min_len] for f in Sp]



        Neu = [np.load(base_path + 'plane' + str(iplane)+'/Fneu.npy') for iplane in range(0,nplanes)]
        if min([n.shape[1] for n in Neu]) != max([n.shape[1] for n in Neu]):
            min_len = min([n.shape[1] for n in Neu])
            Neu = [n[:,:min_len] for n in Neu]

        #print("DOING  HACKY TIMET HING")
        #for f_vals in [Neu,F,Sp]:
        #    for i, plane_vals in enumerate(f_vals):
        #        if i>0:
        #            print(plane_vals.shape)
        #            plane_vals = np.roll(plane_vals,i, axis=1)

                    #plane_vals[:,0:i] = np.nan
        #        f_vals[i] = plane_vals


        Neu = np.vstack(Neu)
        F = np.vstack(F)
        Sp = np.vstack(Sp)

    else:
        iscell = np.load(base_path+'combined/iscell.npy')

        F = np.load(base_path+'combined/F.npy')
        Neu = np.load(base_path+'combined/Fneu.npy')
        Sp = np.load(base_path+'combined/spks.npy')

    Sp = Sp[iscell[:,0].astype('bool'),:]

    F_subt = F[iscell[:,0].astype('bool'),:] - .7*Neu[iscell[:,0].astype('bool'),:]

    return iscell, F, Neu, Sp, F_subt


def get_motion_data(base_path):
    """returns corrXY, xoff, and yoff for all planes

    Input: 
        base_path: path to where ops1.npy can be found
    Output: 
        corrXY, xoff and yoff: nplanes x nframes arrays of the suite2p values
    """

    ops = np.load(base_path + 'ops1.npy')
    motion_props = [[ops[i][prop] for i in range(len(ops))] for prop in ('corrXY','xoff','yoff')]
    if min([f.shape[0] for f in motion_props[0]]) != max([f.shape[0] for f in motion_props[0]]):
            min_len = min([f.shape[0] for f in motion_props[0]])
            motion_props = [[f[:min_len] for f in m] for m in motion_props]
    motion_props = [np.vstack(m) for m in motion_props]

    return [mp for mp in motion_props]


def get_cell_field(base_path, field, only_cells=True):
    """Retrieves a data field for identified cells. Returns a numpy array of the data.

    Inputs:
        base_path (str): path to suite2p folder where 'combined' folder exists
        field (str): data to extract (eg. med, iplane)
        
    KwArgs:
        only_cells (bool): default True, whether to only take info for is_cell values
    """
    cell_stats = np.load(base_path+'combined/stat.npy', allow_pickle=True)
    info = [cell_info[field] for cell_info in cell_stats]
    info = np.array(info)
    iscell = np.load(base_path+'combined/iscell.npy', allow_pickle=True)
    
    if only_cells:
        if len(info.shape)>1:
            info = info[iscell[:,0].astype(bool),:]
        else:
            info = info[iscell[:,0].astype(bool)]

    return info

def get_suite2p_mimg(base_path):
    all_ops = np.load(base_path+'ops1.npy')
    mimg = np.zeros((512,512,3,len(all_ops)))
    
    for iplane in range(len(all_ops)):
        ops = all_ops[iplane]
        mimg_green = ops['meanImg']
        mimg_red = ops['meanImg_chan2']
        mimg[:,:,0,iplane] = mimg_red
        mimg[:,:,1,iplane] = mimg_green 
    return mimg.astype('float32')

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% image manipulations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 """

def scale_array(mimg, ranges):
    mimg=mimg.copy().astype('float32')
    mimg[mimg<ranges[0]] = ranges[0]
    
    mimg -= ranges[0]
    mimg[mimg>ranges[1]-ranges[0]] = ranges[1]-ranges[0]
    mimg /= (ranges[1]-ranges[0])
    return mimg

def scale_mimg(mimgs, ranges):
    mimgs = mimgs.copy()
    chan_conv = [1,0]
    for iplane in range(mimgs.shape[3]):
        for ichan in range(min(mimgs.shape[2],2)):
            if mimgs.shape[2]>1:
                mimg = mimgs[:,:,chan_conv[ichan],iplane]
            else:
                mimg = mimgs[:,:,ichan,iplane]
            mimg[mimg<ranges[ichan][0]] = ranges[ichan][0]
            mimg[mimg>ranges[ichan][1]] = ranges[ichan][1]
            mimg -= ranges[ichan][0]
            #print(mimg.max())
            mimg /= (ranges[ichan][1]-ranges[ichan][0])
            #print(mimg.max())
            if mimgs.shape[2]>1:
                mimgs[:,:, chan_conv[ichan],iplane]=mimg
            else:
                mimgs[:,:,ichan,iplane]=mimg
    return mimgs

def scale_mimg_auto(mimgs,prctiles = [5,95]):
    mimgs=mimgs.copy()
    chan_conv = [1,0]
    for iplane in range(mimgs.shape[3]):
        for ichan in range(2):
            mimg = mimgs[:,:,chan_conv[ichan],iplane]
            ranges = [np.percentile(mimg,prctiles[0]),np.percentile(mimg,prctiles[1])]
            mimg = scale_array(mimg, ranges)
            mimgs[:,:, chan_conv[ichan],iplane]=mimg
    return mimgs


"""pandas manipulations"""

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

def mean_by_trial(df, start, stop, add_columns=[]):
    """given df with cell + trial columns as well as any cols in add_columns,
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

def mean_of_triwise_mean(mean_post,col='power'):
    #gives a single value per stim value, rather than per trial
    mean_diffs = mean_post.groupby(['cell',col]).mean()
    sem_diffs = mean_post.groupby(['cell',col]).df.sem()

    mean_diffs = mean_diffs.reset_index(level=['cell',col])
    mean_diffs['sem'] = sem_diffs.values
    return mean_diffs

"""%%%%%%%%%%%%%%%%%%%%%555 plotting things %%%%%%%%%%%%%%%%%%%%%%%%%%"""
def figsize(w,h):
    plt.gcf().set_size_inches(w,h)
    
def get_maxproj_fullsize(ops):
    mimg = np.zeros((512,512))
    mimg[ops['yrange'][0]:ops['yrange'][1],ops['xrange'][0]:ops['xrange'][1]] = ops['max_proj']
    return mimg


def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    from: https://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
    example: light_seismic = cmap_map(lambda x: x/1.2 + 0.3, matplotlib.cm.seismic)
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)