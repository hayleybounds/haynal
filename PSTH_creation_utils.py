
import pandas as pd
import scipy.io as sio
import numpy as np
from glob import glob
import pickle
from scipy.stats.mstats import zscore
from ScanImageTiffReader import ScanImageTiffReader

def chop_file_vars(starts, lengths, skip_ends, skip_starts):
    """for dropping elems of starts and lengths to ignore first or last acqs
    """
    if skip_ends > 0:
        starts = starts[0:len(starts)-skip_ends]
        lengths = lengths[0:len(starts)-skip_ends]
    if skip_starts > 0:
        starts = starts[skip_starts:]
        lengths = lengths[skip_starts:]
    return starts, lengths


def do_pre_dfof(traces, dfof_method, do_zscore):
    """do fluorescence calculations that should occur before chopping into PSTHS
    This occurs for percentile, rolling_percentile, and z scoring.
    """
    if dfof_method == 'percentile':
        f0=np.nanpercentile(traces,30, axis=1)
        f0 = np.reshape(f0,(f0.shape[0],1))
        traces = (traces-f0)/f0
    if dfof_method == 'rolling_percentile':
        f0s = pd.DataFrame(np.transpose(traces)).rolling(200, min_periods=1,center=True).quantile(.20)
        f0s = np.transpose(f0s.values)
        traces = (traces-f0s)/f0s
    if do_zscore:
        traces = zscore(traces, axis=1)
    return traces
        
def do_dfof_cuts(cut_traces, dfof_method, do_baseline, baseline_n):
    if dfof_method == 'at start':
        for i in range(0, np.shape(cut_traces)[1]):
            for j in range(0, np.shape(cut_traces)[0]):
                f0 = np.mean(traces[i,df_range[0]:df_range[1]])
                cut_traces[j,i,:] = (cut_traces[j,i,:]-f0)/f0
    elif dfof_method == 'by trial':
        for i in range(0, np.shape(cut_traces)[1]):
            for j in range(0, np.shape(cut_traces)[0]):
                f0 = np.mean(cut_traces[j,i,df_range[0]:df_range[1]])
                cut_traces[j,i,:] = (cut_traces[j,i,:]-f0)/f0
    elif 'percentile' not in dfof_method:
        print('no dfof method selected, returning raw cut traces')
    if do_baseline:
        print('subtracting baseline!')
        cut_traces  = cut_traces - np.nanmean(cut_traces[:,:,0:baseline_n],axis=2).reshape((cut_traces.shape[0],
                                                                                       cut_traces.shape[1],
                                                                                       1))
    return cut_traces

def get_traces_segment(real_start, real_end, traces, cell=None):
    """get a segment of fluorescent traces with nan padding
    Rteturns a segment of traces corresponding to a time period for a given cell (if cell is not None) or 
    all cells if cell is None. Pads with nans if a start is before the beginning of traces or an end is after
    the end of traces, but will throw an error if real_start exceeds 2nd dim of traces or real_end is less than 0.
    
    Args:
        real_start (int):  start index in frames of segment
        real_end (int): end index in frames of  segment
        traces (2d np array): array of F values, cells x frames 
        cell (int or array): cell or cells for which the F values should be returned. If none, returns for all cells.
    
    Returns:
    
    """
    if cell==None:
        cell = np.linspace(0, traces.shape[0]-1, traces.shape[0]).astype('int16')
    if real_start < 0 and real_end > 0:
        try:
            cut = np.concatenate((np.tile(np.nan, (len(cell),-1*real_start)), traces[cell,0:real_end]),1)
        except:
            print('shape of traces was: ', traces.shape, ' start was ', real_start, ' end was ', real_end,
                  'padded nans were', np.tile(np.nan, (len(cell),-1*real_start)).shape,
                 ' cut bit was ', traces[cell,0:real_end])
            raise Exception('unexpected error, see printed message for details')
    elif real_end < 0: #this should never happen
        raise Exception('requested a trace segment with end time exceeding length of traces')
    #if real_end longer than available data, but real start isn't
    elif real_end > traces.shape[1] and real_start < traces.shape[1]:
        try:
            cut = np.concatenate((traces[cell,real_start:],np.tile(np.nan, (len(cell),real_end-traces.shape[1]))),1)
        except:
            print('shape of traces was: ', traces.shape, ' start was ', real_start, ' end was ', real_end,
                  'padded nans were', np.tile(np.nan, (len(cell),real_end-traces.shape[1])).shape,
                 'cut bit was ', traces[cell,real_start:].shape)
            raise Exception('unexpected error, see printed message for details')
    elif real_end > traces.shape[1] and real_start >= traces.shape[1]: #this also should not happen
        raise Exception('requested a start time that exceeds the length of traces')
    else:
        cut = traces[cell, real_start:real_end]
    return cut