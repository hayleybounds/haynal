

import numpy as np
import sys
import suite2p
from suite2p.run_s2p import run_s2p
from suite2p.io import save 
import timeit
import os
import shutil
import pdb
import numpy
from importlib import reload
import glob
import json
from ScanImageTiffReader import ScanImageTiffReader


data = '/mnt/modulation/frankenshare/hayley/HB41/190731/'
fast_disk = '/home/hbounds/Suite2pTemp/'

# set your options for running
# overwrites the run_s2p.default_ops
ops = {
        'fast_disk': '',#new changed from fast_disk to this # used to store temporary binary file, defaults to save_path0 (set as a string NOT a list)
        'save_path0': [], # stores results, defaults to first item in data_path
        'delete_bin': False, # whether to delete binary file after processing
        # main settings
        'nplanes' : 3, # each tiff has these many planes in sequence
        'nchannels' : 1, # each tiff has these many channels per plane
        'functional_chan' : 1, # this channel is used to extract functional ROIs (1-based)
        'sparse_mode': False,
        'diameter': 6, # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
        'tau':  1.5, # this is the main parameter for deconvolution
        'fs': 6.25,  # sampling rate (total across planes)
        # output settings
        'save_mat': True, # whether to save output as matlab files
        'combined': True, # combine multiple planes into a single result /single canvas for GUI
        # parallel settings
        'num_workers': 11, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
        'num_workers_roi': 11, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
        # registration settings
        'do_registration': False, # whether to register data
        'nonrigid': False,
        'keep_movie_raw': False,

        'nimg_init': 1000, # subsampled frames for finding reference image
        'batch_size': 500,#3000, # number of frames per batch
        'maxregshift': .15, #0.05, # max allowed registration shift, as a fraction of frame max(width and height)
        'align_by_chan' :  1, # when multi-channel, you can align by non-functional channel (1-based)
        'reg_tif': True, # whether to save registered tiffs
        'reg_tif_chan2':False,
        'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
        'do_bidiphase': True, #whether to compute bidirectional phase offset
        # cell detection settings
        'connected': True, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
        'navg_frames_svd': 5000, # max number of binned frames for the SVD
        'nsvd_for_roi': 2000, # changed 12/9/19 from 1k to 3k max number of SVD components to keep for ROI detection
        'max_iterations': 20, # maximum number of iterations to do cell detection
        'ratio_neuropil': 6., # ratio between neuropil basis size and cell radius
        'ratio_neuropil_to_cell': 3, # minimum ratio between neuropil radius and cell radius
        'tile_factor': 1., # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
        'threshold_scaling': 2., # new with git, the default # adjust the automatically determined threshold by this scalar multiplier
        'max_overlap': 0.75, # cells with more overlap than this get removed during triage, before refinement
        'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
        'outer_neuropil_radius': np.inf, # maximum neuropil radius
        'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
        # deconvolution settings
        'baseline': 'maximin', # baselining mode
        'win_baseline': 60., # window for maximin
        'sig_baseline': 10., # smoothing constant for gaussian filter
        'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
        'neucoeff': .7 # neuropil coefficient

        #new!!
        #'preclassify':0.5,
      }



def process_data(animalid,date,expt_ids,
    raw_base='/mnt/recurrence/frankenshare/ChroME Screens/',
    result_base='/datadrive/ChroME Tests2/',
    fast_disk= '/datadrive/Suite2pTemp2/',
    delete_raw=False,diameter=6,
    db = {}):
#as of 1/3/19 diameter is 8
#    save_path0 = result_base+animalid+'/'+date+'/'+'_'.join(expt_ids)
#    data_path = [raw_base+animalid+'/'+date+'/'+lbl for lbl in expt_ids]
    delete_disk = False
    db = prepare_db(animalid,date,expt_ids,db,raw_base=raw_base,
        result_base=result_base,fast_disk=fast_disk,diameter=diameter)
    print('got that db')

    if delete_disk:
        try:
            shutil.rmtree(fast_disk+'/suite2p')
            print('fast disk contents deleted')
        except:
            print('fast disk location empty')


    print('now to suite2p')
    ops['cut_edges'] = False
    ops['cut_range'] = np.linspace(0,511,512).astype('uint16')

    opsEnd=run_s2p(ops=ops,db=db)
    
    
def prepare_db(animalid,date,expt_ids,db,
    raw_base, result_base, fast_disk, diameter):
    save_path0 = result_base+animalid+'/'+date+'/'+'_'.join(expt_ids)
    fast_disk = fast_disk+animalid+'/'+date+'/'+'_'.join(expt_ids)
    data_path = [raw_base+animalid+'/'+date+'/'+lbl for lbl in expt_ids]
    print(data_path)

    with ScanImageTiffReader(glob.glob(data_path[0]+'/*.tif')[0]) as reader:
        metadata = reader.metadata()

    print('got that metadata')
    # provide an h5 path in 'h5py' or a tiff path in 'data_path'
    # db overwrites any ops (allows for experiment specific settings)
    db['h5py'] = [] # a single h5 file path
    db['h5py_key']= 'data'
    db['look_one_level_down']= False # whether to look in ALL subfolders when searching for tiffs,
    db['save_path0']= save_path0
    db['data_path'] = data_path # a list of folders with tiffs 
                                 # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)                                     
    db['subfolders']= [] # choose subfolders of 'data_path' to look in (optional)
    db['fast_disk']= ''#NEW WITH GET REMOVED fast_disk # string which specifies where the binary file will be stored (should be an SSD)
    channel_pass_1 = metadata.split('channelSave = [')
    if len(channel_pass_1)==1:
        db['nchannels'] = 1
    else:
        db['nchannels']= len(metadata.split('channelSave = [')[1].split(']')[0].split(';'))
    try:
        db['nplanes']= len(metadata.split('hFastZ.userZs = [')[1].split(']')[0].split(' '))
    except:
        db['nplanes']=len(metadata.split('hStackManager.arbitraryZs = [')[1].split(']')[0].split(' '))
    db['diameter']= diameter
    db['fs'] =float(metadata.split('scanVolumeRate = ')[1].split('\n')[0])


    print('returning that db')
    return db

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    process_data(*sys.argv[1:])
