from ScanImageTiffReader import ScanImageTiffReader
import numpy as np
import os
from glob import glob
import pickle
from scipy import interpolate
"""
Utils for interacting with tiffs written by scan image. Getting metadata,
finding lengths, finding aux Triggers, ect.
"""

def get_prev_folder_length(base_path, epoch):
    ops = np.load(base_path + 'ops1.npy')
    lengths = ops[0]['frames_per_folder']
    return sum(lengths[0:epoch-1])


def get_nchannels(file):
    with ScanImageTiffReader(file) as reader:
        metadata = reader.metadata()
    channel_pass_1 = metadata.split('channelSave = [')

    if len(channel_pass_1)==1:
        nchannels = 1
    else:
        nchannels = len(metadata.split('channelSave = [')[1].split(']')[0].split(';'))
    return nchannels

def get_nvols(file):
    with ScanImageTiffReader(file) as reader:
        metadata = reader.metadata()
    #rint(metadata.split('hFastZ.userZs')[1])
    #rint(len(metadata.split('hFastZ.userZs')))
    #print(metadata.split('arbitraryZs')[1])
    
    try:
        if metadata.split('hFastZ.userZs = ')[1][0]=='0':
            return 1
        nvols = len(metadata.split('hFastZ.userZs = [')[1].split(']')[0].split(' '))
    except:
        nvols = len(metadata.split('hStackManager.arbitraryZs = [')[1].split(']')[0].split(' '))
    return nvols

def get_vol_locs(file, opto_cal = [0,3.5,6,7,9.5,11,13], um_cal = [0,15,26,30,45,54,66]):
    """ return the vol z locs, and an attempt at conversion to um """
    interper = interpolate.interp1d(opto_cal,um_cal, fill_value='extrapolate')
    with ScanImageTiffReader(file) as reader:
        metadata = reader.metadata()
    #rint(metadata.split('hFastZ.userZs')[1])
    #rint(len(metadata.split('hFastZ.userZs')))
    #print(metadata.split('arbitraryZs')[1])
    
    try:
        if metadata.split('hFastZ.userZs = ')[1][0]=='0':
            return 1
        vols = metadata.split('hFastZ.userZs = [')[1].split(']')[0].split(' ')
    except:
        vols = metadata.split('hStackManager.arbitraryZs = [')[1].split(']')[0].split(' ')
    vols = np.array(vols).astype('float16')
    return vols, interper(vols)

def get_tif_filename(tifFolder, base_path):
    """just get a single tif from the first folder to run other fxns on"""
    # bc the folders won't necessarily be in suite2p order, we need to get subfolders from base_path
    fold_names = base_path.split('/')[-3].split('_')
    print('these are my folders', fold_names)
    #if no underscore, only one folder
    folders = [tifFolder + fold_names[i] for i in range(0, len(fold_names))]
    tif = glob(folders[0]+'/*.tif')[0]
    return tif

def get_starts_multifolder(tifFolder, base_path, get_stim_starts=True):
    """ get tif lengths and starts for all epochs in base_path, by looking at tiffs in tifFolder
    """
    # bc the folders won't necessarily be in suite2p order, we need to get subfolders from base_path
    fold_names = base_path.split('/')[-3].split('_')
    print('these are my folders', fold_names)
    #if no underscore, only one folder
    folders = [tifFolder + fold_names[i] for i in range(0, len(fold_names))]
    full_file_starts = []
    full_lengths = []
    full_stim_starts = []
    for folder in folders:
        print('looking in folder ', folder)
        print('and going to save to ', base_path + os.path.split(folder)[1] + 'PyFileData.p')
        stim_starts, lengths, file_starts = get_tif_data(folder + '/', base_path + os.path.split(folder)[1] + 'PyFileData.p', get_stim_starts=get_stim_starts)
        full_file_starts.append(file_starts)
        full_lengths.append(lengths)
        full_stim_starts.append(stim_starts)

    return full_file_starts, full_lengths,full_stim_starts

def get_starts_suite2p(tifFolder, base_path, epoch, get_stim_starts=True):
    fold_names = base_path.split('/')[-3].split('_')
    #if no underscore, only one folder
    folders = [tifFolder + fold_names[i] for i in range(0, len(fold_names))]
    folder = folders[epoch-1]
    stim_starts, lengths, file_starts = get_tif_data(folder + '/', base_path + os.path.split(folder)[1] + 'PyFileData.p', get_stim_starts=get_stim_starts)
    return file_starts, lengths, stim_starts

def get_fr(tifFolder):
    """checks fr for the first tif in a folder
    """
    files = glob(tifFolder+'*.tif')
    if len(files)==0:
        print('no tifs found in ', tifFolder)
        return
    
    with ScanImageTiffReader(files[0]) as reader:
        FR = float(reader.metadata().split('scanVolumeRate = ')[1].split('\n')[0])
    return FR

def get_tif_data(tifFolder, picklefile, matfile = None, get_stim_starts = True):
    """if there is not already a picklefile with tif data, combs through tifs to get lengths and stim starts.

    Only works on one folder of tifs at a time, use multifolder to coordinate for python.
    Currently only looks at auxTrigger 1.

    Inputs:
        tifFolder (str): path to folder containing tifs
        picklefile (str): where to look for picklefile for already generated data. If not found, will save new file here.
        matfile (optional, str): not used anymore mostly, but will save data to mat file as well
        get_stim_starts (bool, default True): whether to comb frame by frame for aux Trigger data. Goes faster if false.

    Outputs:
        stim_starts: if get_stim_starts is True, the aux triggers. Otherwise, empty. corrected to be in volume units
        file_lengths: length of each tif file, in volumes.
    """
    files=glob(picklefile)

    if len(files) > 0:
        print('found pickle file')
        with open(files[0],'rb') as f:
            [file_lengths, stim_starts, file_starts] = pickle.load(f)
            return stim_starts, file_lengths, file_starts
    else:
        files = glob(tifFolder+'*.tif')
        if len(files)==0:
            print('no tifs found in ', tifFolder)
            return
        else:
            print(len(files), ' tifs found in ', tifFolder)
        files.sort()
        file_lengths = [];
        stim_on_frame = [];
        file_starts = [0];
        FR=None

        for file in files:
            print(file)
            with ScanImageTiffReader(file) as reader:
                file_lengths.append(reader.shape()[0])
                file_starts.append(file_starts[len(file_starts)-1]+reader.shape()[0])
                if get_stim_starts:
                    stim_on = []

                    for framei in range(0,reader.shape()[0]):
                        string = reader.description(framei).split('auxTrigger')[1].split('[')[1]
                        if string[0]!=']' and string[0]!='-':
                            stim_on.append(framei)
                    if len(stim_on)>2:
                        print('oh no toooo many stims!!')
                    if len(stim_on) > 0:
                        stim_on_frame.append(min(stim_on))
                    else:
                        print('oh no no stim')
                        stim_on_frame.append(np.nan)

                FR = float(reader.metadata().split('scanVolumeRate = ')[1].split('\n')[0])
                
                #int(reader.metadata().split('hBeams.powers = ')[1].split('\n')[0])
                
        nchan = get_nchannels(files[0])
        nvol = get_nvols(files[0])
        file_starts = file_starts[0:len(file_starts)-1]
        if get_stim_starts:

            stim_starts = np.asarray(file_starts) + np.asarray(stim_on_frame)
            stim_starts = (stim_starts/(nchan*nvol)).astype('uint16')
        else:
            stim_starts = []

        file_starts = (np.asarray(file_starts)/(nchan*nvol)).astype('uint16')
        file_lengths = (np.asarray(file_lengths)/(nchan*nvol)).astype('uint16')

        dat = [file_lengths, stim_starts, file_starts]
        with open(picklefile, 'wb') as f:
            pickle.dump(dat, f)
        print('saved data!')
        if matfile is not None:
            saveme = {'fileStart':file_starts+1, 'fileLength':file_lengths,'stimStart':stim_starts, 'FR':FR}
            sio.savemat(matfile,saveme)
        return stim_starts, file_lengths, file_starts
    
def get_mimg_data(f):
    """given a scanimage data file, get the mean image per plane and channels"""
    with ScanImageTiffReader(f) as reader:
        data = reader.data();
    data = data.astype('float32')
    nvols = get_nvols(f)
    nchan = get_nchannels(f)
    mimgs = np.zeros((512,512,3,nvols),'float32')
    count=0
    chan_conv = [1,0] #convert gr to rg
    for iplane in range(nvols):
        for ichan in range(nchan):
            mimg = np.nanmean(data[count:data.shape[0]:nvols*nchan,:,:],axis=0)
            
            mimgs[:,:,chan_conv[ichan],iplane] = mimg
            count += 1
    return mimgs