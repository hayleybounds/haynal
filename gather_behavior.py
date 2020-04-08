
import glob
import os
def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


mouse, day = 'HB67_2_209', '190910'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 5,2
ret_epoch_post = ret_path_post = None
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20190910-15_19_35_HB67_2_209_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/190910/HB67_2_209/HB67_2_209_182_001.mat'
oret_path = tif_folder+'190910_HB67_2_209_ret1_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/190910/HB67_2_209/190910_HB67_2_209_B.mat'



mouse, day = 'HB67_2_209', '191014'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 4,1
ret_epoch_post = ret_path_post = None
drop_loc_only = beh_trial_drop = drop_tri_only = []
ret_trial_drop = [0,1,2,3,4]
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191014-19_16_09_HB67_2_209_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/'+('/').join([day,mouse]) +'/HB67_2_209_225_000.mat'
oret_path = tif_folder+'191014_HB67_2_209_ret0_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191014/HB67_2_209/191014_HB67_2_209_C.mat'


mouse, day = 'HB67_2_209', '191016'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 2,1
ret_epoch_post = ret_path_post = None
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191016-18_50_48_HB67_2_209_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/191016/HB67_2_209/HB67_2_209_238_000.mat'
oret_path = tif_folder+'191016_HB67_2_209_ret0_d238_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191016/HB67_2_209/191016_HB67_2_209_A.mat'
beh_trial_drop = list(np.linspace(0,65,66))






"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""

mouse,day = 'HB67_3_217','190925'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 3,2
ret_epoch_post = ret_path_post = None
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20190925-16_04_48_HB67_3_217_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/190925/HB67_3_217/HB67_3_217_184_002.mat'
oret_path = tif_folder+'190925_HB67_3_217_r2_fov2_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/190925/HB67_3_217/190925_HB67_3_217_A.mat'


mouse, day = 'HB67_3_217', '191008'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 3,2
ret_epoch_post = ret_path_post = None
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191008-17_58_55_HB67_3_217_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/191008/HB67_3_217/HB67_3_217_205_001.mat'
oret_path = tif_folder+'HB67_3_217_ret1_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191008/HB67_3_217/191008_HB67_3_217_A.mat'



mouse, day = 'HB67_3_217', '191017'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 2, 1
ret_epoch_post = ret_path_post = None
ret_trial_drop = beh_trial_drop = drop_tri_only = []
drop_loc_only = [2,3]
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191017-17_56_39_HB67_3_217_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/191017/HB67_3/HB67_3_280_001.mat'
oret_path = tif_folder+'191017_HB67_3_217_ret1_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191017/HB67_3_217/191017_HB67_3_217_A.mat'

    
mouse,day = 'HB67_3_217', '191024'
total_base = '/datadrive/ContrastDetection/'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 2,1
ret_epoch_post = ret_path_post = None
ret_trial_drop = drop_loc_only = beh_trial_drop = []
drop_tri_only = [0,1]
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191024-21_06_14_HB67_3_217_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/191024/HB67_3_217/HB67_3_217_200ish_001.mat'
oret_path = tif_folder+'20191024_HB67_3_217_ret1_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191024/HB67_3_217/191024_HB67_3_217_A.mat'


mouse, day = 'HB67_3_217','191025'
total_base = '/datadrive/ContrastDetection/'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 2,1
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
ret_epoch_post = ret_path_post = None
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191025-15_22_16_HB67_3_217_trials.tsv'
ret_path ='/mnt/recurrence/frankenshare/ChroME Screens/StimData/191025/HB67_3_217/HB67_3_217_200_000.mat'
oret_path = tif_folder+'20191025_HB67_3_217_ret0_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191025/HB67_3_217/191025_HB67_3_217_A.mat'


mouse,day = 'HB67_3_217','191114'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 3,1
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
ret_epoch_post = ret_path_post = None
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191114-17_36_13_HB67_3_217_trials.tsv'
ret_path = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/191114/HB67_3_217/HB67_3_217_200_003.mat'
oret_path = tif_folder+'20191114_HB67_3_217_ret4_fov3_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191114/HB67_3_217/191114_HB67_3_217_A.mat'
   





"""%%%%%%%%%%%%%%%%%%%%% hb 67_4 %%%%%%%%%%%%%%%%%%%%%"""
mouse, day = 'HB67_4_210','190926'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 3,2
ret_epoch_post = 4
ret_trial_drop = drop_loc_only = beh_trial_drop = drop_tri_only = []
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20190926-16_27_11_HB67_4_210_trials.tsv'
ret_path ='/mnt/recurrence/frankenshare/ChroME Screens/StimData/190926/HB67_4_210/HB67_4_210_200_002.mat'
ret_path_post = '/mnt/recurrence/frankenshare/ChroME Screens/StimData/190926/HB67_4_210/HB67_4_210_200_003.mat'
oret_path = tif_folder+'20190926_HB67_4_210_ret2_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/190926/HB67_4_210/190926_HB67_4_210_D.mat'


mouse, day = 'HB67_4_210', '191018'
base_path = '/datadrive/ContrastDetection/'+mouse+'/'+day+'/'
base_path = base_path + get_immediate_subdirectories(base_path)[0]+'/suite2p/'
tif_folder = '/mnt/recurrence/frankenshare/hayley/'+mouse+'/'+day+'/'
beh_epoch, ret_epoch = 3,2
ret_trial_drop = []
drop_loc_only = [71]
beh_trial_drop = list(range(16))
drop_tri_only = [39]
ret_epoch_post = ret_path_post = None
tri_df_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetectionTask/20191018-16_38_47_HB67_4_210_trials.tsv'
ret_path ='/mnt/recurrence/frankenshare/ChroME Screens/StimData/191018/HB67_4/HB67_4_220_003.mat'
oret_path = tif_folder+'20191018_HB67_4_210_ret3_PSTHs.mat'
beh_daq_path = '/mnt/recurrence/frankenshare/hayley/ContrastDetection/191018/HB67_4_210/191018_HB67_4_210_A.mat'

