from suite2p_pipeline_git import *

seagate = '/media/hbounds/Seagate Backup Plus Drive/'
frankenlocal = '/mnt/frankenlocal/Hayley/'

#process_data('HB41', '190806',['ret0','ret1','beh0','beh1','beh2','beh3','beh4','beh5'])
#process_data('HB41','190814',['ret0','beh0'])
#process_data('HB50.1','190703',['ret2','ori004'])

#process_data('HB53.3 1-100 syn 0136','190903',['noise','ret','stim'])
#process_data('HB53.3 1-100 syn 0136','190829',['stim']) #not a real screen day
#process_data('HB53.3 1-100 syn 0136','190805',['noise','stim2']) #not sure which are which fov tbh
#process_data('HB53.3 1-100 syn 0136','190729',['noise','stim'])
#process_data('HB53.3 1-100 syn 0136','190726',['noise'])
#process_data('HB53.3 1-100 syn 0136','190719',['noise','stim'])
#process_data('HB53_3','191001',['noise','stim'])


#process_data('HB53.3 1-100 syn 0136','190715',['noise','stim'])
#process_data('HB53.3 1-100 syn 0136','190910',['stim','stim_lwr','noise'])

#process_data('HB53.2 1-1000 syn 132','190908',['noise'])

#process_data('HB53.2 1-1000 syn 132','190805',['noise'])
#process_data('HB53.2 1-1000 syn 132','190729',['noise','stim'])
#process_data('HB53.2 1-1000 syn 132','190726',['noise'])
#process_data('HB53.2 1-1000 syn 132','190719',['noise','stim'])
#process_data('HB53.2 1-1000 syn 132','190715',['noise'])



#190908 has screen ims but nothing more, retake
#process_data('HB33_3 09705','190729',['noise'])
#process_data('HB33_3 09705','190725', ['noise','stim'])
#process_data('HB33_3 09705','190723',['noise'])
#process_data('HB33_3 09705','190719',['noise']) #noise is weirdly big>


#process_data('HB75_1','191004',['1','2'], raw_base = '/datadrive/tiffData/')
#process_data('HB75_1','191008',['masks','noise','ret'], raw_base = '/datadrive/tiffData/')
#process_data('HB75_1','191008',['maskscptemp','noise','ret'], raw_base = '/datadrive/tiffData/')
#process_data('HB75_1','191017',['ret0','ret1','stim1','stim3'])
#process_data('W20_1_alt','200106',['stim','5n10','1n2'])

#process_data('W20_2','200112',['pcurve','stim3'], raw_base = seagate)

#process_data('I138_3','200110',['stim0','singmult','1spike','5spike'], raw_base = seagate)


#process_data('W20_2','191219',['stim','pcurve','1n2'], raw_base = seagate)


#process_data('HB33_2 09710','190719',['noise'])
#process_data('HB33_2 09710','190725',['noise'])


#process_data('HB81','191109',['ret5','ret6','noise','stim1','stim2'])


#process_data('HB66_2', '191104',['vis','vis2','vis3','stim'])
#process_data('HB66_2', '191109',['noise1','ret0','stim1'])
#process_data('HB66_2','191125',['ret0','ret1','noise2','stim1'])
#process_data('HB66_2', '191203',['ret0','ret1','noise','stim','5n102'])

#process_data('HB66_1','191014',['ret0','stim'])
#process_data('HB66_1','191021',['stim'])
#process_data('HB66_1','191112',['stim','stim5n10','stim5n10p2'],diameter=8)
#process_data('HB66_1','191204',['stim','5n10'], raw_base='/media/hbounds/Seagate Backup Plus Drive/')
#process_data('HB66_1','191205',['stim','5n10'], raw_base = '/datadrive/tiffData/')
#process_data('HB66_1','191213',['stim2','5n102'], raw_base = '/datadrive/tiffData/')
#process_data('HB66_1','200130',['stim'],raw_base='/datadrive/tiffData/')

#process_data('HB66_1','200131',['stim0','stim1'],raw_base='/datadrive/tiffData/')
#DONT HAVE BINARY DAMMIT
#process_data('HB66_1','200204',['stim1','tentarg2'], raw_base ='/mnt/recurrence/frankenshare/hayley/',result_base='/datadrive/newOT/')
#process_data('HB66_1','200207',['stim0','tentargs','tentargs3'], raw_base = '/mnt/recurrence/frankenshare/hayley/', result_base='/datadrive/newOT/')
#process_data('HB66_1','200221',['noise2','stim0','pcurve1','tentargs2'],raw_base = '/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')

#process_data('HB50_1', '191112', ['ret1','ori3'])
#process_data('HB50_1', '191112', ['stim'])
#process_data('HB50_1','191212', ['ret0', 'noise1', 'stim', '5n10'], raw_base = seagate)
#process_data('HB50_1','191219',['ret1','noise'], raw_base = seagate)
#process_data('HB50_1','191219',['ori4','ret5','noise6'],raw_base = seagate)
#process_data('HB50_1','191219',['stim2'], raw_base = seagate)

#process_data('HB50_1','191219_pr',['pic_reg'], raw_base = seagate)
#process_data('HB50_1','191219_pr2',['pic_reg'], raw_base = seagate)

#process_data('HB66_1','191119',['stim'])


#process_data('HB82_1','191203_and_05',['03stim','05noise','05stim'])

#don't have binary fuck
#process_data('HB95_1','200214',['stim','tentarg'],raw_base = '/mnt/recurrence/frankenshare/hayley/',result_base='/datadrive/newOT/')

#do have bin
process_data('HB95_1','200310',['stim','pcurve','lowhighsat','noise','tentarg2'],raw_base = '/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')
#process_data('HB95_1','200316',['stim0','noise2','pcurve','pcurveC','pbaltentarg','pbaltentarghigh'], raw_base='/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')

#process_data('HB80_3','200313e',['noise','stim0','pcurve','tentarg'], raw_base = '/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')

#process_data('HB80_3','200314b',['noise','stim0','stim2','pcurve','tentarg'], raw_base='/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')
#process_data('HB80_3','200315',['noise2','stim3','stim4','pcurve','pcurve2','lowhighsat','tentargs'],raw_base='/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')
#process_data('HB80_3','200316d',['noise2','stim0','pcurve','pbaltentarg'],raw_base='/mnt/frankenlocal/Hayley/',result_base='/datadrive/newOT/')



"""TO RUN CHROME TEST 2
"""
# process_data('HB33_2 09710','190910',['noise','stim'])
# process_data('HB33_2 09710','190729',['noise'])
# process_data('HB33_2 09710','190719',['noise'])
# process_data('HB33_2 09710','190910',['noise','stim'])


"""Contrast Task Data"""


# process_data('HB67_4_210','190909',['ret0','ret1','beh0'],
# 	raw_base = '/datadrive/tiffData/',
# 	result_base = '/datadrive/ContrastDetection/')

# process_data('HB67_4_210','190925',['ret0','ret1','beh0'],
# 	raw_base = '/datadrive/tiffData/',
# 	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_4_210','190926',['ret1','ret2','beh0','ret3'],
#	raw_base = '/mnt/recurrence/frankenshare/hayley/',
#	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_4_210','191018',['ret1',ret3','beh0'],
# 	raw_base = '/mnt/recurrence/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')

# #GARBAGE DAYS, testing new diam?
# process_data('HB67_4_210','191021',['ret1','ret2','ret3','beh2','beh4'],
# 	raw_base = '/mnt/recurrence/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/',
# 	diameter=8)


#unclear if this worked?
# process_data('HB67_4_210','191025',['ret0','ret1', 'ret2', 'beh0', 'beh1'],
# 	raw_base = seagate,
# 	result_base = '/datadrive/ContrastDetection/', diameter=8)




#############HB67.2################

#process_data('HB67_2_209','190909',['ret0','beh0','beh1'],
#	raw_base = '/mnt/recurrence/frankenshare/hayley/',
#	result_base = '/datadrive/ContrastDetection/')


# process_data('HB67_2_209','190910',['ret0','ret1','ret2','beh0','beh1'],
# 	raw_base = '/mnt/modulation/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')

# process_data('HB67_2_209','191014',['ret0','beh0','beh1','beh2'],
# 	raw_base = '/datadrive/tiffData/',
# 	result_base = '/datadrive/ContrastDetection/')

# process_data('HB67_2_209','191016',['ret0','beh0'],
# 	raw_base = '/mnt/recurrence/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_2_209','191114',['ret1','ret2','ret3','beh0','beh1','ori4'],
#	raw_base = '/media/hbounds/Seagate Backup Plus Drive/',
#	result_base = '/datadrive/ContrastDetection/')



#########################HB67.3##################33333

# process_data('HB67_3_217','190905',['ret2','ret3','ret4','ret5','beh0','beh1'],
# 	raw_base = '/mnt/modulation/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_3_217','191114',['ret3','ret4','beh0'],
# 	raw_base = '/datadrive/tiffData/',
# 	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_3_217','190925',['ret1','ret2','beh0'],
#	raw_base = '/mnt/modulation/frankenshare/hayley/',
#	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_3_217','191008',['ret0','ret1','beh1'],
#	raw_base = '/datadrive/tiffData/',
#	result_base = '/datadrive/ContrastDetection/')

#process_data('HB67_3_217','191017',['ret1','beh0'],
# 	raw_base = '/mnt/recurrence/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')


# process_data('HB67_3_217','191024',['ret1','beh0'],
# 	raw_base = '/mnt/recurrence/frankenshare/hayley/',
# 	result_base = '/datadrive/ContrastDetection/')

# process_data('HB67_3_217','191025',['ret0','beh0'],
#  	raw_base = seagate,
#  	result_base = '/datadrive/ContrastDetection/')
