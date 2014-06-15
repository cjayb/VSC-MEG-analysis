# 
#
# License: BSD (3-clause)

import sys
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

import mne
from mne.fiff import Raw, pick_types

import numpy as np
import os, errno

do_postproc_univar = True

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def create_trigger_logic(cond):

    trig_mult_factor = {'CS+': 0, 'CS-': 0}
    if '1a' in cond:
        trig_mult_factor['CS+'] = 1
        trig_mult_factor['CS-'] = 2
    elif '1b' in cond:
        trig_mult_factor['CS+'] = 2
        trig_mult_factor['CS-'] = 1


    trig_logic = {}

    for cond in ['CS+', 'CS-']:
        trig_logic.update({cond: {}})
        for session in ['pre', 'post']:
            trig_logic[cond].update({session: {}})
            for imageType in ['VS','FB']:
                trig_logic[cond][session].update({imageType: {}})
                if imageType == 'VS':
                    trig_base = trig_mult_factor[cond]*100
                elif imageType == 'FB':
                    trig_base = trig_mult_factor[cond]*10
                for emotion in ['Neu','Ang']:
                    if emotion == 'Neu':
                        trig_code = trig_base
                    elif emotion == 'Ang':
                        # This is if we want to treat the position of the target separately!
                        # if imageType == 'VS':
                        #     trig_code = trig_base + np.arange(11,17)
                        # elif imageType == 'FB':
                        #     trig_code = trig_base + 1

                        # This way all Angry faces in the VS cond become 101/201...
                        trig_code = trig_base + 1
                    
                    trig_logic[cond][session][imageType].update({emotion: trig_code})

    return trig_logic

def contrast_logic():
    devsA, devsB = range(111,117), range(211,217)
    contrasts = ['VS','FB']
    clogic = dict(VS={}, FB={})
    clogic['VS'].update({'all': dict(all=[100,200] + devsA + devsB)})
    clogic['VS'].update({'face': dict(faceA=[100] + devsA, faceB=[200] + devsB)})
    clogic['VS'].update({'odd': dict(std=[100,200],dev = devsA + devsB)})
    clogic['FB'].update({'all': dict(all=[10, 20, 11, 21])})
    clogic['FB'].update({'face': dict(faceA=[10, 11], faceB=[20, 21])})
    clogic['FB'].update({'odd': dict(std=[10,20],dev = [11, 21])})

    return clogic

def events_logic(events, contrast):
    devsA, devsB = range(111,117), range(211,217)
    if 'face' == contrast:
        VS_eve = mne.pick_events(events, include=range(100,220))
        FB_eve = mne.pick_events(events, include=range(10,22))
        VS_eve = mne.merge_events(VS_eve, [100]+devsA, 1, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, [200]+devsB, 2, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [10,11], 1, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [20,21], 2, replace_events=True)

    return (VS_eve, FB_eve)

###############################################################################
# Set epoch parameters
tmin, tmax = -0.4, 0.6  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.2, 0.2  # reject trial only if blinks in the 400 ms middle portion!
baseline = (-0.2, 0.)
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13)
filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 35.0, 'highpass': 0.5}

filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])

if do_postproc_univar: # do a couple of "main effects"

    contrasts = contrast_logic() # read in some predefined contrasts

    for subj in ad.analysis_dict.keys():

        # Drop the FFA session for now, deal with it separately, also empty room
        session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys() 
                if ('FFA' not in x and 'empty' not in x)]

        raw_path = ad._scratch_folder + '/filtered/' + filter_params['input_files'] + '/' + filt_dir + '/' + subj
        eve_path = ad._scratch_folder + '/events.fif/' + subj + '/raw'

        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        mkdir_p(evo_path)

        for sesname in session_names:

            if '_1' == sesname[-2:]:
                session = 'pre'
            elif '_2' == sesname[-2:]:
                session = 'post'

            fname = raw_path + '/' + sesname + '_filt.fif'
            raw = Raw(fname, preload=False) 
            events = mne.read_events(eve_path + '/' + sesname + '-eve.fif') 

            picks = pick_types(raw.info, meg=True, eeg=False, stim=True, eog=True, misc=True)
            for trial_type in ['VS','FB']:

                for contrast in ['face','odd']:
                    
                    # here we need to make a method that takes in the events and 
                    # spits out mod'd events and id's  
                    event_id = contrasts[trial_type][contrast] # if dict, each value must be int
                    
                    # Don't do any rejection yet
                    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                                        baseline=baseline, reject=None, preload=True,
                                        reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings
                    # these are always contrasts, so do one minus the other 
                    con_keys = contrasts[trial_type][contrast].keys()
                    # average epochs and get an Evoked dataset.
                    evoked = epochs[con_keys[1]].average() -  \
                                epochs[con_keys[0]].average() 
                    cov = mne.compute_covariance(epochs, tmin=baseline[0], tmax=baseline[1])
                            
                    epo_out = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-epo.fif'
                    evo_out = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-ave.fif'
                    cov_out = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-cov.fif'
                    print evo_out
                    epochs.save(epo_out)  # save evoked data to disk
                    evoked.save(evo_out)  # save evoked data to disk
                    cov.save(cov_out)  # save evoked data to disk
                    
                    

###############################################################################
# View evoked response
# times = 1e3 * epochs.times  # time in miliseconds

# ch_max_name, latency = evoked.get_peak(mode='neg')

# import matplotlib.pyplot as plt
# evoked.plot()

# plt.xlim([times[0], times[-1]])
# plt.xlabel('time (ms)')
# plt.ylabel('Potential (uV)')
# plt.title('EEG evoked potential')

# plt.axvline(latency * 1e3, color='red', 
#             label=ch_max_name, linewidth=2,
#             linestyle='--')
# plt.legend(loc='best')

# plt.show()

