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

###############################################################################
# Set epoch parameters
tmin, tmax = -0.5, 1.5          # This wide a window might also be good for detecting responses?
rej_tmin, rej_tmax = -0.2, 0.2  # reject trial only if blinks in the 400 ms middle portion!
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13)
filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 35.0, 'highpass': 0.5}

filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])


for subj in ad.analysis_dict.keys():

    cond_names = ad.analysis_dict[subj][filter_params['input_files']].keys()

    raw_path = ad._scratch_folder + '/filtered/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
    eve_path = ad._scratch_folder + '/events.fif/' + subj + '/raw'

    evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
    mkdir_p(evo_path)

    trig_logic = {}

    for cond in cond_names:

        if 'VS' in cond: # do the FFA's separately
            if '_1' == cond[-2:]:
                session = 'pre'
            elif '_2' == cond[-2:]:
                session = 'post'

            if not any(trig_logic):
                trig_logic = create_trigger_logic(cond)

            fname = raw_path + '/' + cond + '_filt.fif'
            raw = Raw(fname, preload=False) 
            orig_events = mne.read_events(eve_path + '/' + cond + '-eve.fif') 

            # replace 111-116 with 101 and 211-216 with 201
            events = mne.merge_events(orig_events, np.arange(111,117), 101, replace_events=True)                    
            events = mne.merge_events(events, np.arange(211,217), 201, replace_events=True)                    
            
            picks = pick_types(raw.info, meg=True, eeg=False, stim=True, eog=True, misc=True)
            for stim_type in ['CS+','CS-']:
                for trial_type in ['VS','FB']:

                    #add all events!
                    # first grand-average and covariance!
                    event_id = {'Neu': trig_logic[stim_type][session][trial_type]['Neu'],
                                'Ang': trig_logic[stim_type][session][trial_type]['Ang']}
                    
                    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                                        baseline=(None, 0), reject=reject, preload=True,
                                        reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings
                    evoked = epochs.average()  # average epochs and get an Evoked dataset.
                    cov = mne.compute_covariance(epochs, tmin=None, tmax=0.)
                            
                    evo_out = evo_path + '/' + trial_type + '_' + stim_type + '_' + session + '-ave.fif'
                    cov_out = evo_path + '/' + trial_type + '_' + stim_type + '_' + session + '-cov.fif'
                    print evo_out
                    evoked.save(evo_out)  # save evoked data to disk
                    cov.save(cov_out)  # save evoked data to disk
                    
                    
                    for emotion in ['Neu','Ang']:

                        event_id = {emotion: trig_logic[stim_type][session][trial_type][emotion]}                        
                        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                                            baseline=(None, 0), reject=None, preload=True)

                        # Look at channels that caused dropped events, showing that the subject's
                        # blinks were likely to blame for most epochs being dropped
                        # epochs.drop_bad_epochs()
                        # epochs.plot_drop_log()
                        # epochs.plot()

                        evoked = epochs.average()  # average epochs and get an Evoked dataset.
                            
                        evo_out = evo_path + '/' + trial_type + '_' + stim_type + '_' + session + '_' + emotion + '-ave.fif'
                        print evo_out
                        evoked.save(evo_out)  # save evoked data to disk


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

