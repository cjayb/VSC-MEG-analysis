import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('/projects/MINDLAB2014_MEG-PTSD/misc/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

import mne
import numpy as np
from mne.io import Raw
from mne.preprocessing import read_ica
from mne.report import Report

from VSC_utils import *

#input_files = 'tsss_initial'
#filter_string = '0.5-35.0Hz'
filter_string = filt_dir # from VSC_utils!

ica_estimates_path = ad._scratch_folder + '/ica/' + input_files 
ica_epochs_path = ad._scratch_folder + '/epochs/ica/' + input_files 
corr_eves_path = ad._scratch_folder + '/events.fif/'
ica_excludes_path = ad._misc_folder + '/ica/' + input_files

report_folder = ica_epochs_path + '/report'
mkdir_p(report_folder) # parents will be made too

# for checking the evokeds, use savgol on the unfiltered raw
savgol_hf = 35.

CLOBBER = False

for subj in ad.analysis_dict.keys():
#for subj in ['006_HEN',]:

    eve_folder = corr_eves_path + subj + '/raw'
    ica_folder = ica_estimates_path + '/' + subj

    epochs_folder = ica_epochs_path + '/' + subj
    epochs_folder_filt = epochs_folder + '/' + filter_string

    #ica_check_img_folder = epochs_folder + '/img'

    if not CLOBBER and os.path.isdir(epochs_folder):
        continue
    
    mkdir_p(epochs_folder_filt) # parents will be made too
    #mkdir_p(ica_check_img_folder)

    report = Report(info_fname=None, subjects_dir=None, subject=subj,
                    title='Epoching check with ICA applied', verbose=None)

    
    cond_names = ad.analysis_dict[subj][input_files].keys()
    # sort names alphabetically
    cond_names.sort()

    for cond in cond_names:
        if 'empty' not in cond:
            if 'VS' in cond:
                trial_types = ['VS','FB']
                session_no = cond[-1]
                ica_check_eves = ['stdA','stdB']
                ica_cond = 'VS' + session_no
            elif 'FFA' in cond:
                trial_types = ['FFA',]
                session_no = ''
                ica_check_eves = ['A','B']
                ica_cond = 'FFA'
            
            ica_excludes = load_excludes(ica_excludes_path, subj, ica_cond)
            print 30*'*'
            print 'ICA excludes:', ica_excludes
            print 30*'*'
            raw_path = ad._scratch_folder + '/' + input_files + '/' + subj
            filtered_path = ad._scratch_folder + '/filtered/' + input_files + \
                '/' + filter_string + '/' + subj
            in_fnames = ad.analysis_dict[subj][input_files][cond]['files'] 
            events = mne.read_events(eve_folder + '/' + cond + '-eve.fif')
            eve_dict, id_dict = split_events_by_trialtype(events,condition=cond)
            for fname in in_fnames:

                print 'In: ', fname
                raw = Raw(fname, preload=False) 
                raw_filt = Raw(filtered_path + '/' +cond+'_filt.fif', preload=False) 
                ica = read_ica(ica_folder + '/' + cond + '-ica.fif')

                picks = mne.pick_types(raw.info, meg=True, eog=True)

                for trial_type in trial_types:
                    epochs = mne.Epochs(raw, eve_dict[trial_type], id_dict[trial_type],
                                    tmin, tmax, picks=picks, verbose=False,
                                    baseline=baseline, reject=reject,
                                    preload=True, reject_tmin=rej_tmin,
                                    reject_tmax=rej_tmax) # Check rejection settings
                    ica_check_evoked = epochs[ica_check_eves].average()

                    filtered_epochs = mne.Epochs(raw_filt, eve_dict[trial_type], id_dict[trial_type],
                                    tmin, tmax, picks=picks, verbose=False,
                                    baseline=baseline, reject=reject,
                                    preload=True, reject_tmin=rej_tmin,
                                    reject_tmax=rej_tmax) # Check rejection settings

                    ica_check_evoked = epochs[ica_check_eves].average()
                    ica_check_evoked.savgol_filter(h_freq=savgol_hf)

                    ica_check_evoked_filt = filtered_epochs[ica_check_eves].average()

                    fig = ica.plot_overlay(ica_check_evoked, exclude=ica_excludes)  # plot EOG cleaning
                    #fig.savefig(ica_check_img_folder + '/' +trial_type + session_no + '-savgol.png')
                    report.add_figs_to_section(fig, trial_type + session_no, 
                            section='Savitzky-Golay', scale=None, image_format='png')
                    plt.close(fig)
                    fig = ica.plot_overlay(ica_check_evoked_filt, exclude=ica_excludes)  # plot EOG cleaning
                    #fig.savefig(ica_check_img_folder + '/' +trial_type + session_no + '-rawfilt.png')
                    report.add_figs_to_section(fig, trial_type + session_no, 
                            section='Filtered from raw', scale=None, image_format='png')
                    plt.close(fig)

                    ica.exclude = ica_excludes

                    ica.apply(epochs, copy=False)
                    ica.apply(filtered_epochs, copy=False)

                    print('Resampling epochs...')
                    epochs.resample(epoch_params['rsl'], n_jobs=4, verbose=False) # Trust the defaults here
                    print('Resampling filtered_epochs...')
                    filtered_epochs.resample(epoch_params['rsl'], n_jobs=4, verbose=False) # Trust the defaults here

                    epochs.save(epochs_folder + '/' + trial_type + session_no  + '-epo.fif')
                    filtered_epochs.save(epochs_folder_filt + '/' + trial_type + session_no  + '_filt-epo.fif')

    report.save(fname=report_folder + '/' + subj + '.html', open_browser = False, overwrite = CLOBBER)
