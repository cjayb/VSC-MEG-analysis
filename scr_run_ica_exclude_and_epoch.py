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
import os, errno
from mne.io import Raw
from mne.preprocessing import read_ica
import csv

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def split_events_by_trialtype(events, condition='VS'):
    if 'VS' in condition:
        devsA, devsB = range(111,117), range(211,217)
        VS_eve = mne.pick_events(events, include=range(100,220))
        VS_eve = mne.merge_events(VS_eve, [100], 10, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, [200], 20, replace_events=True)
        # Don't replace the deviants, make a copy instead!
        VS_eve = mne.merge_events(VS_eve, devsA, 11, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, devsB, 21, replace_events=True)

        # This hack is needed to get both 11/21's and 11N/21N's together!
        tmp = mne.pick_events(events, include=devsA+devsB)
        #tmp[:,0] += 1 # add a ms
        VS_eve = np.concatenate((VS_eve, tmp), axis=0)
        VS_eve = VS_eve[np.argsort(VS_eve[:, 0])]

        FB_eve = mne.pick_events(events, include=range(10,22))
        
        eve_dict = dict(VS=VS_eve, FB=FB_eve)

    elif 'FFA' in condition:
        FFA_eve = mne.pick_events(events, include=[100, 150, 200])
        eve_dict = dict(FFA=FFA_eve)

    id_dict = dict(VS=dict(stdA=10, stdB=20, devA=11, devB=21,
            A1=111, A2=112,A3=113,A4=114,A5=115,A6=116,
            B1=211, B2=212,B3=213,B4=214,B5=215,B6=216),
            FB=dict(stdA=10, stdB=20, devA=11, devB=21),
            FFA=dict(A=100, B=200, blur=150))

    return eve_dict, id_dict

def load_exludes(ica_excludes_folder, subj, cond):
    pth = ica_excludes_folder + '/' + subj + '.csv'

    with open(pth, 'rb') as csvfile:
        exreader = csv.reader(csvfile, delimiter=',')
        hdr = exreader.next()
        try:
            colind = hdr.index(cond)
        except ValueError:
            print 'condition must be VS1, VS2 or FFA!'
            raise ValueError

        ica_excludes = []
        for row in exreader:
            ica_excludes += row[colind].split('|')

    return map(int,ica_excludes)

input_files = 'tsss_initial'
filter_string = '0.5-35.0Hz'
ica_estimates_path = ad._scratch_folder + '/ica/' + input_files 
ica_epochs_path = ad._scratch_folder + '/epochs/ica/' + input_files 
corr_eves_path = ad._scratch_folder + '/events.fif/'
ica_excludes_path = ad._misc_folder + '/ica/' + input_files

# Set epoch parameters
tmin, tmax = -0.2, 0.4  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.15, 0.2  # reject trial only if blinks in the 300 ms middle portion!
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13) # compare to standard rejection
#reject = None
baseline = (-0.15, 0.)
# for checking the evokeds, use savgol on the unfiltered raw
savgol_hf = 35.
    
for subj in ad.analysis_dict.keys():
#for subj in ['006_HEN',]:

    eve_folder = corr_eves_path + subj + '/raw'
    ica_folder = ica_estimates_path + '/' + subj

    epochs_folder = ica_epochs_path + '/' + subj
    epochs_folder_filt = epochs_folder + '/' + filter_string
    ica_check_img_folder = epochs_folder + '/img'

    mkdir_p(epochs_folder_filt) # parents will be made too
    mkdir_p(ica_check_img_folder)
    
    cond_names = ad.analysis_dict[subj][input_files].keys()
    # sort names alphabetically
    cond_names.sort()

    for cond in cond_names:
        if 'empty' not in cond:
            if 'VS' in cond:
                trial_types = ['VS','FB']
                ica_check_eves = ['stdA','stdB']
                ica_cond = 'VS' + cond[-1]
            elif 'FFA' in cond:
                trial_types = ['FFA',]
                ica_check_eves = ['A','B']
                ica_cond = 'FFA'
            
            ica_excludes = load_exludes(ica_excludes_path, subj, ica_cond)
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
                    fig.savefig(ica_check_img_folder + '/' + cond + '-savgol.png')
                    plt.close(fig)
                    fig = ica.plot_overlay(ica_check_evoked_filt, exclude=ica_excludes)  # plot EOG cleaning
                    fig.savefig(ica_check_img_folder+ '/' + cond + '-rawfilt.png')
                    plt.close(fig)

                    ica.exclude = ica_excludes

                    ica.apply(epochs, copy=False)
                    ica.apply(filtered_epochs, copy=False)
                    epochs.save(epochs_folder + '/' + cond + '-epo.fif')
                    filtered_epochs.save(epochs_folder_filt + '/' + cond + '_filt-epo.fif')

