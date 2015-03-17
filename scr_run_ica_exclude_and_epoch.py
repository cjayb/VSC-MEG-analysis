import matplotlib
#matplotlib.use('Agg')
matplotlib.use('Qt')

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

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def split_events_by_trialtype(events):
    devsA, devsB = range(111,117), range(211,217)
    VS_eve = mne.pick_events(events, include=range(100,220))
    VS_eve = mne.merge_events(VS_eve, [100], 10, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, [200], 20, replace_events=True)
    # Don't replace the deviants, make a copy instead!
    VS_eve = mne.merge_events(VS_eve, devsA, 11, replace_events=False)
    VS_eve = mne.merge_events(VS_eve, devsB, 21, replace_events=False)
    FB_eve = mne.pick_events(events, include=range(10,22))
    
    FFA_eve = mne.pick_events(events, include=[100, 150, 200])

    eve_dict = dict(VS=VS_eve, FB=FB_eve, FFA=FFA_eve)
    id_dict = dict(stdA=10, stdB=20, devA=11, devB=21,
            A1=111, A2=112,A3=113,A4=114,A5=115,A6=116,
            B1=211, B2=212,B3=213,B4=214,B5=215,B6=216,
            A=100, B=200, blur=150) 

    return eve_dict, id_dict

def load_exludes(ica_excludes_folder, subj, condition):
    import csv
    pth = ica_excludes_folder + '/' + subj + '.csv'

    with open(pth, 'rb') as csvfile:
        exreader = csv.reader(csvfile, delimiter=',')
        hdr = exreader.next()
        try:
            colind = hdr.index(condition)
        except ValueError:
            print 'condition must be VS1, VS2 or FFA!'
            raise ValueError

        ica_excludes = []
        for row in exreader:
            ica_excludes += row[colind].split('|')

    return ica_excludes

input_files = 'tsss_initial'
ica_estimates = ad._scratch_folder + '/ica/' + input_files 
corr_eves = ad._scratch_folder + '/events.fif/'
ica_excludes_folder = ad._misc_folder + '/ica/' + input_files
t_start, t_stop = 40.0, 100. # for displaying sources of 1 minutes, works even for short FFA-sessions..
# Set epoch parameters
tmin, tmax = -0.2, 0.4  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.1, 0.2  # reject trial only if blinks in the 300 ms middle portion!
#reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13) # compare to standard rejection
reject = None
baseline = (-0.1, 0.)
    
#for subj in ad.analysis_dict.keys():
for subj in ['006_HEN',]:

    eve_folder = corr_eves + subj + '/raw'
    ica_folder = ica_estimates + '/' + subj
    img_folder = ica_folder + '/img'

    cond_names = ad.analysis_dict[subj][input_files].keys()
    # sort names so that VS comes before FFA!
    cond_names.sort(reverse=True)
    for cond in cond_names:
        img_prefix = img_folder + '/' + cond
        if 'empty' not in cond:

            ica_excludes = load_exludes(ica_excludes_folder, subj, cond)
            if 'VS' in cond:
                trial_types = ['VS','FB']
                ica_check_eves = ['stdA','stdB']
            elif 'FFA' in cond:
                trial_types = ['FFA',]
                ica_check_eves = ['A','B']
            
            raw_path = ad._scratch_folder + '/' + input_files + '/' + subj
            in_fnames = ad.analysis_dict[subj][input_files][cond]['files'] 
            events = mne.read_events(eve_folder + '/' + cond + '-eve.fif')
            eve_dict, id_dict = split_events_by_trialtype(events)
            for fname in in_fnames:

                print 'In: ', fname
                raw = Raw(fname, preload=False) 
                ica = read_ica(ica_folder + '/' + cond + '-ica.fif')

                picks = mne.pick_types(raw.info, meg=True)

                for trial_type in trial_types:
                    epochs = mne.Epochs(raw, eve_dict[trial_type], id_dict,
                                    tmin, tmax, picks=picks, verbose=False,
                                    baseline=baseline, reject=reject,
                                    preload=False, reject_tmin=rej_tmin,
                                    reject_tmax=rej_tmax) # Check rejection settings
                    ica_check_evoked = epochs[ica_check_eves].average()

                fig = ica.plot_overlay(ica_check_evoked, exclude=ica_excludes)  # plot EOG cleaning
                #fig.savefig(img_prefix + '_eog_evoked_overlay_heog.png')
