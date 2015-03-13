
import matplotlib
matplotlib.use('agg') # force non-interactive plotting
import numpy as np
import os, errno
import json
from copy import deepcopy
from operator import add # for stc reduction operation

machine_name = os.uname()[1].split('.')[0]

if 'isis' in machine_name:
    import sys
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
    import subprocess
    from access import Query
    from analysis_dict import Anadict

    db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
    ad=Anadict(db)
elif 'mba-cjb' in machine_name or 'horus' in machine_name:
    class local_Anadict():
        def __init__(self):
            self._scratch_folder = '/Users/cjb/tmp/VSC-scratch'

    ad = local_Anadict()

import mne
#try:
from mne.io import Raw
from mne import pick_types
from mne.preprocessing import ICA, read_ica

def split_events_by_trialtype(events):
    devsA, devsB = range(111,117), range(211,217)
    VS_eve = mne.pick_events(events, include=range(100,220))
    VS_eve = mne.merge_events(VS_eve, [100], 10, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, [200], 20, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, devsA, 11, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, devsB, 21, replace_events=True)
    FB_eve = mne.pick_events(events, include=range(10,22))

    eve_dict = dict(VS=VS_eve, FB=FB_eve)
    id_dict = dict(stdA=10, stdB=20, devA=11, devB=21) # now same for VS and FB

    return eve_dict, id_dict

# Set epoch parameters
tmin, tmax = -0.4, 0.6  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.2, 0.2  # reject trial only if blinks in the 400 ms middle portion!
baseline = (-0.2, 0.)
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13) # compare to standard rejection

raw_path = ad._scratch_folder + '/tsss_initial/007_SGF'
eve_path = ad._scratch_folder + '/events.fif/007_SGF/raw'

fname = raw_path + '/VS_1a_1_tsss_mc.fif'

raw = Raw(fname, preload=True)

ica = read_ica(raw_path + '/ica_pre.fif')
print 'Excluding', ica.exclude
raw_ica = ica.apply(raw, copy=True)

events = mne.read_events(eve_path + '/VS_1a_1-eve.fif')
picks = pick_types(raw.info, meg=True, eeg=False, stim=True, eog=True, misc=True)
eve_dict, id_dict = split_events_by_trialtype(events)
for trial_type in ['VS']:

    epochs = mne.Epochs(raw, eve_dict[trial_type], id_dict,
                        tmin, tmax, picks=picks, verbose=False,
                        baseline=baseline, reject=reject, preload=True,
                        reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings
    epochs_ica = mne.Epochs(raw_ica, eve_dict[trial_type], id_dict,
                        tmin, tmax, picks=picks, verbose=False,
                        baseline=baseline, reject=None, preload=True)
    # Check the drop_log for these preload'ed epochs: does the drop
    # log indicate the dropped epochs, can they be un-dropped after the fact?
    # Do we in fact have to actively drop them, despite preloading?
    evo= epochs.average()
    evo_ica = epochs_ica.average()

    fig = evo.plot()
    fig.savefig(img_folder + '/evo.png')
    fig = evo_ica.plot()
    fig.savefig(img_folder + '/evo_ica.png')
