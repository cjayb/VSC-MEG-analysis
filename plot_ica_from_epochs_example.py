# This  works poorly, IC's should definitely come from the raw
# data. MNE can't combine EOG traces with MEG into a single
# ICA. First, compute_covariance doesn't calculate non-meeg-
# channels. Second, the Z-scoring/PCA-whitening crashes if
# EOGs part of the data matrix (not sure why)


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
from mne import pick_types, compute_covariance
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
rej_tmin, rej_tmax = -0.2, 0.3  # reject trial only if blinks in the 500 ms middle portion!
baseline = (-0.2, 0.)
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13) # compare to standard rejection
ica_reject = dict(eog=300-6, mag=5e-12, grad=5000e-13) # compare to standard rejection

raw_path = ad._scratch_folder + '/tsss_initial/007_SGF'
img_folder = ad._scratch_folder + '/tsss_initial/007_SGF/img/epochs'
eve_path = ad._scratch_folder + '/events.fif/007_SGF/raw'

fname = raw_path + '/VS_1a_1_tsss_mc.fif'

raw = Raw(fname, preload=True)
picks = pick_types(raw.info, meg=True, eog=True)
n_components = raw.estimate_rank(picks=picks)
print 'Estimated raw to be of rank', n_components
n_max_eog = 3
title = 'Sources related to %s artifacts (red)'

#ica = read_ica(raw_path + '/ica_pre.fif')
#print 'Excluding', ica.exclude
#raw_ica = ica.apply(raw, copy=True)

events = mne.read_events(eve_path + '/VS_1a_1-eve.fif')
eve_dict, id_dict = split_events_by_trialtype(events)
for trial_type in ['VS']:

    epochs = mne.Epochs(raw, eve_dict[trial_type], id_dict,
                        tmin, tmax, picks=picks, verbose=False,
                        baseline=baseline, reject=reject, preload=True,
                        reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings

    picks = pick_types(epochs.info, meg=True, eog=False)

    # cannot compute EOG covariance?!?!!!
    #baseline_cov = compute_covariance(epochs, tmin=None, tmax=0)
    #ica = ICA(n_components=n_components, max_pca_components=None,max_iter=256, noise_cov=baseline_cov)
    ica = ICA(n_components=n_components, max_pca_components=None,max_iter=256)
    ica.fit(epochs, picks=picks, decim = 5, reject=ica_reject)

    eog_inds, scores = ica.find_bads_eog(epochs, ch_name='EOG001,EOG003')
    allscores = np.vstack((scores[0], scores[1]))
    mscores = np.max(np.abs(allscores), axis=0)
    # now scores is 
    show_picks = mscores.argsort()[::-1][:5]
    eog_inds = list(show_picks[:n_max_eog])
    ica.exclude += eog_inds

    fig = ica.plot_scores(scores, exclude=eog_inds, title=title % 'eog')
    fig.savefig(img_folder + '/ica_eog_scores.png')

    fig = ica.plot_sources(epochs, show_picks, exclude=eog_inds, title=title % 'eog')
    fig.savefig(img_folder + '/ica_eog_sources.png')
    fig = ica.plot_components(show_picks, title=title % 'eog', colorbar=True)
    fig.set_size_inches(12.,8.)
    fig.savefig(img_folder + '/ica_eog_components.png')

    # Now exclude what we find
    # epochs_ica = ica.apply(epochs, copy=True)
    #evo_ica = epochs_ica.average()
    evo= epochs.average()

    #megp = pick_types(epochs.info, meg=True, eog=False)
    #eogp = pick_types(epochs.info, meg=False, eog=True)
    #fig = ica.plot_overlay(evo, picks=megp)  # 
    fig = ica.plot_overlay(evo)  # 
    fig.savefig(img_folder + '/ica_evo_overlay_meg.png')
    #fig = ica.plot_overlay(evo, picks=eogp)  # 
    #fig.savefig(img_folder + '/ica_evo_overlay_eog.png')

    #fig = evo.plot()
    #fig.savefig(img_folder + '/epevo.png')
    #fig = evo_ica.plot()
    #fig.savefig(img_folder + '/epevo_ica.png')
