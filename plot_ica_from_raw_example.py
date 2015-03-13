"""
==================================
Compute ICA components on raw data
==================================

ICA is fit to MEG raw data.
The sources matching the ECG and EOG are automatically found and displayed.
Subsequently, artifact detection and rejection quality are assessed.
"""
print(__doc__)

# Authors: Denis Engemann <denis.engemann@gmail.com>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)

import sys
sys.path.append('/projects/MINDLAB2014_MEG-PTSD/misc/stormdb')

import numpy as np
import mne
from mne.io import Raw
from mne.preprocessing import ICA, read_ica
from mne.preprocessing import create_ecg_epochs, create_eog_epochs
from mne.cov import compute_raw_data_covariance
from stormdb.access import Query

import matplotlib as mpl
mpl.use('Agg')

def mkdir_p(pth):                                                                
    import os, errno                                                                                 
    try:                                                                         
        os.makedirs(pth)                                                         
    except OSError as exc:                                                       
        if exc.errno == errno.EEXIST and os.path.isdir(pth):                     
            pass                                                                 
        else:                                                                    
            raise                                                                

###############################################################################
# Setup paths and prepare raw data

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

VERBOSE=True
SAVE=False # NB

db = Query(proj_code=proj_code,verbose=True)
subjects = db.get_subjects()
example_subject = subjects[1][1:] #007
print "Using {:s} as example".format(example_subject)

proj_folder = '/projects/' + proj_code
misc_folder = proj_folder + '/misc'
scratch_folder = proj_folder + '/scratch'
raw_data_folder = scratch_folder + '/tsss_initial/' + \
        example_subject
img_folder = raw_data_folder + '/img'
mkdir_p(img_folder)

raw_fname = raw_data_folder + '/VS_1a_1_tsss_mc.fif'

###############################################################################
# 0) calculate noise covariance from empty room data
#print "calculate noise covariance from empty room data"
#empty_room = Raw(empty_room_fname, preload=True)
#picks = mne.pick_types(empty_room.info, meg=True, eeg=False, eog=False,
#                       stim=False, exclude='bads')
#noise_cov_er = compute_raw_data_covariance(empty_room, tmin=0., tmax=120.,
#        reject=dict(mag=4e-12, grad=4000e-13), picks=picks)
#noise_cov_er.save(empty_room_fname[:-4]+'-cov.fif')
###############################################################################
# 1) Fit ICA model using the FastICA algorithm

# Other available choices are `infomax` or `extended-infomax`
# We pass a float value between 0 and 1 to select n_components based on the
# percentage of variance explained by the PCA components.

raw = Raw(raw_fname, preload=True)

picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=False, ecg=False,
                       stim=False, exclude='bads')
# maximum number of components to reject
n_max_ecg, n_max_eog = 2, 3  # here we expect horizontal EOG components
try:
    ica = read_ica(raw_data_folder + '/ica_pre.fif')
except:
    ica = ICA(n_components=0.95, max_pca_components = 64, method='fastica')
        #noise_cov = noise_cov_er)

    #ica.fit(raw, picks=picks, decim=3, reject=dict(mag=4e-12, grad=4000e-13))
    ica.fit(raw, picks=picks, decim = 5, reject=dict(mag=4e-11, grad=4000e-12))
    # To save an ICA solution you can say:

###############################################################################
# 2) identify bad components by analyzing latent sources.

title = 'Sources related to %s artifacts (red)'

# generate ECG epochs use detection via phase statistics

picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, ecg=True,
                       stim=False, exclude='bads')

# create_ecg_epochs is strange: it strips the channels of anything non M/EEG
# UNLESS picks=None
picks=None
ecg_epochs = create_ecg_epochs(raw, ch_name='ECG002', tmin=-.5, tmax=.5, picks=picks, verbose=True)

# This will work with the above, but uses MASSIVE RAM
# Not sure the ECG quality is good enough for the QRS-detector
ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps', ch_name='ECG002')
# This creates a synthetic ECG from magnetometers, probably better...
ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps', ch_name='ECG002')

fig = ica.plot_scores(scores, exclude=ecg_inds, title=title % 'ecg')
fig.savefig(img_folder + '/ica_ecg_scores.png')

show_picks = np.abs(scores).argsort()[::-1][:5]

fig = ica.plot_sources(raw, show_picks, exclude=ecg_inds, title=title % 'ecg')
fig.savefig(img_folder + '/ica_ecg_sources.png')
fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
fig.set_size_inches(12.,8.)
fig.savefig(img_folder + '/ica_ecg_components.png')

ecg_inds = ecg_inds[:n_max_ecg]
ica.exclude += ecg_inds

# detect EOG by correlation
# First do both together
eog_inds, scores = ica.find_bads_eog(raw, ch_name='EOG001,EOG003')
allscores = np.vstack((scores[0], scores[1]))
mscores = np.max(np.abs(allscores), axis=0)
# now scores is 
show_picks = mscores.argsort()[::-1][:5]
eog_inds = list(show_picks[:n_max_eog])

fig = ica.plot_scores(scores, exclude=eog_inds, title=title % 'eog')
fig.savefig(img_folder + '/ica_eog_scores.png')

fig = ica.plot_sources(raw, show_picks, exclude=eog_inds, title=title % 'eog')
fig.savefig(img_folder + '/ica_eog_sources.png')
fig = ica.plot_components(show_picks, title=title % 'eog', colorbar=True)
fig.set_size_inches(12.,8.)
fig.savefig(img_folder + '/ica_eog_components.png')

# NB: Not using the "combined" EOG components, but the separate ones from below

# Then VEOG
veog_inds, scores = ica.find_bads_eog(raw, ch_name='EOG001')
show_picks = np.abs(scores).argsort()[::-1][:5]
veog_inds = list(show_picks[:n_max_eog])

fig = ica.plot_scores(scores, exclude=veog_inds, title=title % 'veog')
fig.savefig(img_folder + '/ica_veog_scores.png')

fig = ica.plot_sources(raw, show_picks, exclude=veog_inds, title=title % 'veog')
fig.savefig(img_folder + '/ica_veog_sources.png')
fig = ica.plot_components(show_picks, title=title % 'veog', colorbar=True)
fig.set_size_inches(12.,8.)
fig.savefig(img_folder + '/ica_veog_components.png')

#eog_inds = eog_inds[:n_max_eog]
veog_inds = veog_inds[:1] # manually determined, take only first IC
ica.exclude += veog_inds

# Then HEOG
heog_inds, hscores = ica.find_bads_eog(raw, ch_name='EOG003', l_freq=5., h_freq=20.)
show_picks_h = np.abs(hscores).argsort()[::-1][:5]
heog_inds = list(show_picks_h[:n_max_eog])

fig = ica.plot_scores(hscores, exclude=heog_inds, title=title % 'heog')
fig.savefig(img_folder + '/ica_heog_scores.png')

fig = ica.plot_sources(raw, show_picks_h, exclude=heog_inds, title=title % 'heog')
fig.savefig(img_folder + '/ica_heog_sources.png')
fig = ica.plot_components(show_picks_h, title=title % 'heog', colorbar=True)
fig.set_size_inches(12.,8.)
fig.savefig(img_folder + '/ica_heog_components.png')

heog_inds = heog_inds[:2] # manually determined, take only first two ICs
ica.exclude += heog_inds

###############################################################################
# 3) Assess component selection and unmixing quality

# estimate average artifact
ecg_evoked = ecg_epochs.average()
fig = ica.plot_sources(ecg_evoked, exclude=ecg_inds)  # plot ECG sources + selection
fig.savefig(img_folder + '/ica_ecg_evoked_sources.png')
fig = ica.plot_overlay(ecg_evoked, exclude=ecg_inds)  # plot ECG cleaning
fig.savefig(img_folder + '/ica_ecg_evoked_overlay.png')

#eog_evoked = create_eog_epochs(raw, tmin=-.5, tmax=.5, picks=picks).average()
#fig = ica.plot_sources(eog_evoked, exclude=eog_inds)  # plot EOG sources + selection
#fig.savefig(img_folder + '/ica_eog_evoked_sources.png')
#fig = ica.plot_overlay(eog_evoked, exclude=eog_inds)  # plot EOG cleaning
#fig.savefig(img_folder + '/ica_eog_evoked_overlay.png')

tmp=ica.exclude
ica.exclude = []
veog_evoked = create_eog_epochs(raw, ch_name='EOG001', tmin=-.5, tmax=.5, picks=picks).average()
fig = ica.plot_sources(veog_evoked, exclude=veog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_veog_evoked_sources_veog_inds.png')
fig = ica.plot_overlay(veog_evoked, exclude=veog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_veog_evoked_overlay_veog_inds.png')
fig = ica.plot_sources(veog_evoked, exclude=heog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_veog_evoked_sources_heog_inds.png')
fig = ica.plot_overlay(veog_evoked, exclude=heog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_veog_evoked_overlay_heog_inds.png')
fig = ica.plot_sources(veog_evoked, exclude=eog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_veog_evoked_sources_eog_inds.png')
fig = ica.plot_overlay(veog_evoked, exclude=eog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_veog_evoked_overlay_eog_inds.png')

heog_evoked = create_eog_epochs(raw, ch_name='EOG003', tmin=-.5, tmax=.5, picks=picks).average()
fig = ica.plot_sources(heog_evoked, exclude=veog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_heog_evoked_sources_veog_inds.png')
fig = ica.plot_overlay(heog_evoked, exclude=veog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_heog_evoked_overlay_veog_inds.png')
fig = ica.plot_sources(heog_evoked, exclude=heog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_heog_evoked_sources_heog_inds.png')
fig = ica.plot_overlay(heog_evoked, exclude=heog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_heog_evoked_overlay_heog_inds.png')
fig = ica.plot_sources(heog_evoked, exclude=eog_inds)  # plot EOG sources + selection
fig.savefig(img_folder + '/ica_heog_evoked_sources_eog_inds.png')
fig = ica.plot_overlay(heog_evoked, exclude=eog_inds)  # plot EOG cleaning
fig.savefig(img_folder + '/ica_heog_evoked_overlay_eog_inds.png')

ica.exclude=tmp

# check the amplitudes do not change
fig = ica.plot_overlay(raw, start=611.5, stop=612.)  # 
fig.savefig(img_folder + '/ica_raw_overlay.png')

###############################################################################
# Save with information on excludes!
ica.save(raw_data_folder + '/ica_pre.fif')
#
# You can later load the solution by saying:
# >>> from mne.preprocessing import read_ica
# >>> read_ica('my_ica.fif')
#
# Apply the solution to Raw, Epochs or Evoked like this:
#ica.apply(raw, copy=False)
#raw.savefig(raw_data_folder + '/ec_rest_before_tsss_mc_rsl_ica.fif')
