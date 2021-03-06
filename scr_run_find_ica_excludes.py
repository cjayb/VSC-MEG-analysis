import matplotlib
matplotlib.use('Agg')
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
from mne.preprocessing import ICA, read_ica
from mne.preprocessing import create_ecg_epochs, create_eog_epochs

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

input_files = 'tsss_initial'
ica_estimates = ad._scratch_folder + '/ica/' + input_files 
run_cond = ['VS', 'FFA']
t_start, t_stop = 40.0, 100. # for displaying sources of 1 minutes, works even for short FFA-sessions..
n_max_ecg, n_max_eog = 3, 4
    
for subj in ad.analysis_dict.keys():

    ica_folder = ica_estimates + '/' + subj
    img_folder = ica_folder + '/img'
    mkdir_p(img_folder)

    # Reset for each subject
    rank_estimate = None

    cond_names = ad.analysis_dict[subj][input_files].keys()
    # sort names so that VS comes before FFA!
    cond_names.sort(reverse=True)
    for cond in cond_names:
        if 'empty' not in cond:
            
            raw_path = ad._scratch_folder + '/' + input_files + '/' + subj
            in_fnames = ad.analysis_dict[subj][input_files][cond]['files'] 
            for fname in in_fnames:

                img_prefix = img_folder + '/' + cond

                print 'In: ', fname
                raw = Raw(fname, preload=True) # for finding events from raw, must be preloaded
                ica = read_ica(ica_folder + '/' + cond + '-ica.fif')
                # 2) identify bad components by analyzing latent sources.
                title = 'Sources related to %s artifacts (red)'

                # generate ECG epochs use detection via phase statistics

                picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, ecg=True, stim=False, exclude='bads')

                # create_ecg_epochs is strange: it strips the channels of anything non M/EEG
                # UNLESS picks=None
                #picks=None
                # This will work with the above, but uses MASSIVE RAM
                # Not sure the ECG quality is good enough for the QRS-detector
                ecg_inds, scores = ica.find_bads_ecg(raw, method='ctps', ch_name='ECG002', threshold=0.25)
                if len(ecg_inds) < 1:
                    # destroy the ECG channel by changing it to an EMG!
                    raw.info['chs'][1]['kind'] = 302
                    # then try to generate ECG from mags
                    ecg_inds, scores = ica.find_bads_ecg(raw, method='ctps', ch_name='ECG002', threshold=0.25)

                if len(ecg_inds) > 0:
                    fig = ica.plot_scores(scores, exclude=ecg_inds, title=title % 'ecg')
                    fig.savefig(img_prefix + '_ecg_scores.png')

                    show_picks = np.abs(scores).argsort()[::-1][:n_max_ecg]

                    fig = ica.plot_sources(raw, show_picks, exclude=ecg_inds, start=t_start, stop=t_stop, title=title % 'ecg')
                    fig.savefig(img_prefix + '_ecg_sources.png')
                    fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
                    fig.set_size_inches(9.,6.)
                    fig.savefig(img_prefix + '_ecg_components.png')

                    # estimate average artifact
                    ecg_epochs = create_ecg_epochs(raw, ch_name='ECG002', tmin=-.2, tmax=.5, picks=picks, flat='grad', verbose=False)

                    ecg_evoked = ecg_epochs.average()
                    fig = ica.plot_sources(ecg_evoked, exclude=ecg_inds)  # plot ECG sources + selection
                    fig.savefig(img_prefix + '_ecg_evoked_sources.png')
                    fig = ica.plot_overlay(ecg_evoked, exclude=ecg_inds)  # plot ECG cleaning
                    fig.savefig(img_prefix + '_ecg_evoked_overlay.png')


                # detect EOG by correlation
                # do both together
                eog_inds, scores = ica.find_bads_eog(raw, ch_name='EOG001,EOG003')
                allscores = np.vstack((scores[0], scores[1]))
                mscores = np.max(np.abs(allscores), axis=0)
                # now scores is 
                show_picks = mscores.argsort()[::-1][:n_max_eog]
                eog_inds = list(show_picks[:n_max_eog])

                fig = ica.plot_scores(scores, exclude=eog_inds, title=title % 'eog')
                fig.savefig(img_prefix + '_eog_scores.png')

                fig = ica.plot_sources(raw, show_picks, exclude=eog_inds, title=title % 'eog')
                fig.savefig(img_prefix + '_eog_sources.png')
                fig = ica.plot_components(show_picks, title=title % 'eog', colorbar=True)
                fig.set_size_inches(9.,6.)
                fig.savefig(img_prefix + '_eog_components.png')

                # estimate average artifact
                veog_evoked = create_eog_epochs(raw, ch_name='EOG001', tmin=-.2, tmax=.5, picks=picks).average()
                fig = ica.plot_sources(veog_evoked, exclude=eog_inds)  # plot EOG sources + selection
                fig.savefig(img_prefix + '_eog_evoked_sources_veog.png')
                fig = ica.plot_overlay(veog_evoked, exclude=eog_inds)  # plot EOG cleaning
                fig.savefig(img_prefix + '_eog_evoked_overlay_veog.png')

                heog_evoked = create_eog_epochs(raw, ch_name='EOG003', tmin=-.5, tmax=.5, picks=picks).average()
                fig = ica.plot_sources(heog_evoked, exclude=eog_inds)  # plot EOG sources + selection
                fig.savefig(img_prefix + '_eog_evoked_sources_heog.png')
                fig = ica.plot_overlay(heog_evoked, exclude=eog_inds)  # plot EOG cleaning
                fig.savefig(img_prefix + '_eog_evoked_overlay_heog.png')
