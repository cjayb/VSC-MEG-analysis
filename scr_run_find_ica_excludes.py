import sys
sys.path.append('/projects/MINDLAB2014_MEG-PTSD/misc/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

import mne
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
t_start, t_stop = 240.0, 540. # for displaying sources of 5 minutes
    
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
                print 'In: ', fname
                raw = Raw(fname, preload=False)
                ica = read_ica(ica_folder + '/' + cond + '-ica.fif')
                # 2) identify bad components by analyzing latent sources.
                title = 'Sources related to %s artifacts (red)'

                # generate ECG epochs use detection via phase statistics

                picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, ecg=True, stim=False, exclude='bads')

                # create_ecg_epochs is strange: it strips the channels of anything non M/EEG
                # UNLESS picks=None
                picks=None
                ecg_epochs = create_ecg_epochs(raw, ch_name='ECG002', tmin=-.2, tmax=.5, picks=picks, flat='grad', verbose=True)

                # This will work with the above, but uses MASSIVE RAM
                # Not sure the ECG quality is good enough for the QRS-detector
                ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps', ch_name='ECG002')

                fig = ica.plot_scores(scores, exclude=ecg_inds, title=title % 'ecg')
                fig.savefig(img_folder + '/ica_ecg_scores.png')

                show_picks = np.abs(scores).argsort()[::-1][:5]

                fig = ica.plot_sources(raw, show_picks, exclude=ecg_inds, start=t_start, stop=t_stop, title=title % 'ecg')
                fig.savefig(img_folder + '/ica_ecg_sources.png')
                fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
                fig.set_size_inches(12.,8.)
                fig.savefig(img_folder + '/ica_ecg_components.png')
