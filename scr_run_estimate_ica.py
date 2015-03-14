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
outdir_base = ad._scratch_folder + '/' + input_files + '/ica'
run_cond = ['VS', 'FFA']
    
for subj in ad.analysis_dict.keys():

    outdir = outdir_base + '/' + subj
    img_folder = outdir + '/img'
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
                # 1) Fit ICA model using the FastICA algorithm

                # Other available choices are `infomax` or `extended-infomax`
                # We pass a float value between 0 and 1 to select n_components based on the
                # percentage of variance explained by the PCA components.

                raw = Raw(fname, preload=True)
                picks = mne.pick_types(raw.info, meg=True, eeg=False, 
                        eog=False, ecg=False, stim=False, exclude='bads')

                if rank_estimate is None:
                    # estimate the rank only for the second VS task
                    # use 300 seconds
                    rank_estimate = raw.estimate_rank(tstart=240., tstop=540., picks=picks)
                    print 'Estimated raw to be of rank', rank_estimate

                ica = ICA(n_components=rank_estimate, max_pca_components = None, 
                        max_iter=256, method='fastica')

                ica.fit(raw, picks=picks, decim = 5, reject=dict(mag=4e-11, grad=4000e-12))
                # Save with information on excludes!
                ica.save(outdir + '/' + cond + '-ica.fif')
                

