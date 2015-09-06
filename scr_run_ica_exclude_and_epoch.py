import mne
import numpy as np
from mne.io import Raw
from mne.preprocessing import read_ica
from mne.report import Report

from VSC_utils import *

import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from stormdb.access import Query
from analysis_dict import Anadict

db=Query(proj_code)
ad=Anadict(db)

performBandpassFilter = True
filter_string = filt_dir if performBandpassFilter else ''

ica_estimates_path = opj(scratch_folder, 'ica', input_files)
ica_epochs_path = opj(scratch_folder, 'epochs',
                               input_files, 'ica', filter_string)
corr_eves_path = opj(scratch_folder, 'events.fif')
ica_excludes_path = opj(misc_folder, 'ica', input_files)

report_folder = opj(ica_epochs_path, 'report')
mkdir_p(report_folder) # parents will be made too

CLOBBER = False

for subj in ad.analysis_dict.keys():
#for subj in ['006_HEN',]:

    eve_folder = opj(corr_eves_path, subj, 'raw')
    ica_folder = opj(ica_estimates_path, subj)

    epochs_folder = opj(ica_epochs_path, subj)

    #ica_check_img_folder = epochs_folder + '/img'

    if not CLOBBER and os.path.isdir(epochs_folder):
        continue

    mkdir_p(epochs_folder) # parents will be made too

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
                ica_check_eves = dict(VS=\
                    evoked_categories['VS']['stdA'][0] + \
                    evoked_categories['VS']['devA'][0],
                    FB=['stdA', 'stdB'])
                ica_cond = 'VS' + session_no
            elif 'FFA' in cond:
                trial_types = ['FFA',]
                session_no = ''
                ica_check_eves = dict(FFA=['A','B'])
                ica_cond = 'FFA'

            ica_excludes = load_excludes(ica_excludes_path, subj, ica_cond)
            print 30*'*'
            print 'ICA excludes:', ica_excludes
            print 30*'*'
            raw_path = opj(ad._scratch_folder, input_files, subj)
            in_fnames = ad.analysis_dict[subj][input_files][cond]['files']
            events = mne.read_events(opj(eve_folder, cond + '-eve.fif'))
            eve_dict, id_dict = \
                 split_events_by_trialtype(events, condition=cond)
            for fname in in_fnames:

                ica = read_ica(opj(ica_folder, cond + '-ica.fif'))

                print 'In: ', fname
                raw = Raw(fname, preload=performBandpassFilter)
                rep_section_name = ''
                if performBandpassFilter:
                    raw.filter(filter_params['highpass'],
                               filter_params['lowpass'],
                               method='iir', n_jobs=4
                               )

                picks = mne.pick_types(raw.info, meg=True, eog=True)

                for trial_type in trial_types:
                    epochs = mne.Epochs(raw, eve_dict[trial_type],
                                    id_dict[trial_type],
                                    tmin, tmax, picks=picks, verbose=False,
                                    baseline=baseline, reject=reject,
                                    preload=True, reject_tmin=rej_tmin,
                                    reject_tmax=rej_tmax) # Check rejection settings
                    ica_check_evoked = \
                                epochs[ica_check_eves[trial_type]].average()

                    fig = ica.plot_overlay(ica_check_evoked, exclude=ica_excludes)  # plot EOG cleaning
                    #fig.savefig(ica_check_img_folder + '/' +trial_type + session_no + '-savgol.png')
                    report.add_figs_to_section(fig, trial_type + session_no,
                            section=filter_string, scale=None, image_format='png')
                    plt.close(fig)

                    ica.exclude = ica_excludes

                    ica.apply(epochs, copy=False)

                    print('Resampling epochs...')
                    epochs.resample(epoch_params['rsl'],
                                    n_jobs=4, verbose=False)
                                    # Trust the defaults here

                    epochs.save(opj(epochs_folder,
                                    trial_type + session_no + '-epo.fif'))

    report.save(fname=report_folder + '/' + subj + '.html', open_browser = False, overwrite = CLOBBER)
