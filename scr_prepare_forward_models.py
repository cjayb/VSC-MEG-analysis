# -*- coding: utf-8 -*-
"""
Complete forward model by including head-MRI transformation

Input:
scratch/trans/SUB_ID-trans.fif
scratch/fs_subjects_dir/SUB_ID/bem/SUB_ID-head.fif
scratch/fs_subjects_dir/SUB_ID/bem/SUB_ID-5120-bem-sol.fif
scratch/fs_subjects_dir/SUB_ID/bem/SUB_ID-oct-6-src.fif

Output:
scratch/fwd/SUB_ID/CONDITION-fwd.fif

Created on 22 May 2014

@author: cjayb
"""

import sys
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

recon_all_bin = '/opt/local/freesurfer-releases/5.3.0/bin/recon-all'
subjects_dir = ad._scratch_folder + '/fs_subjects_dir'
# check_path_exists(subjects_dir)

VERBOSE=True

# Base it on this
# mne_do_forward_solution --spacing oct-6 --bem 030_WAH-5120-bem-sol.fif
# --trans ../../../trans/030_WAH-trans.fif --meas ../../../tsss_initial/030_WAH/FFA_tsss_mc.fif
# --destdir ../../../fwd/030_WAH --mindist 5 --megonly

# BUT: don't do it on the raw files, since we might set up some SSP or ICA projectors which have
# to be take into account in the FWD model!

params = {'spacing': ' --spacing oct-6 --megonly ',
        'bem': '-5120-bem-sol ',
        'force': False}


for subj in ad.analysis_dict.keys():

    bash_script = ['#!/usr/bin/env bash']
    bash_script.append('source ~/.bashrc')
    bash_script.append('use mne')
    bash_script.append('export SUBJECTS_DIR='+subjects_dir)

    if not any('recon-all' in item for item in ad.analysis_dict[subj].keys()):
        print "Skipping %s due to missing recon-all reconstruction" % subj
        continue

    bash_script.append('export SUBJECT=' + subj)
    fwd_cmd = 'mne_do_forward_solution'
    if params['force']:
        fwd_cmd += ' --overwrite '
    fwd_cmd += params['spacing']
    fwd_cmd += subj + params['bem']

    bash_script.append('' + params['forward_model'])
    bash_script.append('if [[ $? != 0 ]] ; then exit 1; fi')
