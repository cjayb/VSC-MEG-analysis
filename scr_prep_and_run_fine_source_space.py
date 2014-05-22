# -*- coding: utf-8 -*-
"""
Create an extra high-res source space (ico 5)

NB: this fails for 009 (something about a vertex in the map multiple times)
 I've run mne_setup_source_space --spacing 4
 on her, which seems to produce a similar number of dipoles per hemisphere (ca. 10k) as ico 5

NBNB!: I've symlinked the even spacing source space to appear as an ico 5 for 009!!
    This is a bit dodgy, but as I understand the homeomorphism to a sphere is not used (yet) in MNE...?

@author: cjayb
"""

#from database import Query
import sys
sys.path.append('/users/cjb/src/PyCharmProjects')
from stormdb.access import Query

from analysis_dict import Anadict

# from mne import find_events, write_events
# from mne.fiff import Raw, pick_channels

import numpy as np

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

db = Query(proj_code=proj_code,verbose=True)
ad = Anadict(db, verbose=False)

recon_all_bin = '/opt/local/freesurfer-releases/5.3.0/bin/recon-all'
subjects_dir = ad._scratch_folder + '/fs_subjects_dir'
# check_path_exists(subjects_dir)

VERBOSE=True

params = {'source_space': '--ico 5',
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
    #echo $SUBJECT

    cmd = 'mne_setup_source_space ' + params['source_space']
    if params['force']:
        cmd += ' --overwrite'
    bash_script.append(cmd)
    bash_script.append('if [[ $? != 0 ]] ; then exit 2; fi')

    # This assumes the key does not already exist, otherwise it will be overwritten!
    ad.analysis_dict[subj].update({'highres_source_space': {}})

    # Start with a fresh copy of the defaults
    cur_params = params.copy()

    ad.analysis_dict[subj]['highres_source_space'].update({'params': cur_params})
    ad.analysis_dict[subj]['highres_source_space'].update({'command': bash_script})

ad.apply_bash_script('highres_source_space', n_processes=6)