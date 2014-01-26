# -*- coding: utf-8 -*-
"""
Make the individual source spaces and prepare forward models.
Generate high-res head surface for coregistration.
Next step: complete forward model by including head-MRI transformation.

Output: 
scratch/fs_subjects_dir/SUB_ID/bem/watershed/*
scratch/fs_subjects_dir/SUB_ID/bem/SUB_ID_head.fif
scratch/fs_subjects_dir/SUB_ID/bem/SUB_ID_***-sol.fif

Created on Sun 26 Jan 2014

@author: cjayb
"""

from database import Query
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

bash_script = ['#!/usr/bin/env bash']
bash_script.append('export SUBJECTS_DIR='+subjects_dir)

params = {'source_space': '--ico -6', 
          'forward_model': '--homog --surf --ico 4',
          'morph_maps': True,
          'highres_head': True,
          'force': False}


for subj in ad.analysis_dict.keys():

    if not any('recon-all' in item for item in anadict.analysis_dict[subj].keys()):
        print "Skipping %s due to missing recon-all reconstruction" % subj
        continue

    bash_script.append('export SUBJECT=' + subj)
    #echo $SUBJECT
    bash_script.append('mne_watershed_bem --overwrite')
	
    cmd = '''
	cd ${SUBJECTS_DIR}/${SUBJECT}/bem
	ln -s watershed/${SUBJECT}_inner_skull_surface ${SUBJECT}-inner_skull.surf
	ln -s watershed/${SUBJECT}_outer_skin_surface ${SUBJECT}-outer_skin.surf
	ln -s watershed/${SUBJECT}_outer_skull_surface ${SUBJECT}-outer_skull.surf
	cd ''' + self._project_folder
    bash_script_append(cmd)
    
    cmd = 'mne_setup_source_space ' + params['source_space']
    if params['force']:
        cmd += ' --overwrite'
    bash_script_append(cmd)
	
	# Prepare for forward computation
    cmd = 'mne_setup_forward_model ' + params['forward_model']
    bash_script_append(cmd)
	
	# Generate morph maps for morphing between daniel and fsaverage
    cmd = 'mne_make_morph_maps --from ${SUBJECT} --to fsaverage'
    bash_script_append(cmd)

    cmd = '''
    cd ${SUBJECTS_DIR}/${SUBJECT}/bem
    head=${SUBJECTS_DIR}/${SUBJECT}/bem/${SUBJECT}-head.fif
    head_low=${SUBJECTS_DIR}/${SUBJECT}/bem/${SUBJECT}-head-lowres.fif
    if [ -e $head ]; then
      printf '\nmoving existing head surface %s\n' $head
      mv $head $head_low
    fi
    ${MNE_PYTHON}/bin/mne_make_scalp_surfaces.py -s ${SUBJECT} -o
    head_medium=${SUBJECTS_DIR}/${SUBJECT}/bem/${SUBJECT}-head-medium.fif
    printf '\nlinking %s as main head surface\n' $head_medium
    ln -s $head_medium $head
    '''
    bash_script_append(cmd)
    
    # This assumes the key does not already exist, otherwise it will be overwritten!
    ad.analysis_dict[subj].update({'watershed_bem': {}})
    
    # Start with a fresh copy of the defaults
    wb_params = params.copy()
    
    ad.analysis_dict[subj].update({'params': params})
        

#ad.apply_freesurfer('recon-all_initial', fake=False, verbose=True, n_processes=5)
