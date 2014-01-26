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

db = Query(proj_code=proj_code,verbose=True)
anadict = Anadict(db, verbose=False)    

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

recon_all_bin = '/opt/local/freesurfer-releases/5.3.0/bin/recon-all'
subjects_dir = ad._scratch_folder + '/fs_subjects_dir'
check_path_exists(subjects_dir)

VERBOSE=True

bash_script = ['#!/usr/bin/env bash']
bash_script.append('export SUBJECTS_DIR='+subjects_dir)

params = {'source_space': '--ico -6', 
          'forward_model': '--homog --surf --ico 4',
          'morph_maps': True,
          'highres_head': True,
          'force': False}

# Run this if T1 images not yet attached
#ad.attach_T1_images(db, verbose=False, save=True)

for subj in ad.analysis_dict.keys():

    if not 'T1' in ad.analysis_dict[subj]:
        print "Skipping %s due to missing T1" % subj
        continue

    bash_script.append('export SUBJECT=' + subj
	#echo $SUBJECT
	bash_script.append('mne_watershed_bem --overwrite')
	
    cmd = '''
	cd ${SUBJECTS_DIR}/${SUBJECT}/bem
	ln -s watershed/${SUBJECT}_inner_skull_surface ${SUBJECT}-inner_skull.surf
	ln -s watershed/${SUBJECT}_outer_skin_surface ${SUBJECT}-outer_skin.surf
	ln -s watershed/${SUBJECT}_outer_skull_surface ${SUBJECT}-outer_skull.surf
	cd ''' + self._project_folder)
	bash_script_append(cmd)
    
    cmd = 'mne_setup_source_space ' + params['source_space']
    cmd += ' --overwrite' if params['force']
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
    ad.analysis_dict[subj].update({'recon-all_initial': {}})
    
    cur_dict = ad.analysis_dict[subj]['recon-all_initial']

    # Start with a fresh copy of the defaults
    # needs to be here in case mf_params was altered (empty_room)
    fs_params = fs_params_defaults.copy()
    fs_params['input_file'] = ad.analysis_dict[subj]['T1']['files'][0] # first DICOM file
    
    cur_dict.update({'fs_params': fs_params})
        

ad.apply_freesurfer('recon-all_initial', fake=False, verbose=True, n_processes=5)