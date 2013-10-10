# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:18:16 2013

@author: cjb
"""
from mindlab_dicomdb_access.database import Query
from analysis_dict import Anadict

import os
import errno

def check_path_exists(chkpath):
    
    try: 
        os.makedirs(chkpath)        
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(chkpath):
            pass
        else:
            raise


proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'
db=Query(proj_code)
ad=Anadict(db)


recon_all_bin = '/opt/local/freesurfer-releases/5.3.0/bin/recon-all'
subjects_dir = ad._scratch_folder + '/fs_subjects_dir'
check_path_exists(subjects_dir)

fs_params_defaults = {'input_file': None, 'use_gpu': True, 'num_threads': 8,
                        'fs_bin': recon_all_bin, 'subjects_dir': subjects_dir,
                        'fs_args': '-all', 'force': False}

# Run this if T1 images not yet attached
#ad.attach_T1_images(db, verbose=False, save=True)

for subj in ad.analysis_dict.keys():

    if not 'T1' in ad.analysis_dict[subj]:
        print "Skipping %s due to missing T1" % subj
        continue

    # This assumes the key does not already exist, otherwise it will be overwritten!
    ad.analysis_dict[subj].update({'recon-all_initial': {}})
    
    cur_dict = ad.analysis_dict[subj]['recon-all_initial']

    # Start with a fresh copy of the defaults
    # needs to be here in case mf_params was altered (empty_room)
    fs_params = fs_params_defaults.copy()
    fs_params['input_file'] = ad.analysis_dict[subj]['T1']['files'][0] # first DICOM file
    
    cur_dict.update({'fs_params': fs_params})
        

ad.apply_freesurfer('recon-all_initial', fake=False, verbose=True, n_processes=5)