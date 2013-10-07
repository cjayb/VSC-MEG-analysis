# -*- coding: utf-8 -*-
"""
Fill in the dictionary for freesurfer brain extraction and labeling

@author: cjb
"""

from database import Query
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

VERBOSE=True

db = Query(proj_code=proj_code,verbose=True)
anadict = Anadict(db, verbose=False)    

recon_all_cmd = '/opt/local/freesurfer-releases/5.3.0/bin/recon-all'
subjects_dir = anadict._scratch_folder + 'fs_subjects_dir'
check_path_exists(subjects_dir)

fs_params_defaults = {'input_file': None, 
                    'logfile': None, 'use_gpu': True, 'num_threads': 4}


for subj in anadict.analysis_dict.keys():

    # This assumes the key does not already exist, otherwise it will be overwritten!
    anadict.analysis_dict[subj].update({'recon-all_initial': {}})
    
    cur_dict = anadict.analysis_dict[subj]['recon-all_initial']

    # Start with a fresh copy of the defaults
    # needs to be here in case mf_params was altered (empty_room)
    fs_params = fs_params_defaults.copy()
    fs_params['input_file'] = anadict.analysis_dict[subj]['T1']['files'][0] # first DICOM file
    cur_dict.update({'fs_params': fs_params})
        
#anadict.save('Added process dictionary for initial recon-all run.')    