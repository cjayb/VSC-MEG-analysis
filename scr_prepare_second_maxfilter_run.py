# -*- coding: utf-8 -*-
"""
Fill in the dictionary for initial maxfiltering:
- tsss_mc
- autobad on
- estimate origin (y<70mm)

Created on Thu Sep 26 10:25:28 2013

@author: cjb
"""

from database import Query
from analysis_dict import Anadict
from maxfilter_cfin import fit_sphere_to_headshape

from mne.fiff import Raw
from mne.utils import set_log_level as mne_set_log_level

#from sys import exit as sysexit
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

db = Query(proj_code=proj_code,verbose=True)
anadict = Anadict(db, verbose=False)    


for subj in anadict.analysis_dict.keys():

    try:
        tsss_initial_dict = anadict.analysis_dict[subj]['tsss_initial']
    except:
        raise Exception("analysis_missing")

    # This assumes the key does not already exist...
    anadict.analysis_dict[subj].update({'tsss_second': {}})
    
    for task in tsss_initial_dict.keys():

        anadict.analysis_dict[subj]['tsss_second'].update({task: {}})
        
        cur_dict = anadict.analysis_dict[subj]['tsss_second'][task]
        
        task_output_files = []
        task_mf_params = [] #NB: this is a list!!

        for ii_pars,init_pars in enumerate(sorted(tsss_initial_dict[task]['mf_params'])):
            fnum_raw = "%02d" % ii_pars

            # Start with a fresh copy of the defaults
            mf_params = init_pars.copy()
            
            output_folder = anadict._scratch_folder + '/tsss_second/' + subj
            check_path_exists(output_folder)
            if len() > 1:
                output_name_base = output_folder + '/'+task+ '-' + fnum_raw
            else:
                output_name_base = output_folder + '/'+task
            
            if not 'empty' in task:
                mf_params['output_file'] = output_name_base + '_tsss_mc.fif'
                mf_params['mv_hp'] = output_name_base + '_tsss_mc.pos'
                mf_params['logfile'] = output_name_base + '_tsss_mc.log'
            else:
                mf_params['output_file'] = output_name_base + '_tsss.fif'
                mf_params['mv_hp'] = None
                mf_params['logfile'] = output_name_base + '_tsss.log'
            
            # Since both task_input and task_output_files are lists, they
            # will now remain ordered 1-to-1
            task_output_files.append(mf_params['output_file'])
            task_mf_params.append(mf_params.copy())
            
        cur_dict.update({'files': task_output_files})
        cur_dict.update({'mf_params': task_mf_params})
        
#anadict.save('Added process dictionary for second tSSS run.')    