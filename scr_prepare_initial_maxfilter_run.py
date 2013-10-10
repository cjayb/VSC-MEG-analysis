# -*- coding: utf-8 -*-
"""
Fill in the dictionary for initial maxfiltering:
- tsss_mc
- autobad on
- estimate origin (y<70mm)

Created on Thu Sep 26 10:25:28 2013

@author: cjb
"""

from mindlab_dicomdb_access.database import Query
from analysis_dict import Anadict
from maxfilter_cfin import fit_sphere_to_headshape

from mne.fiff import Raw

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

VERBOSE=True
SAVE=False

db = Query(proj_code=proj_code,verbose=True)
anadict = Anadict(db, verbose=False)    

mx_cmd = anadict._misc_folder + '/bin/maxfilter-2.2.15'
cal_db = anadict._misc_folder + '/databases/sss_cal.dat'
ctc_db = anadict._misc_folder + '/databases/ct_sparse.fif'

mf_params_defaults = {'input_file': None, 'output_file': None,
             'autobad': 'on', 'st': True, 'movecomp': True,
             'st_corr': 0.96, 'st_buflen': 16, 'mv_hp': False,
             'origin_head': [0,0,40], 'radius_head': None,
             'bad': None, 'hpicons': True, 'linefreq': 50.,
             'cal': cal_db, 'ctc': ctc_db,
             'force': True, 'verbose': True, 'maxfilter_bin': mx_cmd,
             'logfile': None}


for subj in anadict.analysis_dict.keys():

    in_fname = anadict.analysis_dict[subj]['raw']['FFA']['files'][0] # Any file with HPI will do!
    raw = Raw(in_fname)
    radius_head, origin_head, origin_devive = fit_sphere_to_headshape(raw.info,ylim=0.070,verbose=VERBOSE)
    raw.close()        
    
    # This assumes the key does not already exist...
    anadict.analysis_dict[subj].update({'tsss_initial': {}})
    
    for task in anadict.analysis_dict[subj]['raw'].keys():

        # Start with a fresh copy of the defaults
        # needs to be here in case mf_params was altered (empty_room)
        mf_params = mf_params_defaults.copy()
        # Needs to be reset here if previous task was "empty_room"
        mf_params['movecomp'] = mf_params_defaults['movecomp']
        mf_params['hpicons'] = mf_params_defaults['hpicons']
        mf_params['origin_head'] = origin_head
        mf_params['radius_head'] = radius_head    

        anadict.analysis_dict[subj]['tsss_initial'].update({task: {}})
        cur_dict = anadict.analysis_dict[subj]['tsss_initial'][task]
        
        task_input_files = anadict.analysis_dict[subj]['raw'][task]['files']

        task_output_files = []
        task_mf_params = [] #NB: this is a list!!

        for ii_raw,raw_name in enumerate(sorted(task_input_files)):
            fnum_raw = "%02d" % ii_raw

            mf_params['input_file'] = raw_name
            
            output_folder = anadict._scratch_folder + '/tsss_initial/' + subj
            check_path_exists(output_folder)
            if len(task_input_files) > 1:
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
                mf_params['movecomp'] = False
                mf_params['hpicons'] = False
                mf_params['origin_head'] = False # Must be False, if None, will try to estimate it!
                mf_params['radius_head'] = False
            
            # Since both task_input and task_output_files are lists, they
            # will now remain ordered 1-to-1
            task_output_files.append(mf_params['output_file'])
            task_mf_params.append(mf_params.copy())
            
        cur_dict.update({'files': task_output_files})
        cur_dict.update({'mf_params': task_mf_params})

if SAVE:        
    anadict.save('Added process dictionary for initial tSSS run.')    