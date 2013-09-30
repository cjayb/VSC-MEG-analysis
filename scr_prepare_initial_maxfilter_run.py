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
from maxfilter_cfin import fit_sphere_to_headshape, apply_maxfilter

from mne.fiff import Raw
from mne.utils import set_log_level as mne_set_log_level

from sys import exit as sysexit
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
mne_set_log_level(verbose=False)

db = Query(proj_code=proj_code,verbose=True)
anadict = Anadict(db, verbose=False)    

mx_cmd = anadict._misc_folder + '/bin/maxfilter-2.2.15'
cal_db = anadict._misc_folder + '/databases/sss_cal.dat'
ctc_db = anadict._misc_folder + '/databases/ct_sparse.fif'

mf_params_defaults = {'input_file': None, 'output_file': None,
             'autobad': 'on', 'tsss': True, 'movecomp': True,
             'st_corr': 0.96, 'st_buflen': 16, 'mv_hp': False,
             'origin_head': [0,0,40], 'radius_head': None,
             'bad': [], 'hpicons': True, 'linefreq': 50.,
             'cal': cal_db, 'ctc': ctc_db,
             'overwrite': True, 'verbose': True, 'maxfilter_bin': mx_cmd,
             'logfile': None}


for subj in anadict.analysis_dict.keys():

    mf_params = mf_params_defaults.copy()

    in_fname = anadict.analysis_dict[subj]['raw']['FFA']['files'][0] # Any file with HPI will do!
    raw = Raw(in_fname)
    r, o_head, o_dev = fit_sphere_to_headshape(raw.info,ylim=0.070)
    raw.close()        
    
    mf_params['origin_head'] = o_head
    mf_params['radius_head'] = r    
    
    anadict.analysis_dict[subj].update({'tsss_initial': {}})
    
    for task in anadict.analysis_dict[subj]['raw'].keys():

        anadict.analysis_dict[subj]['tsss_initial'].update({task: {}})
        cur_dict = anadict.analysis_dict[subj]['tsss_initial'][task]
        
        task_input_files = anadict.analysis_dict[subj]['raw'][task]['files']

        task_output_files = []
        task_mf_params = []

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
                mf_params['movecomp'] = False
                mf_params['origin_head'] = False
                mf_params['logfile'] = output_name_base + '_tsss.log'
            
            task_mf_params.append(mf_params.copy()) # NB: .copy is important here!           
            
            # Since both task_input and task_output_files are lists, they
            # will now remain ordered 1-to-1
            task_output_files.append(mf_params['output_file'])
                    
        cur_dict.update({'files': task_output_files})
        cur_dict.update({'mf_params': task_mf_params})
        
anadict.save('Added process dictionary for initial tSSS run.')    