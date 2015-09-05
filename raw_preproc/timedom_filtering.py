"""
Apply to raw-like data
"""
import mne
#try:
from mne.io import Raw
from mne import pick_types
from mne.report import Report
import numpy as np

# Dirty hack to get the parent dir in the path!
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.curdir, os.path.pardir)))
###

# get basic stuff like mkdir_p and some defaults
from VSC_utils import *
# Redefine these though
filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 40.0, 'highpass': 1.0}

filt_dir = '{hi:.1f}-{lo:.1f}Hz'.format(hi=filter_params['highpass'],
                                        lo=filter_params['lowpass'])

for subj in db.get_subjects():
#############################################


    outdir = outdir_base + '/' + subj
    mkdir_p(outdir)

    cond_names = ad.analysis_dict[subj][filter_params['input_files']].keys()
    for cond in cond_names:
        for run_cond in filter_params['conditions']:
            if run_cond in cond:

                in_fnames = ad.analysis_dict[subj][filter_params['input_files']][cond]['files']
                for fname in in_fnames:
                    print 'In: ', fname
                    raw = mne.fiff.Raw(fname, preload=True)

                    raw.filter(filter_params['highpass'], filter_params['lowpass'],
                                method=filter_params['method'], n_jobs=filter_params['n_jobs'],
                                l_trans_bandwidth=filter_params['l_trans_bandwidth'])

                    out_fname = outdir + '/' + cond + '_filt' + '.fif'
                    print 'Out:', out_fname
                    raw.save(out_fname, format='single')
