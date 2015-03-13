import sys
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

import mne
import os

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise
    
filter_params = {'input_files': 'tsss_initial', 'conditions': ['VS','FFA'],
                 'lowpass': 35.0, 'highpass': 0.5, 'method': 'iir', 'n_jobs': 6,
                 'l_trans_bandwidth': 0.4}


filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])
outdir_base = ad._scratch_folder + '/filtered/' + filter_params['input_files'] + '/' + filt_dir
#mkdir_p(ad._scratch_folder + '/filtered/' + '/' + filter_params['input_files']) + '/' + filt_dir 


for subj in ad.analysis_dict.keys():

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
                    
