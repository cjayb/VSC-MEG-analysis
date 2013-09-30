# -*- coding: utf-8 -*-
"""
test maxfilter_cfin.py

Created on Fri Sep 27 22:30:51 2013

@author: cjb
"""

from database import Query
from analysis_dict import Anadict
from maxfilter_cfin import apply_maxfilter

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

db = Query(proj_code=proj_code)
anadict = Anadict(db)    

subj = '023_Q4V'
task = 'FFA'

out_fname = db._scratch_folder + '/tmp/' + subj + '.' + task + '_tsss_mc.fif'
out_hp = db._scratch_folder + '/tmp/' + subj + '.' + task + '_tsss_mc.pos'

mx_cmd = '/projects/' + proj_code + '/misc/bin/maxfilter-2.2.15'
cal_db = '/projects/' + proj_code + '/misc/databases/sss_cal.dat'
ctc_db = '/projects/' + proj_code + '/misc/databases/ct_sparse.fif'
logfile= db._scratch_folder + '/tmp/' + subj + '.' + task + '_tsss_mc.log'

#set_log_file(logfile)

apply_maxfilter(anadict.analysis_dict[subj][task]['raw']['files'][0], out_fname, 
                    origin=None, frame='head',
                    bad=None, autobad='on', skip=None, force=False,
                    st=True, st_buflen=16.0, st_corr=0.96, mv_trans=None,
                    mv_comp=True, mv_headpos=False, mv_hp=out_hp,
                    mv_hpistep=None, mv_hpisubt=None, mv_hpicons=True,
                    linefreq=None, cal=cal_db, ctc=ctc_db, mx_args='',
                    overwrite=True, verbose=True, maxfilter_cmd=mx_cmd,
                    logfile=logfile)
                    
                    