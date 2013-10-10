# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:18:16 2013

@author: cjb
"""
from mindlab_dicomdb_access.database import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

return_codes = ad.apply_maxfilter('tsss_second', force=None, verbose=False, fake=True,
                   n_processes=4, subj_list=None)
#foo!