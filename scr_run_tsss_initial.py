# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:18:16 2013

@author: cjb
"""
from database import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

ad.apply_maxfilter('tsss_initial')
#foo!