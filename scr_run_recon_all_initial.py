# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:18:16 2013

@author: cjb
"""
from mindlab_dicomdb_access.database import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

ad.attach_T1_images(db)
#foo!