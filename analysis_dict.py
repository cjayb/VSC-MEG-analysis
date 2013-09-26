# -*- coding: utf-8 -*-
"""
Create a dictionary with relevant files

Created on Thu Sep 26 12:15:11 2013

@author: cjb
"""

#class Anadict_item():
#    def __init__(self, verbose=False):
#        
#        self.

import pickle
import os
import subprocess as subp

class Anadict():
    
    def __init__(self, db, verbose=False):
            
        try: 
            with open(db.analysis_dict_name, "rb"):
                if verbose: 
                    print 'Reading from existing analysis dictionary:'
                    print db.analysis_dict_name
                    
                self.analysis_dict = pickle.load(open(db.analysis_dict_name, "rb"))                
                    
        except IOError:
            print 'No analysis dictionary found in'
            print db.analysis_dict_name
            print 'Creating a new dictionary'
        
            db.get_subjects(verbose=False)
            
            self.analysis_dict = {}
            for subj in db.subjects:
                meg_study = db.get_studies(subj, modality='MEG',unique=True)
                series = db.get_series(subj, meg_study, 'MEG')
                
                for ser in series:
                    meg_file = db.get_files(subj, meg_study, 'MEG',ser)
            
                    self.analysis_dict.update({subj: {ser: {"raw": meg_file}}})
                        

            pickle.dump(self.analysis_dict, open(db.analysis_dict_name, "wb"))
            curdir = os.getcwd()
            os.chdir(db.analysis_dict_folder)
            cmd = 'git add %s' % db.analysis_dict_name
            subp.call(cmd, shell=True)
            cmd = 'git commit -m \"Analysis dict initialised.\"'
            subp.call(cmd, shell=True)
            os.chdir(curdir)

        self.analysis_dict_folder = db.analysis_dict_folder
        self.analysis_dict_name = db.analysis_dict_name
