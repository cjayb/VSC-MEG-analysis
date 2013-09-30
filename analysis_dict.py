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
import errno
import subprocess as subp

class Anadict():
    
    def __init__(self, db, verbose=False):
            
        self.proj_code = db.proj_code
        self._project_folder = '/projects/' + self.proj_code
        self._raw_folder = '/projects/' + self.proj_code + '/raw'
        self._scratch_folder = '/projects/' + self.proj_code + '/scratch'
        self._misc_folder = '/projects/' + self.proj_code + '/misc'
        self._result_folder = '/projects/' + self.proj_code + '/result'
        self._doc_folder = '/projects/' + self.proj_code + '/doc'
        self._scripts_folder = '/projects/' + self.proj_code + '/scripts'
        #self.subjects = None

        self.analysis_dict_folder = self._misc_folder + '/analysis_dict'
        self.analysis_dict_name   = self.analysis_dict_folder + '/analysis_dict.pickle'        

        try: 
            with open(self.analysis_dict_name, "rb"):
                if verbose: 
                    print 'Reading from existing analysis dictionary:'
                    print self.analysis_dict_name
                    
                self.analysis_dict = pickle.load(open(self.analysis_dict_name, "rb"))                
                    
        except IOError:
            print 'No analysis dictionary found in'
            print self.analysis_dict_name
            print 'Creating a new dictionary'

            self._make_ad_folder() # Check that folder exists, initalize git
        
            subjects = db.get_subjects(verbose=False)
            
            self.analysis_dict = {}
            for subj in subjects:
                
                self.analysis_dict.update({subj: {}})                
                
                meg_study = db.get_studies(subj, modality='MEG',unique=True)
                series = db.get_series(subj, meg_study, 'MEG',verbose=verbose) # This is a 2D list with [series_number, series_name]
                
                ser_dict = {}
                for ser in series:
                    task_name = ser[0]
                    task_num = ser[1]
                    meg_file_names = db.get_files(subj, meg_study, 'MEG',task_num) 
                    
                    ser_dict.update({task_name : {'files' : meg_file_names}})
                    #ser_dict.update({ser[0]: {"raw": {"files": meg_file_names}}}) # meg_file_names is a list, possibly of just one
                    
                self.analysis_dict[subj].update({'raw': ser_dict})
            
            self.save('Analysis dict initialised.')
                        

#            pickle.dump(self.analysis_dict, open(db.analysis_dict_name, "wb"))
#            curdir = os.getcwd()
#            os.chdir(db.analysis_dict_folder)
#            cmd = 'git add %s' % db.analysis_dict_name
#            subp.call(cmd, shell=True)
#            cmd = 'git commit -m \"Analysis dict initialised.\"'
#            subp.call(cmd, shell=True)
#            os.chdir(curdir)


    def _commit_to_git(self, commit_message):
        curdir = os.getcwd()
        os.chdir(self.analysis_dict_folder)
        cmd = 'git add %s' % self.analysis_dict_name
        subp.call(cmd, shell=True)
        cmd = 'git commit -m \"' + commit_message + '\"'
        subp.call(cmd, shell=True)
        os.chdir(curdir)

    def _make_ad_folder(self):
        
        try: 
            os.makedirs(self.analysis_dict_folder)        
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(self.analysis_dict_folder):
                pass
            else:
                raise
        curdir = os.getcwd()
        os.chdir(self.analysis_dict_folder)
        subp.call('git init', shell=True)
        os.chdir(curdir)
   
    def save(self, commit_message, verbose=False):
        
        pickle.dump(self.analysis_dict, open(self.analysis_dict_name, "wb"))
        self._commit_to_git(commit_message)