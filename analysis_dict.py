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

from maxfilter_cfin import build_maxfilter_cmd, apply_maxfilter_cmd

import logging
import time

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

    def attach_T1_images(self, db, sequence_name='t1_mprage_3D_sag', 
                         verbose=False, save=False):
        """
        Sets up a "T1" key in the dictionary (overwrites if already exists)
        """

        if verbose:
            logging.basicConfig(level=logging.INFO)
        else:
            logging.basicConfig(level=logging.ERROR)
            

        subjects = db.get_subjects(verbose=False)
        for subj in subjects:
            if subj in self.analysis_dict:
                logging.info('Subject %s already present, augmenting', subj)
            else:
                logging.info('Adding subject %s to dictionary', subj)
                self.analysis_dict.update({subj: {}})
                    
            cur_ana_dict = self.analysis_dict[subj]
            
            mr_study = db.get_studies(subj, modality='MR',unique=True)
            if not mr_study is None:
                cur_ana_dict.update({'T1': {}})
                series = db.get_series(subj, mr_study, 'MR',verbose=verbose) # This is a 2D list with [series_name, series_number]
                for ser in series:
                    if sequence_name in ser:
                        T1_file_names = db.get_files(subj, mr_study, 'MR',ser[1]) 
                        cur_ana_dict['T1'].update({'files':T1_file_names})
        if save:            
            self.save('T1 images attached to each subject.')    

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
        
    def apply_maxfilter(self, analysis_name, force=None, verbose=False, fake=False):
        '''
        Apply a maxfilter-analysis that's already in the dictionary.
        
        Arguments:
            force: Re-set parameter "force" to (supercedes what's in the dictionary)
        '''
        log_name = self._scratch_folder + "/" + analysis_name + '/analysis.log'
        logging.basicConfig(filename=log_name, level=logging.INFO)

        timestr = time.asctime(time.localtime(time.time()))
        logging.info('Analysis log for %s, started on %s' % (analysis_name, timestr))        
        
        # Check that input files exist etc
        for subj in self.analysis_dict.keys():
            try:
                cur_ana_dict = self.analysis_dict[subj][analysis_name]
            except KeyError:
                logging.error('Subject %s is missing the analysis \"%s\"' % (subj, analysis_name))
                raise Exception("subject_missing")

            logging.info('Entering subject %s' % subj)

            for task in cur_ana_dict.keys():
                
                cur_mf_params = cur_ana_dict[task]['mf_params']

                # Loop over all files (usually just 1) in this task                
                for ii_mfp,mfp in enumerate(cur_mf_params):

                    if not os.path.isfile(mfp['input_file']):
                        logging.error("Following input file does not exist!")
                        logging.error(mfp['input_file'])
                        raise Exception("input_file_missing")

                    if (os.path.isfile(mfp['output_file']) and not mfp['force']):
                        if force is None:
                            logging.error("Output file exists, but option: force is False")
                            logging.error(mfp['output_file'])
                            logging.error("Set force-argument to override...")
                            raise Exception("output_file_exists")
                        elif force is True:
                            mfp['force'] = True

                    mf_cmd = build_maxfilter_cmd(mfp['input_file'], mfp['output_file'], origin=mfp['origin_head'], frame='head',
                                                 bad=mfp['bad'], autobad=mfp['autobad'], force=mfp['force'],
                                                 st=mfp['st'], st_buflen=mfp['st_buflen'], st_corr=mfp['st_corr'], 
                                                 movecomp=mfp['movecomp'], mv_hp=mfp['mv_hp'], hpicons=mfp['hpicons'],
                                                 linefreq=mfp['linefreq'], cal=mfp['cal'], ctc=mfp['ctc'],
                                                 verbose=mfp['verbose'], maxfilter_bin=mfp['maxfilter_bin'],
                                                 logfile=mfp['logfile'])

                    logging.info('Initiating Maxfilter with following command')
                    logging.info(mf_cmd)
                    time_then = time.time()
                    apply_maxfilter_cmd(mf_cmd)
                    time_now = time.time()
                    mfp.update({'command': mf_cmd}) # This will simply overwrite an existing "command"
                    logging.info('Task %s for subject %s completed in %.0f seconds' % (task, subj, time_now-time_then))
        
        self.save('Analysis name: %s completed and committed to git' % analysis_name)

    def apply_freesurfer(self, analysis_name, force=None, verbose=False, fake=False):
        '''
        Apply a freesurfer-analysis that's already in the dictionary.
        
        Arguments:
            force: Re-set parameter "force" to (supercedes what's in the dictionary)
        '''
        if verbose:
            log_level=logging.INFO
        else:
            log_level=logging.ERROR
                
        #log_name = self._scratch_folder + "/" + analysis_name + '/analysis.log'
        #logging.basicConfig(filename=log_name, level=logging.INFO)
        logging.basicConfig(level=log_level)

        # Check that input files exist etc
        for subj in self.analysis_dict.keys():
            try:
                cur_ana_dict = self.analysis_dict[subj][analysis_name]
            except KeyError:
                logging.error('Subject %s is missing the analysis \"%s\"' % (subj, analysis_name))
                raise Exception("subject_missing")

            logging.info('Entering subject %s' % subj)
            
            cur_fs_params = cur_ana_dict[task]['fs_params']

            fs_cmd = 'SUBJECTS_DIR=' + cur_fs_params['subjects_dir'] + ' ' \
                    + cur_fs_params['fs_bin'] + ' ' + cur_fs_params['fs_args']
                    
            if cur_fs_params['input_file']:
                fs_cmd += ' -i ' + cur_fs_params['input_file']
            if cur_fs_params['use_gpu']:
                fs_cmd += ' -use-gpu '
            if cur_fs_params['num_threads']:
                fs_cmd += ' -openmp %d ' % cur_fs_params['num_threads']

            logging.info('Running command:')
            logging.info(fs_cmd)
            
            if not fake:
                st = os.system(cmd)
                if st != 0:
                    raise RuntimeError('MaxFilter returned non-zero exit status %d' % st)
        
            logger.info('[done]')
