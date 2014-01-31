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
import tempfile

from maxfilter_cfin import build_maxfilter_cmd, apply_maxfilter_cmd

import logging
import sys
import time

import subprocess
import multiprocessing

logger = logging.getLogger('anadict')
logger.propagate=False
try:
    if not logger.handlers:
        stdout_stream = logging.StreamHandler(sys.stdout)
        logger.addHandler(stdout_stream)
except:
    pass

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
                         verbose=False, save=False, force=False, subj_list=None):
        """
        Sets up a "T1" key in the dictionary (overwrites if already exists)
        """

        if verbose:
            log_level=logging.INFO
        else:
            log_level=logging.ERROR            
        logger.setLevel(log_level)

        if subj_list is None:            
            logger.info('Attaching T1 images to each subject.')
            subjects = db.get_subjects(verbose=False)
        else:
            if not type(subj_list) is list: 
                logger.error('subj_list must be a list!')
                raise Exception("input_error")
                
            logger.info('Attaching T1 images to subjects: %s' % ', '.join(subj_list))
            subjects = subj_list
        
        for subj in subjects:
            if subj in self.analysis_dict:
                logger.info('Subject %s already present, augmenting', subj)
            else:
                logger.info('Adding subject %s to dictionary', subj)
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
            if subj_list is None:            
                self.save('T1 images attached to each subject.')    
            else:
                self.save('T1 images attached to subjects: %s' % ', '.join(subj_list))    
                
            

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
        

    def apply_maxfilter(self, analysis_name, force=None, verbose=False, fake=False,
                        n_processes=1, subj_list=None):
        '''
        Apply a maxfilter-analysis that's already in the dictionary.
        
        Parameters:
            analysis_name: A dictionary key defining a maxfilter-analysis
        
        Arguments:
            n_processes: Integer defining how many to run in parallel
                         (on individual cores of the present computer, NOT clusterized yet)
            
            subj_list:   default is None, meaning apply analysis to all subjects in dictionary
                         if set, should be a list of subject names in the dictionary
                         (non-existent name leads to Exception)
                         
            force: Re-set parameter "force" to (supercedes what's in the dictionary)
        '''
        if verbose:
            log_level=logging.INFO
        else:
            log_level=logging.ERROR

        log_name = self._scratch_folder + "/" + analysis_name + '/analysis.log'
        logfile_stream = logging.FileHandler(log_name)
        logger.addHandler(logfile_stream)

        logger.setLevel(log_level)

        timestr = time.asctime(time.localtime(time.time()))
        logger.info('Analysis log for %s, started on %s' % (analysis_name, timestr))        
        
        if subj_list is None:
            subj_list = self.analysis_dict.keys()
        else:
            if not type(subj_list) is list: 
                logger.error('subj_list must be a list!')
                raise Exception("input_error")
            for subj in subj_list:
                if not subj in self.analysis_dict.keys():
                    logger.error('Subject %s unknown' % subj)
                    raise Exception("unknown_subject")

        if not fake:
            pool = multiprocessing.Pool(processes=n_processes)

        all_cmds=[]

        # Check that input files exist etc
        for subj in subj_list:
            try:
                cur_ana_dict = self.analysis_dict[subj][analysis_name]
            except KeyError:
                logger.error('Subject %s is missing the analysis \"%s\"' % (subj, analysis_name))
                raise Exception("analysis_missing")

            logger.info('Entering subject %s' % subj)


            for task in cur_ana_dict.keys():
                
                cur_mf_params = cur_ana_dict[task]['mf_params']

                # Loop over all files (usually just 1) in this task                
                for ii_mfp,mfp in enumerate(cur_mf_params):

                    if not os.path.isfile(mfp['input_file']):
                        logger.error("Following input file does not exist!")
                        logger.error(mfp['input_file'])
                        raise Exception("input_file_missing")

                    if (os.path.isfile(mfp['output_file']) and not mfp['force']):
                        if force is None:
                            logger.error("Output file exists, but option: force is False")
                            logger.error(mfp['output_file'])
                            logger.error("Set force-argument to override...")
                            raise Exception("output_file_exists")
                        elif force is True:
                            mfp['force'] = True

                    # NB! This would not work if any one of the parameters in mfp were not defined !
                    # Meaning: this is a very hard-coded way to define a "generic" maxfilter run...
                    mf_cmd = build_maxfilter_cmd(mfp['input_file'], mfp['output_file'], origin=mfp['origin_head'], frame='head',
                                                 bad=mfp['bad'], autobad=mfp['autobad'], force=mfp['force'],
                                                 st=mfp['st'], st_buflen=mfp['st_buflen'], st_corr=mfp['st_corr'], 
                                                 movecomp=mfp['movecomp'], mv_hp=mfp['mv_hp'], hpicons=mfp['hpicons'],
                                                 linefreq=mfp['linefreq'], cal=mfp['cal'], ctc=mfp['ctc'],
                                                 verbose=mfp['verbose'], maxfilter_bin=mfp['maxfilter_bin'],
                                                 logfile=mfp['logfile'])

                    logger.info('Initiating Maxfilter with following command')
                    logger.info(mf_cmd)
                    
                    #time_then = time.time() #This won't make sense in parallel
                    
                    #apply_maxfilter_cmd(mf_cmd)
                    all_cmds.append(mf_cmd)
                    
                    #time_now = time.time() #This won't make sense in parallel

                    mfp.update({'command': mf_cmd}) # This will simply overwrite an existing "command"                    
                    
                    #logging.info('Task %s for subject %s completed in %.0f seconds' % (task, subj, time_now-time_then))

        if subj_list is None:
            save_msg='Maxfilter run %s completed for all subjects.' % analysis_name
        else:
            save_msg='Maxfilter run %s completed for subjects: %s.' % \
                            (analysis_name, ', '.join(subj_list))
        
        if not fake:
            return_codes = pool.map(_parallel_task,all_cmds)
            pool.close()
            pool.join()

            self.save(save_msg)
            return return_codes
            
        elif verbose:
            print "The following would execute, if this were not a FAKE run:"
            for cmd in all_cmds:
                print "%s" % cmd
                
            print """self.save('%s')""" % save_msg
        
        logger.removeHandler(logfile_stream)
        

    def apply_freesurfer(self, analysis_name, force=None, verbose=False, fake=False,
                         n_processes=1, subj_list=None):
        '''
        Apply a freesurfer-analysis that's already in the dictionary.
        
        Parameters:
            analysis_name: A dictionary key defining a freesurfer-analysis
        
        Arguments:
            n_processes: Integer defining how many to run in parallel
                         (on individual cores of the present computer, NOT clusterized yet)
            
            subj_list:   default is None, meaning apply analysis to all subjects in dictionary
                         if set, should be a list of subject names in the dictionary
                         (non-existent name leads to Exception)
                         
            force: Re-set parameter "force" to (supercedes what's in the dictionary)
        '''
        if verbose:
            log_level=logging.INFO
        else:
            log_level=logging.ERROR
                
        logger.setLevel(log_level)

        if subj_list is None:
            subj_list = self.analysis_dict.keys()
        else:
            if not type(subj_list) is list: 
                logger.error('subj_list must be a list!')
                raise Exception("input_error")

            for subj in subj_list:
                if not subj in self.analysis_dict.keys():
                    logger.error('Subject %s unknown' % subj)
                    raise Exception("unknown_subject")

        if not fake:
            pool = multiprocessing.Pool(processes=n_processes)

        all_cmds=[]

        # Check that input files exist etc
        for subj in subj_list:
            try:
                cur_ana_dict = self.analysis_dict[subj][analysis_name]
            except KeyError:
                logger.info('Subject %s is missing the analysis \"%s\"' % (subj, analysis_name))
                logger.info('Skipping subject...')
                continue                                
                #raise Exception("subject_missing")

            logger.info('Entering subject %s' % subj)
            cur_fs_params = cur_ana_dict['fs_params']

            fs_cmd = 'SUBJECTS_DIR=' + cur_fs_params['subjects_dir'] + ' ' \
                    + cur_fs_params['fs_bin'] + ' ' + cur_fs_params['fs_args'] \
                    + ' -s ' + subj 
                    
            if cur_fs_params['input_file']:
                fs_cmd += ' -i ' + cur_fs_params['input_file']
            if cur_fs_params['use_gpu']:
                fs_cmd += ' -use-gpu '
            if cur_fs_params['num_threads']:
                fs_cmd += ' -openmp %d ' % cur_fs_params['num_threads']

            if os.path.exists(cur_fs_params['subjects_dir'] + '/' + subj) and not force:
                logger.info('Subject %s appears to be done, skipping (use force to overwrite)' % subj)
                
                # This is a bugfix, remove from final master
                # if not 'fs_cmd' in cur_fs_params:
                #     cur_fs_params.update({'fs_cmd': fs_cmd})
                #     logger.info('Adding fs_cmd to previously performed run (DEBUG ONLY)')
                #     #############                
                    
                continue

            cur_fs_params.update({'fs_cmd': fs_cmd})

            # DEBUG
            #fs_cmd = './reveal_pid.sh'

            logger.info('Running command:')
            logger.info(fs_cmd)
            
            all_cmds.append(fs_cmd)
            #st = os.system(fs_cmd)
            #if st != 0:
            #    raise RuntimeError('Freesurfer returned non-zero exit status %d' % st)
                
            logger.info('[done]')

        if subj_list is None:
            save_msg='Freesurfer run %s completed for all subjects.' % analysis_name
        else:
            save_msg='Freesurfer run %s completed for subjects: %s.' % \
                            (analysis_name, ', '.join(subj_list))
            
        
        if not fake:
            return_codes = pool.map(_parallel_task,all_cmds)
            pool.close()
            pool.join()

            self.save(save_msg)
            return return_codes
            
        elif verbose:
            print "The following would execute, if this were not a FAKE run:"
            for cmd in all_cmds:
                print "%s" % cmd
                
            print """self.save('%s')""" % save_msg
            
    def apply_bash_script(self, analysis_name, force=None, verbose=False, fake=False,
                          n_processes=1, subj_list=None):
        '''
        Apply a bash script that's already in the dictionary.
        
        Parameters:
            analysis_name: A dictionary key defining a bash based-analysis
        
        Arguments:
            n_processes: Integer defining how many to run in parallel
                         (on individual cores of the present computer, NOT clusterized yet)
            
            subj_list:   default is None, meaning apply analysis to all subjects in dictionary
                         if set, should be a list of subject names in the dictionary
                         (non-existent name leads to Exception)
                         
            force: Re-set parameter "force" to (supercedes what's in the dictionary)
        '''
        if verbose:
            log_level=logging.INFO
        else:
            log_level=logging.ERROR
                
        logger.setLevel(log_level)

        if subj_list is None:
            subj_list = self.analysis_dict.keys()
        else:
            if not type(subj_list) is list: 
                logger.error('subj_list must be a list!')
                raise Exception("input_error")

            for subj in subj_list:
                if not subj in self.analysis_dict.keys():
                    logger.error('Subject %s unknown' % subj)
                    raise Exception("unknown_subject")

        if not fake:
            pool = multiprocessing.Pool(processes=n_processes)

        all_cmds = []
        all_script_files = []

        # Check that input files exist etc
        for subj in subj_list:
            try:
                cur_ana_dict = self.analysis_dict[subj][analysis_name]
            except KeyError:
                logger.info('Subject %s is missing the analysis \"%s\"' % (subj, analysis_name))
                logger.info('Skipping subject...')
                continue                                
                #raise Exception("subject_missing")

            logger.info('Entering subject %s' % subj)

### EDIT FROM HERE ON!
            script_file = tempfile.NamedTemporaryFile(delete=False)
            script_log_name = script_file.name + '.log'
            script_file.write('\n'.join(cur_ana_dict['command']))
            script_file.close() # file is not immediately deleted because we
                                # used delete = False
            all_script_files.append(script_file.name)
            all_script_files.append(script_log_name)
            all_cmds.append('sh ' + script_file.name + '  2>&1 > ' + script_log_name)

            # if os.path.exists(cur_fs_params['subjects_dir'] + '/' + subj) and not force:
            #     logger.info('Subject %s appears to be done, skipping (use force to overwrite)' % subj)
                

        if subj_list is None:
            save_msg='Script %s completed for all subjects.' % analysis_name
        else:
            save_msg='Script %s completed for subjects: %s.' % \
                            (analysis_name, ', '.join(subj_list))
            
        
        if not fake:
            return_codes = pool.map(_parallel_task,all_cmds)
            pool.close()
            pool.join()

            if any(return_codes):
                print "Some subprocesses didn't complete!"
                print "Return codes: ", return_codes
                print "Dying here, temp files not deleted, see below:"
                print '\n'.join(all_cmds)
                raise
            else:
                self.save(save_msg)
            
        elif verbose:
            print "The following would execute, if this were not a FAKE run:"
            for cmd in all_cmds:
                print "%s" % cmd
                
            print """self.save('%s')""" % save_msg
        
        # cleanup
        for fname in all_script_files:
            os.unlink(fname)

def _parallel_task(command):
    """
        General purpose method to submit Unix executable-based analyses (e.g.
        maxfilter and freesurfer) to the shell.
        
        Parameters:
            command:    The command to execute (single string)
                                        
        Returns:        The return code (shell) of the command
    """
    #proc = subprocess.Popen([fs_cmd],stdout=subprocess.PIPE, shell=True)
    proc = subprocess.Popen([command], shell=True)
    
    proc.communicate()
    return proc.returncode
