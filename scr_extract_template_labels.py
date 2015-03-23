"""
Strangely, mri_annotation2labels doesn't allow selecting which label to extract!
Created on Wed 11 Mar 2015

@author: cjayb
"""

from database import Query
from analysis_dict import Anadict

# from mne import find_events, write_events
# from mne.fiff import Raw, pick_channels

import numpy as np
import multiprocessing
import subprocess

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

db = Query(proj_code=proj_code,verbose=True)
ad = Anadict(db, verbose=False)    

fs_bin = '/opt/local/freesurfer-releases/5.3.0/bin'
subjects_dir = ad._scratch_folder + '/fs_subjects_dir'

# Arno Klein, Jason Tourville. Frontiers in Brain Imaging Methods. 
# 6:171. DOI: 10.3389/fnins.2012.00171 
template = 'aparc.DKTatlas40'

n_processes = 4


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

VERBOSE=True
fake = False

if not fake:
    pool = multiprocessing.Pool(processes=n_processes)

all_cmds=[]
for subj in db.get_subjects():
    if len(subj) == 8:
        subj = subj[1:]
    bash_script = ['#!/usr/bin/env bash']
    bash_script.append('source ~/.bashrc')
    bash_script.append('use mne')
    bash_script.append('export SUBJECTS_DIR=' +subjects_dir)
    bash_script.append('export TEMPLATE=' + template)

    if not any('recon-all' in item for item in ad.analysis_dict[subj].keys()):
        print "Skipping %s due to missing recon-all reconstruction" % subj
        continue

    bash_script.append('export SUBJECT=' + subj)
    bash_script.append('mkdir -p ${SUBJECTS_DIR}/${SUBJECT}/label/DKT40_labels')

    cmd = '''
mri_annotation2label --annotation ${TEMPLATE} --subject ${SUBJECT} --hemi lh --outdir ${SUBJECTS_DIR}/${SUBJECT}/label/DKT40_labels
mri_annotation2label --annotation ${TEMPLATE} --subject ${SUBJECT} --hemi rh --outdir ${SUBJECTS_DIR}/${SUBJECT}/label/DKT40_labels
cd ${SUBJECTS_DIR}/${SUBJECT}/label
rm *fusiform.label
ln -s DKT40_labels/*fusiform.label .
'''
    bash_script.append(cmd)

    all_cmds.append('\n'.join(bash_script))

if not fake:
    return_codes = pool.map(_parallel_task,all_cmds)
    pool.close()
    pool.join()
elif VERBOSE:
    print "The following would execute, if this were not a FAKE run:"
    for cmd in all_cmds:
        print "%s" % cmd

