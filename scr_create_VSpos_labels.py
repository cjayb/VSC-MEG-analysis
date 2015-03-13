"""
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
# check_path_exists(subjects_dir)
template_dir = ad._misc_folder + '/templates'
polar_templ  = template_dir + '/mh.V1.poltmp.sym.mgh'
eccen_templ  = template_dir + '/mh.V1.ecctmp.sym.mgh'
VSLUT  = template_dir + '/VSLUT.txt'
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
for subj in ad.analysis_dict.keys():
    bash_script = ['#!/usr/bin/env bash']
    bash_script.append('source ~/.bashrc')
    bash_script.append('use mne')
    bash_script.append('export SUBJECTS_DIR=' +subjects_dir)

    if not any('recon-all' in item for item in ad.analysis_dict[subj].keys()):
        print "Skipping %s due to missing recon-all reconstruction" % subj
        continue

    bash_script.append('export SUBJECT=' + subj)
    bash_script.append('export POLTMP=' + polar_templ)
    bash_script.append('export ECCTMP=' + eccen_templ)
    bash_script.append('export VSLUT=' + VSLUT)
    bash_script.append('export MGH_OUT=${SUBJECTS_DIR}/${SUBJECT}/label')

    # This takes the longest
    bash_script.append('surfreg --s '+subj+' --t fsaverage_sym --lh')
    bash_script.append('if [[ $? != 0 ]] ; then exit 1; fi')
    bash_script.append('surfreg --s '+subj+' --t fsaverage_sym --lh --xhemi')
    bash_script.append('if [[ $? != 0 ]] ; then exit 2; fi')
    
    cmd = '''
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $SUBJECT --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval ${POLTMP} --tval ${MGH_OUT}/lh.V1.poltmp.mgh
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $SUBJECT --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval ${ECCTMP} --tval ${MGH_OUT}/lh.V1.ecctmp.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${SUBJECT}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval ${POLTMP} --tval ${MGH_OUT}/rh.V1.poltmp.mgh
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${SUBJECT}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval ${ECCTMP} --tval ${MGH_OUT}/rh.V1.ecctmp.mgh
''' + ad._project_folder
    bash_script.append(cmd)


# remember: lh has positions 6,1 & 2, where 0deg=up and pos 1=90 deg!
# remember: rh has positions 5,4 & 3, where 0deg=up and pos 4=90 deg!
    cmd = '''
cd ${MGH_OUT} 
mkdir -p VSpos
mri_binarize --i lh.V1.poltmp.mgh --min 0.0 --max 60.0 --o VSpos/lh.V1.VSpos6.nii
mri_binarize --i lh.V1.poltmp.mgh --min 60.0 --max 120.0 --o VSpos/lh.V1.VSpos1.nii
mri_binarize --i lh.V1.poltmp.mgh --min 120.0 --max 180.0 --o VSpos/lh.V1.VSpos2.nii
mri_binarize --i rh.V1.poltmp.mgh --min 0.0 --max 60.0 --o VSpos/rh.V1.VSpos5.nii
mri_binarize --i rh.V1.poltmp.mgh --min 60.0 --max 120.0 --o VSpos/rh.V1.VSpos4.nii
mri_binarize --i rh.V1.poltmp.mgh --min 120.0 --max 180.0 --o VSpos/rh.V1.VSpos3.nii

mri_binarize --i lh.V1.ecctmp.mgh --min 3. --max 20.0 --o VSpos/lh.V1.VSecc.nii
mri_binarize --i rh.V1.ecctmp.mgh --min 3. --max 20.0 --o VSpos/rh.V1.VSecc.nii

fscalc VSpos/lh.V1.VSpos1.nii mul VSpos/lh.V1.VSecc.nii --o VSpos/lh.V1.VS1.nii
fscalc VSpos/lh.V1.VSpos2.nii mul VSpos/lh.V1.VSecc.nii mul 2 --o VSpos/lh.V1.VS2.nii
fscalc VSpos/lh.V1.VSpos6.nii mul VSpos/lh.V1.VSecc.nii mul 6 --o VSpos/lh.V1.VS6.nii
fscalc VSpos/rh.V1.VSpos3.nii mul VSpos/rh.V1.VSecc.nii mul 3 --o VSpos/rh.V1.VS3.nii
fscalc VSpos/rh.V1.VSpos4.nii mul VSpos/rh.V1.VSecc.nii mul 4 --o VSpos/rh.V1.VS4.nii
fscalc VSpos/rh.V1.VSpos5.nii mul VSpos/rh.V1.VSecc.nii mul 5 --o VSpos/rh.V1.VS5.nii

fscalc VSpos/lh.V1.VS1.nii add VSpos/lh.V1.VS2.nii add VSpos/lh.V1.VS6.nii --o VSpos/lh.V1.VSpos.nii
fscalc VSpos/rh.V1.VS3.nii add VSpos/rh.V1.VS4.nii add VSpos/rh.V1.VS5.nii --o VSpos/rh.V1.VSpos.nii

# Note that --o assumes we want it to go into label
mris_seg2annot --seg VSpos/lh.V1.VSpos.nii --ctab ${VSLUT} --s ${SUBJECT} --h lh --o lh.V1.VSpos.annot
mris_seg2annot --seg VSpos/rh.V1.VSpos.nii --ctab ${VSLUT} --s ${SUBJECT} --h rh --o rh.V1.VSpos.annot

mri_annotation2label --annotation V1.VSpos --subject ${SUBJECT} --hemi lh --outdir ${SUBJECTS_DIR}/${SUBJECT}/label
mri_annotation2label --annotation V1.VSpos --subject ${SUBJECT} --hemi rh --outdir ${SUBJECTS_DIR}/${SUBJECT}/label
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

