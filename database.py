"""
=========================
Methods to interact with the MINDLab subject database
=========================

"""
# Author: Chris Bailey <cjb@cfin.au.dk>
#
# License: BSD (3-clause)


import subprocess as subp
from sys import exit as sysexit

def _wget_error_handling(stdout):
    # Do error checking of wget output: will start with "error:"
    if 'error' in stdout:
        print "Something is wrong, database answers as follows (dying...):"
        print stdout
        return -1
        
    return 0

def _wget_system_call(url):
    cmd = 'wget -qO - test ' + url
    pipe = subp.Popen(cmd, stdout=subp.PIPE,stderr=subp.PIPE,shell=True)
    output, stderr = pipe.communicate()
    #output = subp.call([cmd,opts], shell=True)

    if _wget_error_handling(output) < 0:
        sysexit(-1)        

    return output

def get_subjects(mindlablogin, proj_code, subj_type='included'):
    
    if subj_type == 'all': #Doesn't work yet!
        scode = 'subjectswithcode'
    elif subj_type == 'included':
        scode = 'subjectswithcode'
    elif subj_type == 'excluded':
        scode = 'excludedsubjectswithcode'
    
    url = 'http://hyades00.pet.auh.dk/modules/Ptdb/extract/'+ scode + '?' + mindlablogin + '\\&projectCode=' + proj_code
    output = _wget_system_call(url)

    # Split at '\n'
    subj_list = output.split('\n')
    # Remove any empty entries!
    subj_list = [x for x in subj_list if x]
    
    return subj_list

def get_studies(mindlablogin,proj_code,subj_ID, modality=None):    
    
    url = 'http://hyades00.pet.auh.dk/modules/Ptdb/extract/studies?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID 
    output = _wget_system_call(url)

    # Split at '\n'
    stud_list = output.split('\n')
    # Remove any empty entries!
    stud_list = [x for x in stud_list if x]

    if modality:
        for study in stud_list:
            url = 'http://hyades00.pet.auh.dk/modules/Ptdb/extract/modalities?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID + '\\&study=' + study
            output = _wget_system_call(url).split('\n')[0]
            #print output, '==', modality
            if output == modality:
                return study
                ### NB!! This only matches first hit! If subject contains several studies with this modality, 
                ### only first one is returned... Fixme
    else:    
        return stud_list

def get_series(mindlablogin,proj_code,subj_ID, study, modality):    
    
    url = 'http://hyades00.pet.auh.dk/modules/Ptdb/extract/series?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID + '\\&study=' + study+ '\\&modality=' + modality
    output = _wget_system_call(url)

    # Split at '\n'
    series_list = output.split('\n')
    # Remove any empty entries!
    series_list = [x for x in series_list if x]

    return series_list


mindlablogin = 'templogin=272@db16e3511df0c8488b9e2b38e84ff121510cb665'
proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

subj_list=get_subjects(mindlablogin, proj_code, subj_type='included')
print "Subjects\n", subj_list

stud_list=get_studies(mindlablogin, proj_code, subj_ID=subj_list[0],modality='MEG')
print "Subject: ", subj_list[0], 'has this study with MEG modality present'
print stud_list

series_list = get_series(mindlablogin,proj_code,subj_ID=subj_list[0], study=stud_list, modality='MEG')
print "in which there are the following series:"
print series_list