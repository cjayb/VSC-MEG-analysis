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
from getpass import getuser, getpass
import os

class Query():
    """
    Query object for communicating with the MINDLab "DICOM" database
    """

    def __init__(self, proj_code, mindlablogin='~/.mindlabdblogin', username=None, verbose=None):

        try: 
            with open(os.path.expanduser(mindlablogin)):
                if verbose: print 'Reading login credentials from ' + mindlablogin
                f = open(os.path.expanduser(mindlablogin))
                self._login_code = f.readline()
                f.close()
        except IOError:
            print 'Login credentials not found, please enter them here'
            print 'WARNING: This might not work if you\'re in an IDE (e.g. spyder)!'
            if username:
                usr = username                
            else:
                usr = getuser()
            
            prompt = 'User \"%s\", please enter your password: ' % usr
            pwd = getpass(prompt)
            
            url = 'login/username/' + usr +  '/password/' + pwd
            output = _wget_system_call(url)
            if _wget_error_handling(output) < 0:
                sysexit(-1)
            else:
                self._login_code = output
                #mindlablogin='~/.mindlabdblogin'
                print "Code generated, writing to %s" % mindlablogin
                fout = open(os.path.expanduser(mindlablogin),'w')
                fout.write(self._login_code)
                fout.close()
                os.chmod(os.path.expanduser(mindlablogin), 0400)
            
                

def _wget_error_handling(stdout):
    # Do error checking of wget output: will start with "error:"
    if 'error' in stdout:
        print "Something is wrong, database answers as follows (dying...):"
        print stdout
        return -1
        
    return 0

def _wget_system_call(url):
    server = 'http://hyades00.pet.auh.dk/modules/Ptdb/extract/'
    cmd = 'wget -qO - test ' + server + url
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
    
    url = scode + '?' + mindlablogin + '\\&projectCode=' + proj_code
    output = _wget_system_call(url)

    # Split at '\n'
    subj_list = output.split('\n')
    # Remove any empty entries!
    subj_list = [x for x in subj_list if x]
    
    return subj_list

def get_studies(mindlablogin,proj_code,subj_ID, modality=None, unique=True):    
    
    url = 'studies?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID 
    output = _wget_system_call(url)

    # Split at '\n'
    stud_list = output.split('\n')
    # Remove any empty entries!
    stud_list = [x for x in stud_list if x]

    if modality:
        for study in stud_list:
            url = 'modalities?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID + '\\&study=' + study
            output = _wget_system_call(url).split('\n')
            #print output, '==', modality
            
            for entry in output:              
                if entry == modality:
                    if unique:
                        return study
                        ### NB!! This only matches first hit! If subject contains several studies with this modality, 
                        ### only first one is returned... Fixme
                    else:
                        # must re-write code a bit to accommodate the existence of
                        # several studies containing the desired modality...
                        print "Error: non-unique modalities not implemented yet!"
                        sysexit(-1)
    else:    
        return stud_list

def get_series(mindlablogin,proj_code,subj_ID, study, modality):    
    
    url = 'series?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID + '\\&study=' + study+ '\\&modality=' + modality
    output = _wget_system_call(url)

    # Split at '\n'
    series_list = output.split('\n')
    # Remove any empty entries!
    series_list = [x for x in series_list if x]

    # Series are in a list with [name, number], need to combine them to fit the 
    # actual DB structure (00N.name). @Lars: why does wget return this when 
    # what we actually want is the following?!...

    for ii,entry in enumerate(series_list):
        tmp = entry.split(' ')
        filenum = '%03d.' % int(tmp[1])
        series_list[ii] = filenum + tmp[0] 

    return series_list

def get_files(mindlablogin,proj_code,subj_ID, study, modality, series):    
    
    url = 'files?' + mindlablogin + '\\&projectCode=' + proj_code + '\\&subjectNo=' + subj_ID + '\\&study=' + study+ '\\&modality=' + modality+ '\\&serieNo=' + series
    output = _wget_system_call(url)

    # Split at '\n'
    file_list = output.split('\n')
    # Remove any empty entries!
    file_list = [x for x in file_list if x]

    return file_list

if __name__ == '__main__':
    
    #test code

    mindlablogin = 'templogin=272@db16e3511df0c8488b9e2b38e84ff121510cb665'
    proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'
    
    subj_list=get_subjects(mindlablogin, proj_code, subj_type='included')
    print "Subjects\n", subj_list
    
    stud_list=get_studies(mindlablogin, proj_code, subj_ID=subj_list[1],modality='MEG')
    print "Subject: ", subj_list[1], 'has this study with MEG modality present'
    print stud_list
    
    series_list = get_series(mindlablogin,proj_code,subj_ID=subj_list[1], study=stud_list, modality='MEG')
    print "in which there are the following series:"
    print series_list
    
    print 'The first series (' + series_list[0] + ') contains the files:'
    file_list = get_files(mindlablogin,proj_code,subj_ID=subj_list[1], study=stud_list, modality='MEG', series=series_list[0])
    print file_list

    Q = Query(proj_code=proj_code)