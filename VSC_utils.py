import sys, os, errno
import mne
import csv
import numpy as np

# this is the sss-run used as "raw"
input_files = 'tsss_initial'

# Set epoch parameters
tmin, tmax = -0.3, 0.4  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.15, 0.2  # reject trial only if blinks in the 300 ms middle portion!
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13) # compare to standard rejection
#reject = None
baseline = (-0.15, 0.)

# This defines what has been done in a scr_run-file filtering the tsss'd data
# should really go into a module, along with other defaults and a couple of
# utility functions (esp. mkdir_p)
filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 35.0, 'highpass': 0.5}

filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])

epoch_params = {'rsl': 250}

fwd_params = {'spacing': 'oct-6',
        'bem': '-5120-bem-sol.fif ',
        'others': ' --megonly --mindist 5 ',
        'force': True}


def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def split_events_by_trialtype(events, condition='VS'):
    if 'VS' in condition:
        devsA, devsB = range(111,117), range(211,217)
        VS_eve = mne.pick_events(events, include=range(100,220))
        VS_eve = mne.merge_events(VS_eve, [100], 10, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, [200], 20, replace_events=True)
        # Don't replace the deviants, make a copy instead!
        VS_eve = mne.merge_events(VS_eve, devsA, 11, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, devsB, 21, replace_events=True)

        # This hack is needed to get both 11/21's and 11N/21N's together!
        tmp = mne.pick_events(events, include=devsA+devsB)
        #tmp[:,0] += 1 # add a ms
        VS_eve = np.concatenate((VS_eve, tmp), axis=0)
        VS_eve = VS_eve[np.argsort(VS_eve[:, 0])]

        FB_eve = mne.pick_events(events, include=range(10,22))
        
        eve_dict = dict(VS=VS_eve, FB=FB_eve)

    elif 'FFA' in condition:
        FFA_eve = mne.pick_events(events, include=[100, 150, 200])
        eve_dict = dict(FFA=FFA_eve)

    id_dict = dict(VS=dict(stdA=10, stdB=20, devA=11, devB=21,
            A1=111, A2=112,A3=113,A4=114,A5=115,A6=116,
            B1=211, B2=212,B3=213,B4=214,B5=215,B6=216),
            FB=dict(stdA=10, stdB=20, devA=11, devB=21),
            FFA=dict(A=100, B=200, blur=150))

    return eve_dict, id_dict

def load_excludes(ica_excludes_folder, subj, cond):
    pth = ica_excludes_folder + '/' + subj + '.csv'

    with open(pth, 'rb') as csvfile:
        exreader = csv.reader(csvfile, delimiter=',')
        hdr = exreader.next()
        try:
            colind = hdr.index(cond)
        except ValueError:
            print 'condition must be VS1, VS2 or FFA!'
            raise ValueError

        ica_excludes = []
        for row in exreader:
            ica_excludes += row[colind].split('|')

    # remove emptys
    ica_excludes = filter(len, ica_excludes)
    return map(int,ica_excludes)

