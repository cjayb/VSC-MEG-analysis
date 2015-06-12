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

epoch_params = {'rsl': 250, 'savgol_hf': 20.}

evoked_categories = dict(
        VS =  dict(face=(['stdB','devB'], ['stdA','devA']),
                 oddA1 =(['A1'],['stdA']),oddB1 =(['B1'],['stdB']),
                 oddA2 =(['A2'],['stdA']),oddB2 =(['B2'],['stdB']),
                 oddA3 =(['A3'],['stdA']),oddB3 =(['B3'],['stdB']),
                 oddA4 =(['A4'],['stdA']),oddB4 =(['B4'],['stdB']),
                 oddA5 =(['A5'],['stdA']),oddB5 =(['B5'],['stdB']),
                 oddA6 =(['A6'],['stdA']),oddB6 =(['B6'],['stdB']),
                 odd1  =(['A1','B1'],['stdA','stdB']),
                 odd2  =(['A2','B2'],['stdA','stdB']),
                 odd3  =(['A3','B3'],['stdA','stdB']),
                 odd4  =(['A4','B4'],['stdA','stdB']),
                 odd5  =(['A5','B5'],['stdA','stdB']),
                 odd6  =(['A6','B6'],['stdA','stdB']),
                 stdA=(['stdA'],),devA=(['devA'],),
                 stdB=(['stdB'],),devB=(['devB'],),
                 A1 =(['A1'],),B1 =(['B1'],),
                 A2 =(['A2'],),B2 =(['B2'],),
                 A3 =(['A3'],),B3 =(['B3'],),
                 A4 =(['A4'],),B4 =(['B4'],),
                 A5 =(['A5'],),B5 =(['B5'],),
                 A6 =(['A6'],),B6 =(['B6'],)
                 ),
        N2pc =  dict(
                 diff  =(['A1','A2','A3','B1','B2','B3'], ['A4','A5','A6','B4','B5','B6']),
                 diffA =(['A1','A2','A3'], ['A4','A5','A6']),
                 diffB =(['B1','B2','B3'], ['B4','B5','B6']),
                 devLH  =(['A1','A2','A3','B1','B2','B3'], ['stdA','stdB']),
                 devRH  =(['A4','A5','A6','B4','B5','B6'], ['stdA','stdB']),
                 ),
        FB =  dict(face=(['stdB','devB'], ['stdA','devA']),
                 odd =(['devA','devB'],  ['stdA','stdB']),
                 stdA=(['stdA'],),devA=(['devA'],),
                 stdB=(['stdB'],),devB=(['devB'],)
                 ),
        FFA = dict(diff=(['A','B'], ['blur']),
                 face=(['A','B'],),
                 blur=(['blur'],),
                 )
        )

fwd_params = {
        'spacing': 'ico5', # following Khan et al. (2013)
        #'spacing': 'oct-6',
        'bem': '-5120-bem-sol.fif ',
        'others': ' --megonly --mindist 5 ',
        'mindist': 5.,
        'force': True}

inv_params = dict(loose=0.2, depth=0.8,
        limit_depth_chs=True,
        fixed=False)

def mkdir_p(pth):

    try: 
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def file_exists(pth):
    return os.path.exists(pth)

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


# From J-R King, June 2015
# ad hoc: Scaled Logistic Regression with probabilistic output
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

class force_predict(object):
    def __init__(self, clf, mode='predict_proba', axis=0):
        self._mode = mode
        self._axis = axis
        self._clf = clf
 
    def fit(self, X, y, **kwargs):
        self._clf.fit(X, y, **kwargs)
        self._copyattr()
 
    def predict(self, X):
        if self._mode == 'predict_proba':
            return self._clf.predict_proba(X)[:, self._axis]
        elif self._mode == 'decision_function':
            distances = self._clf.decision_function(X)
            if len(distances.shape) > 1:
                return distances[:, self._axis]
            else:
                return distances
        else:
            return self._clf.predict(X)
 
    def get_params(self, deep=True):
        return dict(clf=self._clf, mode=self._mode, axis=self._axis)
 
    def _copyattr(self):
        for key, value in self._clf.__dict__.iteritems():
            self.__setattr__(key, value)

# Area Under the Curve Scorer:
def auc_scorer(y_true, y_pred):
    le = LabelBinarizer()
    y_true = le.fit_transform(y_true)
    return roc_auc_score(y_true, y_pred)

