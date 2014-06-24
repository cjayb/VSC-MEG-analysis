import sys
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
from access import Query
from analysis_dict import Anadict

db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
ad=Anadict(db)

import mne

filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 35.0, 'highpass': 0.5}

filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])

fwd_params = {'spacing': 'oct-6',
        'bem': '-5120-bem-sol.fif ',
        'others': ' --megonly --mindist 5 ',
        'force': True}

subj = '007_SGF'
evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
opr_path = ad._scratch_folder + '/operators/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj

session = 'pre'
trial_type = 'VS' # use the VS basline covariance

evo_file = evo_path + '/' + trial_type + '_' + session + '-avg.fif'
cov_file = evo_path + '/' + trial_type + '_' + session + '-cov.fif'
fwd_file = opr_path + '/' + trial_type + '_' + session + \
        '-' + fwd_params['spacing'] + '-fwd.fif'

# Load data
stdA = mne.read_evokeds(evo_file, condition='stdA')
devA = mne.read_evokeds(evo_file, condition='devA')
oddA = devA - stdA

Leff = 1./(1./stdA.nave + 1./devA.nave)

fwd_opr = mne.read_forward_solution(fwd_file, surf_ori=True)

noise_cov = mne.read_cov(cov_file)

# regularize noise covariance, assume nave not involved!
noise_cov = mne.cov.regularize(noise_cov, oddA.info,
        mag=0.05, grad=0.05, proj=True)

snr = 3.0
lambda2 = 1.0 / snr ** 2
method = 'MNE'

inv_opr = mne.minimum_norm.make_inverse_operator(oddA.info, 
        fwd_opr, noise_cov, loose=0.2, depth=0.8)
stc = mne.minimum_norm.apply_inverse(oddA, inv_opr, lambda2, method, pick_ori=None)

oddA.nave = Leff
inv_opr_Leff = mne.minimum_norm.make_inverse_operator(oddA.info, 
        fwd_opr, noise_cov, loose=0.2, depth=0.8)

stc_Leff = mne.minimum_norm.apply_inverse(oddA, inv_opr_Leff, lambda2, method, pick_ori=None)

brain = stc.plot(subject=subj, surface='inflated', hemi='lh', fmin=0., fmid=1e-10, fmax=3e-10)
brain.set_time(90)
brain_Leff = stc_Leff.plot(subject=subj, surface='inflated', hemi='lh', fmin=0., fmid=1e-10, fmax=3e-10)
brain_Leff.set_time(90)
