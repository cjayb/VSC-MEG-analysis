# After initial round of analysis in Matlab (Fieldtrip, J-R King), most sensor-
# level results were replicated in python.
#   * there were differences in the csXoddXsession interaction, even though
#     the csXodd's for pre and post seemed similar...?
#   * JRK used "robust averaging", whereas in python we reject

# Eye movements are likely to still be an issue in this dataset, at least after
# about 150 msec.
#   * the events for stdA, devA, stdB, devB were simply rejected on the basis
#     of large amplitudes; this will not remove saccades from the data,
#   * I think we have to do some SSP or ICA-based cleaning up of the data
#     (frontal sources seen a lot in the source estimates)
# The interaction contrasts didn't look very inspiring in source space
#   * might want to play around with the SNR and nave-parameters, but unless
#     there is good reason to believe there's something there (sensor-level),
#     it's probably a waste of time.
#
# License: BSD (3-clause)
import matplotlib
matplotlib.use('agg') # force non-interactive plotting
import numpy as np
import os, errno
import json
from copy import deepcopy
from operator import add # for stc reduction operation

machine_name = os.uname()[1].split('.')[0]

if 'isis' in machine_name:
    import sys
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
    from access import Query
    from analysis_dict import Anadict

    db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
    ad=Anadict(db)
elif 'mba-cjb' in machine_name or 'horus' in machine_name:
    class local_Anadict():
        def __init__(self):
            self._scratch_folder = '/Users/cjb/tmp/VSC-scratch'

    ad = local_Anadict()

import mne
try:
    from mne.io import Raw, pick_types
    from mne import read_evoked
except:
    from mne.fiff import Raw, pick_types, read_evoked

from viz_cjb import plot_evoked_topomap

do_epoching = False
do_evokeds = False
do_forward_solutions_evoked = False
do_inverse_operators_evoked = False

# perhaps doing sensor-level contrasts is a bad idea after all:
# don't we risk source cancellation when e.g. doing "odd"?
# Then we're localizing differences...
do_evokeds_to_source_estimates = False

do_morph_evokedSEs_to_fsaverage = False
do_average_morphed_evokedSEs = False
do_grandaverage_CScontrasts = False

do_sourcelevel_rmanova_stclustering = True

# These are obsolete since do_evokeds_to_source_estimates will do the 
# simple contrasts as well
do_source_level_contrasts = False
do_morph_contrasts_to_fsaverage = False
do_average_morph_maps = False

# These are mainly to allow comparison to Matlab code
do_sensor_level_contrasts = False
do_sensor_level_contrast_images_across = False

do_sensor_level_contrast_to_sources = False
do_sensor_level_contrast_to_sources_to_fsaverage = False

def mkdir_p(pth):

    try:
        os.makedirs(pth)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(pth):
            pass
        else:
            raise

def create_trigger_logic(cond):

    trig_mult_factor = {'CS+': 0, 'CS-': 0}
    if '1a' in cond:
        trig_mult_factor['CS+'] = 1
        trig_mult_factor['CS-'] = 2
    elif '1b' in cond:
        trig_mult_factor['CS+'] = 2
        trig_mult_factor['CS-'] = 1


    trig_logic = {}

    for cond in ['CS+', 'CS-']:
        trig_logic.update({cond: {}})
        for session in ['pre', 'post']:
            trig_logic[cond].update({session: {}})
            for imageType in ['VS','FB']:
                trig_logic[cond][session].update({imageType: {}})
                if imageType == 'VS':
                    trig_base = trig_mult_factor[cond]*100
                elif imageType == 'FB':
                    trig_base = trig_mult_factor[cond]*10
                for emotion in ['Neu','Ang']:
                    if emotion == 'Neu':
                        trig_code = trig_base
                    elif emotion == 'Ang':
                        # This is if we want to treat the position of the target separately!
                        # if imageType == 'VS':
                        #     trig_code = trig_base + np.arange(11,17)
                        # elif imageType == 'FB':
                        #     trig_code = trig_base + 1

                        # This way all Angry faces in the VS cond become 101/201...
                        trig_code = trig_base + 1

                    trig_logic[cond][session][imageType].update({emotion: trig_code})

    return trig_logic

def contrast_logic():
    devsA, devsB = range(111,117), range(211,217)
    contrasts = ['VS','FB']
    clogic = dict(VS={}, FB={})
    clogic['VS'].update({'all': dict(all=[100,200] + devsA + devsB)})
    clogic['VS'].update({'face': dict(faceA=[100] + devsA, faceB=[200] + devsB)})
    clogic['VS'].update({'odd': dict(std=[100,200],dev = devsA + devsB)})
    clogic['FB'].update({'all': dict(all=[10, 20, 11, 21])})
    clogic['FB'].update({'face': dict(faceA=[10, 11], faceB=[20, 21])})
    clogic['FB'].update({'odd': dict(std=[10,20],dev = [11, 21])})

    return clogic

def events_logic(events, contrast):
    devsA, devsB = range(111,117), range(211,217)
    con_dict = dict(VS={}, FB={})
    VS_eve = mne.pick_events(events, include=range(100,220))
    FB_eve = mne.pick_events(events, include=range(10,22))
    if 'face' == contrast:
        VS_eve = mne.merge_events(VS_eve, [100]+devsA, 1, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, [200]+devsB, 2, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [10,11], 1, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [20,21], 2, replace_events=True)

        con_dict = dict(faceA=1, faceB=2)
        con_names = ['faceB','faceA'] # [0] - [1]

    elif 'odd' == contrast:
        VS_eve = mne.merge_events(VS_eve, [100,200], 1, replace_events=True)
        VS_eve = mne.merge_events(VS_eve, devsA + devsB, 2, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [10,20], 1, replace_events=True)
        FB_eve = mne.merge_events(FB_eve, [11,21], 2, replace_events=True)

        con_dict = dict(std=1, dev=2)
        con_names = ['dev','std'] # [0] - [1]

    eve_dict = dict(VS=VS_eve, FB=FB_eve)
    return eve_dict, con_dict, con_names

def split_events_by_trialtype(events):
    devsA, devsB = range(111,117), range(211,217)
    VS_eve = mne.pick_events(events, include=range(100,220))
    VS_eve = mne.merge_events(VS_eve, [100], 10, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, [200], 20, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, devsA, 11, replace_events=True)
    VS_eve = mne.merge_events(VS_eve, devsB, 21, replace_events=True)
    FB_eve = mne.pick_events(events, include=range(10,22))

    eve_dict = dict(VS=VS_eve, FB=FB_eve)
    id_dict = dict(stdA=10, stdB=20, devA=11, devB=21) # now same for VS and FB

    return eve_dict, id_dict
###############################################################################
# Set epoch parameters
tmin, tmax = -0.4, 0.6  # no need to take more than this, wide enough to see eyemov though
rej_tmin, rej_tmax = -0.2, 0.2  # reject trial only if blinks in the 400 ms middle portion!
baseline = (-0.2, 0.)
reject = dict(eog=150e-6, mag=4e-12, grad=4000e-13)
rsl_fs = 250 # Resample epochs
filter_params = {'input_files': 'tsss_initial',
                 'lowpass': 35.0, 'highpass': 0.5}

filt_dir = '%.1f-%.1fHz' % (filter_params['highpass'], filter_params['lowpass'])

fwd_params = {'spacing': 'oct-6',
        'bem': '-5120-bem-sol.fif ',
        'others': ' --megonly --mindist 5 ',
        'force': True}

if do_epoching:
    for subj in ad.analysis_dict.keys():

        # Drop the FFA session for now, deal with it separately, also empty room
        session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                if ('FFA' not in x and 'empty' not in x)]

        raw_path = ad._scratch_folder + '/filtered/' + filter_params['input_files'] + '/' + filt_dir + '/' + subj
        eve_path = ad._scratch_folder + '/events.fif/' + subj + '/raw'

        epo_path = ad._scratch_folder + '/epochs/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        mkdir_p(epo_path)

        for sesname in session_names:

            if '_1' == sesname[-2:]:
                session = 'pre'
            elif '_2' == sesname[-2:]:
                session = 'post'

            fname = raw_path + '/' + sesname + '_filt.fif'
            raw = Raw(fname, preload=False)
            events = mne.read_events(eve_path + '/' + sesname + '-eve.fif')
            picks = pick_types(raw.info, meg=True, eeg=False, stim=True, eog=True, misc=True)
            eve_dict, id_dict = split_events_by_trialtype(events)
            for trial_type in ['VS','FB']:

                print('Extracting %s (%s) epochs for %s' % (trial_type, session, subj))
                epochs = mne.Epochs(raw, eve_dict[trial_type], id_dict,
                                    tmin, tmax, picks=picks, verbose=False,
                                    baseline=baseline, reject=reject, preload=True,
                                    reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings
                # Check the drop_log for these preload'ed epochs: does the drop
                # log indicate the dropped epochs, can they be un-dropped after the fact?
                # Do we in fact have to actively drop them, despite preloading?

                print('Resampling...')
                epochs.resample(rsl_fs, n_jobs=6, verbose=False) # Trust the defaults here

                epo_out = epo_path + '/' + trial_type + '_' + session + '-epo.fif'
                epochs.save(epo_out)  # save epochs to disk


if do_evokeds: # do a couple of "main effects"

    import matplotlib.pyplot as plt
    clim_all = dict(mag=[-400, 400], grad=[0, 80])
    clim_con = dict(mag=[-125, 125], grad=[0, 25])
    topo_times = np.arange(0.0, 0.210,0.020)


    evoked_categories = dict(all=(['stdA','stdB','devA','devB'], ),
            face=(['stdB','devB'], ['stdA','devA']),
            odd=(['devA','devB'],  ['stdA','stdB']),
            face_std=(['stdB'], ['stdA']),
            face_dev=(['devB'], ['devA']),
            odd_A=(['devA'], ['stdA']),
            odd_B=(['devB'], ['stdB']),
            stdA=(['stdA'],),devA=(['devA'],),stdB=(['stdB'],),devB=(['devB'],))
#                for categ in evoked_categories:
#                    evoked_cur = epochs[categ].average()
#                    evo_out = evo_path + '/' + trial_type + '_' + session + '_' + categ + '-ave.fif'
#                    evoked_cur.save(evo_out)

    #for subj in ['007_SGF']:
    for subj in ad.analysis_dict.keys():

        epo_path = ad._scratch_folder + '/epochs/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        img_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj + '/img'
        mkdir_p(img_path) # evo_path is also written

        for session in ['pre','post']:
            for trial_type in ['VS','FB']:
                fname = epo_path + '/' + trial_type + '_' + session + '-epo.fif'
                epochs = mne.read_epochs(fname)

                evokeds = []
                for categ in evoked_categories.keys():
                    if len(evoked_categories[categ]) == 2:
                        evokeds.append(epochs[evoked_categories[categ][0]].average() - \
                                epochs[evoked_categories[categ][1]].average())
                        evokeds[-1].comment = categ
                    else:
                        evokeds.append(epochs[evoked_categories[categ][0]].average())
                        evokeds[-1].comment = categ

                cov_all = mne.compute_covariance(epochs, tmin=baseline[0], tmax=baseline[1]) # same covariance for all contrasts
                figs = mne.viz.plot_cov(cov_all, epochs.info, show=False)
                figs[0].savefig(img_path + '/' + trial_type + '_' + session + '_all_covMAT.png')
                figs[1].savefig(img_path + '/' + trial_type + '_' + session + '_all_covSVD.png')

                for e in evokeds:
                    if e.comment == 'all':
                        e.plot_image(clim=clim_all, show=False)
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_allERF_time.png')
                        plot_evoked_topomap(e,topo_times, show=False, vmin=[clim_all['grad'][0],clim_all['mag'][0]], vmax=[clim_all['grad'][1],clim_all['mag'][1]])
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_allERF_topo.png')
                    elif e.comment == 'face':
                        e.plot_image(clim=clim_con, show=False)
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_faceERF_time.png')
                        plot_evoked_topomap(e,topo_times, show=False, vmin=[clim_con['grad'][0],clim_con['mag'][0]], vmax=[clim_con['grad'][1],clim_con['mag'][1]])
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_faceERF_topo.png')
                    elif e.comment == 'odd':
                        e.plot_image(clim=clim_con, show=False)
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_oddERF_time.png')
                        plot_evoked_topomap(e,topo_times, show=False, vmin=[clim_con['grad'][0],clim_con['mag'][0]], vmax=[clim_con['grad'][1],clim_con['mag'][1]])
                        plt.savefig(img_path + '/' + trial_type + '_' + session + '_oddERF_topo.png')

                evo_out= evo_path + '/' + trial_type + '_' + session + '-avg.fif'
                mne.write_evokeds(evo_out, evokeds)  # save evoked data to disk

                cov_out = evo_path + '/' + trial_type + '_' + session + '-cov.fif'
                print cov_out
                cov_all.save(cov_out)  # save covariance data to disk

if do_forward_solutions_evoked:

    # check that 'T1' is attached to subject first, assume then MR preproc OK
    for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        fwd_cmd = 'mne_do_forward_solution'
        if fwd_params['force']:
            fwd_cmd += ' --overwrite '
        fwd_cmd += ' --subject ' + subj
        fwd_cmd += ' --spacing ' + fwd_params['spacing']
        fwd_cmd += ' --bem ' + subj + fwd_params['bem']
        fwd_cmd += fwd_params['others']

        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        fwd_path = ad._scratch_folder + '/operators/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        mkdir_p(fwd_path)

        trans_fif = ad._scratch_folder + '/trans/' + subj + '-trans.fif'
        fwd_cmd += ' --mri ' + trans_fif

        for session in ['pre','post']:
            for trial_type in ['VS']: # only take one, the FWD model only depends on
                                      # any possible projections applied
                evo_file = evo_path + '/' + trial_type + '_' + session + '-avg.fif'

                cmd = fwd_cmd
                # A bit hairy: since all the categories will have the same SSP
                # etc applied, we can make a single fwd operator for all
                cmd += ' --meas ' + evo_file # contains all the evokeds
                fwd_out = fwd_path + '/' + trial_type + '_' + session + \
                                    '-' + fwd_params['spacing'] + '-fwd.fif'
                cmd += ' --fwd ' + fwd_out
                print cmd
                st = os.system(cmd)
                if st != 0:
                    raise RuntimeError('mne_do_forward_solution returned with error %d' % st)

                # create a link for the FB trials
                fwd_link = fwd_path + '/FB_' + session + \
                                    '-' + fwd_params['spacing'] + '-fwd.fif'
                os.symlink(fwd_out, fwd_link)


if do_inverse_operators_evoked:
    # We'll use the covariance estimate from the VS baseline for both VS and FB
    # This means that the inverse operators will be identical (just like the fwd)

    # check that 'T1' is attached to subject first, assume then MR preproc OK
    for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        opr_path = ad._scratch_folder + '/operators/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj

        for session in ['pre','post']:

            trial_type = 'VS' # use the VS basline covariance

            evo_file = evo_path + '/' + trial_type + '_' + session + '-avg.fif'
            cov_file = evo_path + '/' + trial_type + '_' + session + '-cov.fif'
            fwd_file = opr_path + '/' + trial_type + '_' + session + \
                    '-' + fwd_params['spacing'] + '-fwd.fif'
            inv_file = opr_path + '/' + trial_type + '_' + session + \
                    '-' + fwd_params['spacing'] + '-inv.fif'
            inv_link = opr_path + '/FB_' + session + \
                    '-' + fwd_params['spacing'] + '-inv.fif'

            # Load data
            evoked = mne.read_evokeds(evo_file, condition='all')
            fwd_opr = mne.read_forward_solution(fwd_file, surf_ori=True)
            noise_cov = mne.read_cov(cov_file)

            # regularize noise covariance
            noise_cov = mne.cov.regularize(noise_cov, evoked.info,
                    mag=0.05, grad=0.05, proj=True)

            inv_opr = mne.minimum_norm.make_inverse_operator(evoked.info,
                    fwd_opr, noise_cov, loose=0.2, depth=0.8)

            mne.minimum_norm.write_inverse_operator(inv_file, inv_opr)
            os.symlink(inv_file, inv_link)

if do_evokeds_to_source_estimates:

    # Project to source space, while cropping to [-100, 200]
#    snr = 3.0
#    lambda2 = 1.0 / snr ** 2
    time_range = (-0.100, 0.200)
    methods = ['MNE','dSPM']
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)

    do_evoked_contrasts = {'stdA': True,'stdB': True,'devA': True,'devB': True,
            'all': True, 'face': True, 'odd': True,
            'face_std': True, 'face_dev': True, 'odd_A': True, 'odd_B': True}

    SNRs = {'stdA': 3.,'stdB': 3.,'devA': 3.,'devB': 3., 'all': 3.,
            'face': 1., 'odd': 1., 'face_std': 1., 'face_dev': 1., 'odd_A': 1., 'odd_B': 1.}

    # check that 'T1' is attached to subject first, assume then MR preproc OK
    for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        opr_path = ad._scratch_folder + '/operators/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        mkdir_p(stc_path)

        for session in ['pre','post']:
            for trial_type in ['VS','FB']:
                evo_file = evo_path + '/' + trial_type + '_' + session + '-avg.fif'
                inv_file = opr_path + '/' + trial_type + '_' + session + \
                        '-' + fwd_params['spacing'] + '-inv.fif'

                print 20*'#'
                print 'Doing %s -> %s_%s...' % (subj, session, trial_type)
                print 20*'#'

                print 'Loading inverse operator...'
                inv_opr = mne.minimum_norm.read_inverse_operator(inv_file, verbose=False)
                for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
                    print 'Applying methods to condition', cond
                    evoked = mne.read_evokeds(evo_file, condition=cond, verbose=False)
                    lambda2 = 1. / SNRs[cond] ** 2.
                    for method in methods:
                        #print 'Applying inverse with method:', method
                        stc = mne.minimum_norm.apply_inverse(evoked, inv_opr,
                                lambda2, method, pick_ori=ori_sel, verbose=False)
                        stc.crop(tmin=time_range[0], tmax=time_range[1]) # CROP

                        # Save result in stc files
                        stc_file = stc_path + '/' + trial_type + '_' + session + \
                                '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                        stc.save(stc_file, verbose=False)

if do_morph_evokedSEs_to_fsaverage:

    do_evoked_contrasts = {'stdA': True,'stdB': True,'devA': True,'devB': True,
            'all': False, 'face': False, 'odd': False,
            'face_std': False, 'face_dev': False, 'odd_A': False, 'odd_B': False}
    methods = ['MNE','dSPM']
    trial_types = ['VS', 'FB']

    # This seems very hacky, but will have to try to under
    # stand later...
    #vertices_to = [np.arange(10242), np.arange(10242)]
    #subject_to = 'fsaverage'
    vertices_to = [np.arange(642), np.arange(642)] # given the below, this might as well say "3" (for grade), right?
    subject_to = 'fsaverage3' #this has to be the fsa3 morph, using the fsa gives doesn't seem to work...
    # Allowed values are: 2 (162 locations), 3 (642 locations), and 4 (2562 locations)

    for subject_from in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subject_from
        morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to
        mkdir_p(morph_stc_path)

        for trial_type in trial_types:
            print 20*'#'
            print '%s -> %s' % (subject_from, trial_type)
            print 20*'#'
            for method in methods:
                for session in ['pre','post']:
                    for key in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
                        stc_file = stc_path + '/' + trial_type + '_' + session + '-' + \
                                fwd_params['spacing'] + '_' + key + '_' + method
                        stc_from = mne.read_source_estimate(stc_file)
                        print 'Morphing to', subject_to
                        stc_to = mne.morph_data(subject_from, subject_to,
                                stc_from, grade=vertices_to, n_jobs=4, verbose=False)

                        stc_file = morph_stc_path + '/' + subject_from + '_' + trial_type + \
                                '_' + session + '-' + fwd_params['spacing'] + '_' + key + '_' + method
                        stc_to.save(stc_file, verbose=False)

if do_average_morphed_evokedSEs:

    trial_types = ['VS', 'FB']
    methods = ['MNE','dSPM']
    
    subject_to = 'fsaverage3'

    do_evoked_contrasts = {'stdA': True,'stdB': True,'devA': True,'devB': True,
            'all': False, 'face': False, 'odd': False,
            'face_std': False, 'face_dev': False, 'odd_A': False, 'odd_B': False}

    morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to
    for trial_type in trial_types:
        for method in methods:
            for key in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
                stc_ave = {'pre': [], 'post': []}
                for session in ['pre','post']:
                    stc_list = []
                    for subject_from in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:
                        stc_file = morph_stc_path + '/' + subject_from + '_' + trial_type + \
                                '_' + session + '-' + fwd_params['spacing'] + '_' + key + '_' + method
                        stc_list.append(mne.read_source_estimate(stc_file))

                    stc_ave[session] = reduce(add, stc_list)
                    stc_ave[session] /= len(stc_list)
                    #stc_ave[session] = stc_list[0]
                    #for ii in range(1,len(stc_list)):
                    #    stc_ave[session] = (float(ii)*stc_ave[session] + stc_list[ii]) / (float(ii)+1.)

                    stc_file = morph_stc_path + '/AVG_' + trial_type + '_' + session + \
                            '-' + fwd_params['spacing'] + '_' + key + '_' + method
                    print 'Saving average %s -> %s -> %s' % (trial_type, key, method)
                    stc_ave[session].save(stc_file, verbose=False)

                stc_ave_avg = ( stc_ave['post'] + stc_ave['pre'] ) / 2.0
                stc_ave_dif = stc_ave['post'] - stc_ave['pre']

                stc_file = morph_stc_path + '/AVG_' + trial_type + '_both' + \
                        '-' + fwd_params['spacing'] + '_' + key + '_' + method
                stc_ave_avg.save(stc_file, verbose=False)
                stc_file = morph_stc_path + '/AVG_' + trial_type + '_diff' + \
                        '-' + fwd_params['spacing'] + '_' + key + '_' + method
                stc_ave_dif.save(stc_file, verbose=False)


if do_grandaverage_CScontrasts:

    trial_types = ['VS', 'FB']
    methods = ['MNE','dSPM']

    subject_to = 'fsaverage3'

    # This isn't doing anything...
    # do_CScontrasts = {'oddXsession': True}

    morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to
    mkdir_p(morph_stc_path + '/CS')
    for trial_type in trial_types:
        for method in methods:
            stc_CSp_diff = []
            stc_CSm_diff = []
            stc_interaction = []
            for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:
                # Drop the FFA session for now, deal with it separately, also empty room
                session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                        if ('FFA' not in x and 'empty' not in x)]

                if '1a' in session_names[0]:
                    CScode = {'CS+': 'A', 'CS-': 'B'}
                elif '1b' in session_names[0]:
                    CScode = {'CS+': 'B', 'CS-': 'A'}

                stc_file = morph_stc_path + '/' + subj + '_' + trial_type + '_post' + \
                        '-' + fwd_params['spacing'] + '_%s_' + method
                stc_CSp_post = mne.read_source_estimate(stc_file % ('dev'+CScode['CS+'])) - \
                        mne.read_source_estimate(stc_file % ('std'+CScode['CS+']))
                stc_CSm_post = mne.read_source_estimate(stc_file % ('dev'+CScode['CS-'])) - \
                        mne.read_source_estimate(stc_file % ('std'+CScode['CS-']))

                stc_file = morph_stc_path + '/' + subj + '_' + trial_type + '_pre' + \
                        '-' + fwd_params['spacing'] + '_%s_' + method
                stc_CSp_pre = mne.read_source_estimate(stc_file % ('dev'+CScode['CS+'])) - \
                        mne.read_source_estimate(stc_file % ('std'+CScode['CS+']))
                stc_CSm_pre = mne.read_source_estimate(stc_file % ('dev'+CScode['CS-'])) - \
                        mne.read_source_estimate(stc_file % ('std'+CScode['CS-']))

                stc_CSp_diff.append(stc_CSp_post - stc_CSp_pre)
                stc_CSm_diff.append(stc_CSm_post - stc_CSm_pre)
                stc_interaction.append(stc_CSp_diff[-1] - stc_CSm_diff[-1])

            # math out the grand average
            stc_ave = reduce(add, stc_CSp_diff)
            stc_ave /= len(stc_CSp_diff)
            #stc_ave = stc_CSp_diff[0]
            #for ii in range(1,len(stc_CSp_diff)):
            #    stc_ave = (float(ii)*stc_ave + stc_CSp_diff[ii]) / (float(ii)+1.)
            stc_file = morph_stc_path + '/CS/' + trial_type + \
                    '-' + fwd_params['spacing'] + '_odd_CSp_diff_' + method
            stc_ave.save(stc_file, verbose=False)

            stc_ave = reduce(add, stc_CSm_diff)
            stc_ave /= len(stc_CSm_diff)
            #stc_ave = stc_CSm_diff[0]
            #for ii in range(1,len(stc_CSm_diff)):
            #    stc_ave = (float(ii)*stc_ave + stc_CSm_diff[ii]) / (float(ii)+1.)
            stc_file = morph_stc_path + '/CS/' + trial_type + \
                    '-' + fwd_params['spacing'] + '_odd_CSm_diff_' + method
            stc_ave.save(stc_file, verbose=False)

            stc_ave = reduce(add, stc_interaction)
            stc_ave /= len(stc_interaction)
            #stc_ave = stc_interaction[0]
            #for ii in range(1,len(stc_interaction)):
            #    stc_ave = (float(ii)*stc_ave + stc_interaction[ii]) / (float(ii)+1.)
            stc_file = morph_stc_path + '/CS/' + trial_type + \
                    '-' + fwd_params['spacing'] + '_oddXsession_' + method
            stc_ave.save(stc_file, verbose=False)

if do_sourcelevel_rmanova_stclustering:

    from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                             grade_to_tris)
    from mne.stats import (spatio_temporal_cluster_test, f_threshold_twoway_rm,
                                   f_twoway_rm, summarize_clusters_stc)

    trial_types = ['FB']
    methods = ['MNE']

    subject_to = 'fsaverage3'
    morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to

    for trial_type in trial_types:
        for method in methods:
            conditions = []
            for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

                session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                        if ('FFA' not in x and 'empty' not in x)]

                if '1a' in session_names[0]:
                    CScode = {'CS+': 'A', 'CS-': 'B'}
                elif '1b' in session_names[0]:
                    CScode = {'CS+': 'B', 'CS-': 'A'}

                #stc_file = morph_stc_path + '/' + subj + '_' + trial_type + '_post' + \
                #        '-' + fwd_params['spacing'] + '_%s_' + method
                #stc_CSp_post = mne.read_source_estimate(stc_file % ('dev'+CScode['CS+'])) - \
                #        mne.read_source_estimate(stc_file % ('std'+CScode['CS+']))
                #stc_CSm_post = mne.read_source_estimate(stc_file % ('dev'+CScode['CS-'])) - \
                #        mne.read_source_estimate(stc_file % ('std'+CScode['CS-']))

                #stc_file = morph_stc_path + '/' + subj + '_' + trial_type + '_pre' + \
                #        '-' + fwd_params['spacing'] + '_%s_' + method
                #stc_CSp_pre = mne.read_source_estimate(stc_file % ('dev'+CScode['CS+'])) - \
                #        mne.read_source_estimate(stc_file % ('std'+CScode['CS+']))
                #stc_CSm_pre = mne.read_source_estimate(stc_file % ('dev'+CScode['CS-'])) - \
                #        mne.read_source_estimate(stc_file % ('std'+CScode['CS-']))

                #conditions.append([stc_CSp_post.crop(0, None), stc_CSm_post.crop(0, None),
                #        stc_CSp_pre.crop(0, None), stc_CSm_pre.crop(0, None)])

                stc_file = morph_stc_path + '/' + subj + '_' + trial_type + '_pre' + \
                        '-' + fwd_params['spacing'] + '_%s_' + method
                conditions.append(# factor A: dev vs std# factor B: identity
                        [mne.read_source_estimate(stc_file % 'devA').crop(0.05,None), 
                            mne.read_source_estimate(stc_file % 'devB').crop(0.05,None),
                            mne.read_source_estimate(stc_file % 'stdA').crop(0.05,None),
                            mne.read_source_estimate(stc_file % 'stdB').crop(0.05,None)]) 

            # we'll only consider the right hemisphere in this example.
            n_vertices_sample, n_times = conditions[0][0].rh_data.shape
            n_subjects = len(conditions)
            tmin = conditions[0][0].tmin
            tstep = conditions[0][0].tstep

            X = np.empty((n_vertices_sample, n_times, n_subjects, 4))
            for jj, subj in enumerate(conditions):
                for ii, condition in enumerate(subj):
                    X[:, :, jj, ii] = condition.rh_data[:, :]

            #    Now we need to prepare the group matrix for the ANOVA statistic.
            #    To make the clustering function work correctly with the
            #    ANOVA function X needs to be a list of multi-dimensional arrays
            #    (one per condition) of shape: samples (subjects) x time x space

            X = np.transpose(X, [2, 1, 0, 3])  # First we permute dimensions
            # finally we split the array into a list a list of conditions
            # and discard the empty dimension resulting from the split using numpy squeeze
            X = [np.squeeze(x) for x in np.split(X, 4, axis=-1)]

            ###############################################################################
            # Prepare function for arbitrary contrast

            # As our ANOVA function is a multi-purpose tool we need to apply a few
            # modifications to integrate it with the clustering function. This
            # includes reshaping data, setting default arguments and processing
            # the return values. For this reason we'll write a tiny dummy function.

            # We will tell the ANOVA how to interpret the data matrix in terms of
            # factors. This is done via the factor levels argument which is a list
            # of the number factor levels for each factor.
            factor_levels = [2, 2]

            # Finally we will pick the interaction effect by passing 'A:B'.
            # (this notation is borrowed from the R formula language)
            effects = 'A:B'  # Without this also the main effects will be returned.
            # Tell the ANOVA not to compute p-values which we don't need for clustering
            return_pvals = False

            # a few more convenient bindings
            n_times = X[0].shape[1]
            n_conditions = 4


            # A stat_fun must deal with a variable number of input arguments.
            def stat_fun(*args):
                # Inside the clustering function each condition will be passed as
                # flattened array, necessitated by the clustering procedure.
                # The ANOVA however expects an input array of dimensions:
                # subjects X conditions X observations (optional).
                # The following expression catches the list input, swaps the first and the
                # second dimension and puts the remaining observations in the third
                # dimension.
                data = np.squeeze(np.swapaxes(np.array(args), 1, 0))
                data = data.reshape(n_subjects, n_conditions,  # generalized if buffer used
                data.size / (n_subjects * n_conditions))
                return f_twoway_rm(data, factor_levels=factor_levels, effects=effects,
                        return_pvals=return_pvals)[0]
                #  drop p-values (empty array).
                # Note. for further details on this ANOVA function consider the
                # corresponding time frequency example.

            ###############################################################################
            # Compute clustering statistic
            #    To use an algorithm optimized for spatio-temporal clustering, we
            #    just pass the spatial connectivity matrix (instead of spatio-temporal)

            subject_to = 'fsaverage3'
            source_space = grade_to_tris(3)
            magic_vertno = 642 # or 642 or 10242
            # as we only have one hemisphere we need only need half the connectivity
            #rh_source_space = source_space[source_space[:, 0] < 10242]
            rh_source_space = source_space[source_space[:, 0] < magic_vertno]
            fsave_vertices = [np.array([]), np.arange(magic_vertno)]  # left hemisphere is empty
            n_vertices_fsave = fsave_vertices[1].shape

            print('Computing connectivity.')
            connectivity = spatial_tris_connectivity(rh_source_space)

            #    Now let's actually do the clustering. Please relax, on a small
            #    notebook and one single thread only this will take a couple of minutes ...
            pthresh = 0.0005
            f_thresh = f_threshold_twoway_rm(n_subjects, factor_levels, effects, pthresh)

            #    To speed things up a bit we will ...
            n_permutations = 100  # ... run fewer permutations (reduces sensitivity)

            print('Clustering.')
            T_obs, clusters, cluster_p_values, H0 = clu = \
                        spatio_temporal_cluster_test(X, connectivity=connectivity, n_jobs=6,
                                threshold=f_thresh, stat_fun=stat_fun,
                                n_permutations=n_permutations,
                                buffer_size=None)
            #    Now select the clusters that are sig. at p < 0.05 (note that this value
            #    is multiple-comparisons corrected).
            good_cluster_inds = np.where(cluster_p_values < 0.05)[0]

            print('Visualizing clusters.')

            #    Now let's build a convenient representation of each cluster, where each
            #    cluster becomes a "time point" in the SourceEstimate
            stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=tstep,
                    vertno=fsave_vertices,
                    subject=subject_to)

            #    Let's actually plot the first "time point" in the SourceEstimate, which
            #    shows all the clusters, weighted by duration

            # The brighter the color, the stronger the interaction between
            # stimulus modality and stimulus location

            brain = stc_all_cluster_vis.plot(subject_to, 'inflated', 'rh',
                    time_label='Duration significant (ms)', time_viewer=True)

            brain.set_data_time_index(0)
            brain.scale_data_colormap(fmin=0, fmid=10, fmax=20, transparent=True)
            brain.show_view('lateral')
            #brain.save_image('cluster-rh.png')
            brain.show_view('medial')


#########################
# Don't do these, use the evokeds as long as possible!
if do_source_level_contrasts:

    methods = ['MNE','dSPM']
    conds = ['stdA', 'stdB', 'devA', 'devB']
    do_contrasts = {'odd_pre': True, 'oddXsession': True, 'csXsession_dev': True, 'csXoddXsession': True}

    for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        stc_path_out = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj + '/sourceCon'
        mkdir_p(stc_path_out)

        # Drop the FFA session for now, deal with it separately, also empty room
        session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                if ('FFA' not in x and 'empty' not in x)]

        if '1a' in session_names[0]:
            CS_cond = dict(stdCSp='stdA',stdCSm='stdB',devCSp='devA',devCSm='devB')
        elif '1b' in session_names[0]:
            CS_cond = dict(stdCSp='stdB',stdCSm='stdA',devCSp='devB',devCSm='devA')

        for trial_type in ['VS','FB']:
            stc = dict(pre={}, post={})
            for method in methods:
                for session in ['pre','post']:
                    for cond in CS_cond.keys():
                        stc_file = stc_path + '/' + trial_type + '_' + session + \
                                '-' + fwd_params['spacing'] + '_' + CS_cond[cond] + '_' + method

                        stc[session].update({cond: mne.read_source_estimate(stc_file)})

                for key in [k for k in do_contrasts.keys() if do_contrasts[k]]:
                    if key == 'csXoddXsession':
                        stc_con = ( ( stc['post']['devCSp'] - stc['post']['stdCSp'] ) - \
                                ( stc['post']['devCSm'] - stc['post']['stdCSm'] )) - \
                                ( ( stc['pre']['devCSp'] - stc['pre']['stdCSp'] ) - \
                                ( stc['pre']['devCSm'] - stc['pre']['stdCSm'] ) )
                    elif key == 'csXsession_dev':
                        stc_con = ( stc['post']['devCSp'] - stc['post']['devCSm'] ) - \
                                ( stc['pre']['devCSp'] - stc['pre']['devCSm'] )
                    elif key == 'odd_pre':
                        stc_con = ( stc['pre']['devCSp'] - stc['pre']['stdCSp'] + \
                                    stc['pre']['devCSm'] - stc['pre']['stdCSm'] ) / 2.
                    elif key == 'oddXsession':
                        stc_con = ( stc['post']['devCSp'] - stc['post']['stdCSp'] + \
                                    stc['post']['devCSm'] - stc['post']['stdCSm'] ) - \
                                  ( stc['pre']['devCSp'] - stc['pre']['stdCSp'] + \
                                    stc['pre']['devCSm'] - stc['pre']['stdCSm'] )

                    stc_file = stc_path_out + '/' + trial_type + '-' + fwd_params['spacing'] + \
                            '_' + key + '_' + method
                    stc_con.save(stc_file)

if do_morph_contrasts_to_fsaverage:

    do_contrasts = {'odd_pre': True, 'oddXsession': True, 'csXsession_dev': True, 'csXoddXsession': True}
    methods = ['MNE','dSPM']

    # This seems very hacky, but will have to try to under
    # stand later...
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_to = 'fsaverage'

    for subject_from in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subject_from + '/sourceCon'
        morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to + '/sourceCon'
        mkdir_p(morph_stc_path)

        for trial_type in ['VS','FB']:
            print 20*'#'
            print '%s -> %s' % (subject_from, trial_type)
            print 20*'#'
            for method in methods:
                for key in [k for k in do_contrasts.keys() if do_contrasts[k]]:
                    stc_file = stc_path + '/' + trial_type + '-' + fwd_params['spacing'] + \
                            '_' + key + '_' + method
                    stc_from = mne.read_source_estimate(stc_file)
                    stc_to = mne.morph_data(subject_from, subject_to,
                            stc_from, grade=vertices_to, n_jobs=6)

                    stc_file = morph_stc_path + '/' + subject_from + '_' + trial_type + \
                            '-' + fwd_params['spacing'] + '_' + key + '_' + method
                    stc_to.save(stc_file)

if do_average_morph_maps:

    methods = ['MNE','dSPM']
    morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to + '/sourceCon'
    for trial_type in ['VS','FB']:
        for method in methods:
            for key in [k for k in do_contrasts.keys() if do_contrasts[k]]:
                stc_list = []
                for subject_from in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:
                    stc_file = morph_stc_path + '/' + subject_from + '_' + trial_type + \
                            '-' + fwd_params['spacing'] + '_' + key + '_' + method
                    stc_list.append(mne.read_source_estimate(stc_file))

                stc_ave = stc_list[0]
                for ii in range(1,len(stc_list)):
                    stc_ave = (float(ii)*stc_ave + stc_list[ii]) / (float(ii)+1.)

                stc_file = morph_stc_path + '/AVG_' + trial_type + \
                        '-' + fwd_params['spacing'] + '_' + key + '_' + method
                stc_ave.save(stc_file)

##########
# once-for-all run for FFA-localizer to source contrasts!
##########

if do_FFA_SEs:

    # do_epoching
    for subj in ad.analysis_dict.keys():

        raw_path = ad._scratch_folder + '/filtered/' + filter_params['input_files'] + '/' + filt_dir + '/' + subj
        eve_path = ad._scratch_folder + '/events.fif/' + subj + '/raw'

        epo_path = ad._scratch_folder + '/epochs/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        img_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj + '/img'

        session_name = 'FFA'
        fname = raw_path + '/' + sesname + '_filt.fif'
        raw = Raw(fname, preload=False)
        events = mne.read_events(eve_path + '/' + sesname + '-eve.fif')
        picks = pick_types(raw.info, meg=True, eeg=False, stim=True, eog=True, misc=True)
        FFA_eve = mne.merge_events(events, [100, 200], 100, replace_events=True)
        id_dict = dict(face=100, blur=150)

        print('Extracting %s (%s) epochs for %s' % (trial_type, session, subj))
        epochs = mne.Epochs(raw, FFA_eve, id_dict,
                tmin, tmax, picks=picks, verbose=False,
                baseline=baseline, reject=reject, preload=True,
                reject_tmin=rej_tmin, reject_tmax=rej_tmax) # Check rejection settings
        # Check the drop_log for these preload'ed epochs: does the drop
        # log indicate the dropped epochs, can they be un-dropped after the fact?
        # Do we in fact have to actively drop them, despite preloading?

        print('Resampling...')
        epochs.resample(rsl_fs, n_jobs=6, verbose=False) # Trust the defaults here

        epo_out = epo_path + '/' + session_name + '-epo.fif'
        epochs.save(epo_out)  # save epochs to disk

        print('Making evokeds...')
        evokeds = []
        for categ in evoked_categories.keys():
            if len(evoked_categories[categ]) == 2:
                evokeds.append(epochs[evoked_categories[categ][0]].average() - \
                        epochs[evoked_categories[categ][1]].average())
                evokeds[-1].comment = categ
            else:
                evokeds.append(epochs[evoked_categories[categ][0]].average())
                evokeds[-1].comment = categ

        cov_all = mne.compute_covariance(epochs, tmin=baseline[0], tmax=baseline[1]) # same covariance for all contrasts
        figs = mne.viz.plot_cov(cov_all, epochs.info, show=False)
        figs[0].savefig(img_path + '/' + trial_type + '_' + session + '_all_covMAT.png')
        figs[1].savefig(img_path + '/' + trial_type + '_' + session + '_all_covSVD.png')


if do_sensor_level_contrasts:
    # This makes no sense: I might as well use the evoked contrasts to image
    # The Leff-calculation can perhaps be deffered to the higher-level cont's?

    conds = ['stdA', 'stdB', 'devA', 'devB']

    for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

        evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj

        # Drop the FFA session for now, deal with it separately, also empty room
        session_names = [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                if ('FFA' not in x and 'empty' not in x)]

        if '1a' in session_names[0]:
            CS_cond = dict(stdCSp='stdA',stdCSm='stdB',devCSp='devA',devCSm='devB')
        elif '1b' in session_names[0]:
            CS_cond = dict(stdCSp='stdB',stdCSm='stdA',devCSp='devB',devCSm='devA')

        for trial_type in ['VS','FB']:
            evo = dict(pre={}, post={})
            Leff = dict()
            evokeds = []
            for session in ['pre','post']:
                evo_file = evo_path + '/' + trial_type + '_' + session + '-avg.fif'
                for cond in CS_cond.keys():
                    evo[session].update({cond: mne.read_evokeds(evo_file, condition=CS_cond[cond])})

            # csXoddXsession = deepcopy(evo['post']['devCSp']) # get info stuff
            csXoddXsession = ( ( evo['post']['devCSp'] - evo['post']['stdCSp'] ) - \
                    ( evo['post']['devCSm'] - evo['post']['stdCSm'] )) - \
                    ( ( evo['pre']['devCSp'] - evo['pre']['stdCSp'] ) - \
                    ( evo['pre']['devCSm'] - evo['pre']['stdCSm'] ) )
            L = 1./( 1./evo['post']['devCSp'].nave + 1./evo['post']['stdCSp'].nave + \
                    1./evo['post']['devCSm'].nave + 1./evo['post']['stdCSm'].nave + \
                    1./evo['pre']['devCSp'].nave + 1./evo['pre']['stdCSp'].nave + \
                    1./evo['pre']['devCSm'].nave + 1./evo['pre']['stdCSm'].nave )
            evokeds.append(csXoddXsession)
            evokeds[-1].comment = 'csXoddXsession'
            Leff.update({'csXoddXsession': L})

            # Angry face vs. neutral, pre-conditioning
            odd_pre = ( ( evo['pre']['devCSp'] - evo['pre']['stdCSp'] ) + \
                    ( evo['pre']['devCSm'] - evo['pre']['stdCSm'] ) )
            evokeds.append(odd_pre)
            evokeds[-1].comment = 'odd_pre'
            L = 1./( 1./evo['pre']['devCSp'].nave + 1./evo['pre']['stdCSp'].nave + \
                    1./evo['pre']['devCSm'].nave + 1./evo['pre']['stdCSm'].nave )
            Leff.update({'odd_pre': L})
            # Angry face vs. neutral, post-conditioning
            odd_post = ( ( evo['post']['devCSp'] - evo['post']['stdCSp'] ) + \
                    ( evo['post']['devCSm'] - evo['post']['stdCSm'] ) )
            evokeds.append(odd_post)
            evokeds[-1].comment = 'odd_post'
            L = 1./( 1./evo['post']['devCSp'].nave + 1./evo['post']['stdCSp'].nave + \
                    1./evo['post']['devCSm'].nave + 1./evo['post']['stdCSm'].nave )
            Leff.update({'odd_post': L})
            ##########################################

            csXodd_pre = ( ( evo['pre']['devCSp'] - evo['pre']['stdCSp'] ) - \
                    ( evo['pre']['devCSm'] - evo['pre']['stdCSm'] ) )
            evokeds.append(csXodd_pre)
            evokeds[-1].comment = 'csXodd_pre'
            L = 1./( 1./evo['pre']['devCSp'].nave + 1./evo['pre']['stdCSp'].nave + \
                    1./evo['pre']['devCSm'].nave + 1./evo['pre']['stdCSm'].nave )
            Leff.update({'csXodd_pre': L})

            # csXodd_post = deepcopy(evo['post']['devCSp'])
            csXodd_post = ( ( evo['post']['devCSp'] - evo['post']['stdCSp'] ) - \
                    ( evo['post']['devCSm'] - evo['post']['stdCSm'] ) )
            evokeds.append(csXodd_post)
            evokeds[-1].comment = 'csXodd_post'
            L = 1./( 1./evo['post']['devCSp'].nave + 1./evo['post']['stdCSp'].nave + \
                    1./evo['post']['devCSm'].nave + 1./evo['post']['stdCSm'].nave )
            Leff.update({'csXodd_post': L})

            cs_dev_pre = evo['pre']['devCSp'] -  evo['pre']['devCSm']
            evokeds.append(cs_dev_pre)
            evokeds[-1].comment = 'cs_dev_pre'
            L = 1./( 1./evo['pre']['devCSp'].nave + 1./evo['pre']['devCSm'].nave )
            Leff.update({'cs_dev_pre': L})

            con_file = evo_path + '/' + trial_type + '-contrasts.fif'
            mne.write_evokeds(con_file, evokeds)
            f = open(evo_path + '/' + trial_type + '-contrasts.Leff', 'w')
            json.dump(Leff, f)
            f.close()

if do_sensor_level_contrast_images_across:

    do_contrasts = {'odd_pre': True, 'odd_post': True, 'csXoddXsession': True,'csXodd_pre': True, 'csXodd_post': True}
    import matplotlib.pyplot as plt
    clim_con = dict(mag=[-50, 50], grad=[0, 15])
    topo_times = np.arange(-0.030, 0.250, 0.030)

    img_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/across/img'
    mkdir_p(img_path)
    group_con_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/across'

    for trial_type in ['VS','FB']:
        group_evokeds = []
        for contrast in [k for k in do_contrasts.keys() if do_contrasts[k]]:
        # for contrast in ['csXoddXsession','csXodd_pre','csXodd_post']:
            for ii,subj in enumerate([x for x in sorted(ad.analysis_dict.keys()) if 'T1' in ad.analysis_dict[x].keys()]):
                evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
                con_file = evo_path + '/' + trial_type + '-contrasts.fif'
                evo = mne.read_evokeds(con_file, condition=contrast)
                #f = open(evo_path + '/' + trial_type + '-contrasts.Leff', 'r')
                #Leff = json.load(f)
                #f.close()
                try:
                    evo_avg.data = ( float(ii)*evo_avg.data + evo.data ) / float(ii+1)
                except:
                    evo_avg = evo # first subject = 006_HEN when sorted!

            evo_avg.plot_image(clim=clim_con, show=False)
            plt.savefig(img_path + '/' + trial_type + '-' + contrast + '_time.png')
            plot_evoked_topomap(evo_avg,topo_times, show=False, vmin=[clim_con['grad'][0],clim_con['mag'][0]], vmax=[clim_con['grad'][1],clim_con['mag'][1]])
            plt.savefig(img_path + '/' + trial_type + '-' + contrast + '_topo.png')

            evo_avg.comment = contrast
            group_evokeds.append(evo_avg.copy())
            # group_evokeds[-1].comment = contrast
            evo_avg = []

        mne.write_evokeds(group_con_path + '/SENSCON-' + trial_type + '.fif', group_evokeds)

if do_sensor_level_contrast_to_sources:

    snr = 1.0
    lambda2 = 1.0 / snr ** 2
    methods = ['MNE','dSPM']
    ori_sel = 'normal'

    #contrast_prefix = 'csXodd_'
    contrast_prefix = 'odd_'

    do_contrasts = {'odd_pre': True, 'odd_post': True, 'csXoddXsession': True,'csXodd_pre': True, 'csXodd_post': True}

    for trial_type in ['VS','FB']:
        for session in ['pre','post']:
            contrast = contrast_prefix + session
            for ii,subj in enumerate([x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]):
                evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
                opr_path = ad._scratch_folder + '/operators/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
                con_file = evo_path + '/' + trial_type + '-contrasts.fif'
                cov_file = evo_path + '/' + trial_type + '_' + session + '-cov.fif'
                fwd_file = opr_path + '/' + trial_type + '_' + session + \
                        '-' + fwd_params['spacing'] + '-fwd.fif'
                stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj

                evo = mne.read_evokeds(con_file, condition=contrast)

# Don't do Leff correction
#                f = open(evo_path + '/' + trial_type + '-contrasts.Leff', 'r')
#                Leff = json.load(f)
#                f.close()
#                evo.nave = Leff[contrast] # Try forcing this!

                noise_cov = mne.read_cov(cov_file)
                # regularize noise covariance
                noise_cov = mne.cov.regularize(noise_cov, evo.info,
                        mag=0.05, grad=0.05, proj=True)

                fwd_opr = mne.read_forward_solution(fwd_file, surf_ori=True)

                inv_opr = mne.minimum_norm.make_inverse_operator(evo.info,
                        fwd_opr, noise_cov, loose=0.2, depth=0.8)

                for method in methods:
                    # if nave is set to Leff, this will apply it
                    stc = mne.minimum_norm.apply_inverse(evo, inv_opr,
                            lambda2, method, pick_ori=ori_sel)

                    # Save result in stc files
                    stc_file = stc_path + '/SENSCON_' + trial_type + \
                            '-' + fwd_params['spacing'] + '_' + contrast + '_' + method
                    stc.save(stc_file)

if do_sensor_level_contrast_to_sources_to_fsaverage:

    methods = ['MNE','dSPM']

    #contrast_prefix = '_csXodd_'
    contrast_prefix = '_odd_'

    # This seems very hacky, but will have to try to under
    # stand later...
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_to = 'fsaverage'

    morph_stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/morph-' + subject_to

    for trial_type in ['VS','FB']:
        for method in methods:
            stc_list = {'pre': [], 'post': [], 'diff': []}
            for subject_from in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:

                stc_path = ad._scratch_folder + '/estimates/' + filt_dir + '/' + filter_params['input_files'] + '/' + subject_from
                stc_file_pre = stc_path + '/SENSCON_' + trial_type + '-' + \
                        fwd_params['spacing'] + contrast_prefix + 'pre_' + method
                stc_file_post = stc_path + '/SENSCON_' + trial_type + '-' + \
                        fwd_params['spacing'] + contrast_prefix + 'post_' + method
                stc_pre = mne.read_source_estimate(stc_file_pre)
                stc_post = mne.read_source_estimate(stc_file_post)

                stc_pre_avg = mne.morph_data(subject_from, subject_to,
                        stc_pre, grade=vertices_to, n_jobs=6)
                stc_post_avg = mne.morph_data(subject_from, subject_to,
                        stc_post, grade=vertices_to, n_jobs=6)

                stc_diff = stc_post_avg - stc_pre_avg
                stc_list['pre'].append(stc_pre_avg)
                stc_list['post'].append(stc_post_avg)
                stc_list['diff'].append(stc_diff)

            for stc_key in stc_list.keys():
                stc_ave = stc_list[stc_key][0]
                for ii in range(1,len(stc_list[stc_key])):
                    stc_ave = (float(ii)*stc_ave + stc_list[stc_key][ii]) / (float(ii)+1.)

                stc_file = morph_stc_path + '/SENSCON_' + trial_type + \
                        '-' + fwd_params['spacing'] + contrast_prefix + \
                        stc_key + '_' + method
                stc_ave.save(stc_file)

if False:
    from mne.viz import plot_image_epochs

    import matplotlib.pyplot as plt
    subj = '007_SGF'

    evo_path = ad._scratch_folder + '/evoked/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
    session = 'pre'
    contrast = 'face'
    trial_type = 'FB'

    evo_in = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-ave.fif'

    evoked = read_evoked(evo_in)
    mag_picks = pick_types(evoked.info, meg='mag', eeg=False, ref_meg=False,
                   exclude='bads')
    vmin, vmax = evoked.data[mag_picks,:].min(), evoked.data[mag_picks].max()
    ax1 = plt.subplot2grid((3, 20), (0, 0), colspan=9, rowspan=2)
    im = plt.imshow(evoked.data[mag_picks,:],
                    extent=[1e3 * evoked.times[0], 1e3 * evoked.times[-1],
                            0, len(mag_picks)],
                    aspect='auto', origin='lower',
                    vmin=vmin, vmax=vmax)
    ax2 = plt.subplot2grid((3, 20), (2, 0), colspan=9, rowspan=1)
    #if colorbar:
    #    ax3 = plt.subplot2grid((3, 10), (0, 9), colspan=1, rowspan=3)
    ax1.set_title('Magnetometers')
    ax1.set_ylabel('Channel')
    ax1.axis('auto')
    ax1.axis('tight')
    ax1.axvline(0, color='m', linewidth=3, linestyle='--')
    #ax2.plot(1e3 * evoked.times, scalings[ch_type] * evoked.data[i])
    ax2.plot(1e3 * evoked.times, evoked.data[mag_picks,:].T)
    ax2.set_xlabel('Time (ms)')
    #ax2.set_ylabel(units[ch_type])
    #ax2.set_ylim([vmin, vmax])
    ax2.axvline(0, color='m', linewidth=3, linestyle='--')
    # if colorbar:
     #    plt.colorbar(im, cax=ax3)
      #   tight_layout(fig=this_fig)


    grad_picks = pick_types(evoked.info, meg='grad', eeg=False, ref_meg=False,
               exclude='bads')
    vmin, vmax = evoked.data[grad_picks,:].min(), evoked.data[grad_picks].max()
    ax1 = plt.subplot2grid((3, 20), (0, 10), colspan=9, rowspan=2)
    im = plt.imshow(evoked.data[grad_picks,:],
                    extent=[1e3 * evoked.times[0], 1e3 * evoked.times[-1],
                            0, len(grad_picks)],
                    aspect='auto', origin='lower',
                    vmin=vmin, vmax=vmax)
    ax2 = plt.subplot2grid((3, 20), (2, 10), colspan=9, rowspan=1)
    #if colorbar:
    #    ax3 = plt.subplot2grid((3, 10), (0, 9), colspan=1, rowspan=3)
    ax1.set_title('Gradiometers')
    ax1.set_ylabel('Channel')
    ax1.axis('auto')
    ax1.axis('tight')
    ax1.axvline(0, color='m', linewidth=3, linestyle='--')
    #ax2.plot(1e3 * evoked.times, scalings[ch_type] * evoked.data[i])
    ax2.plot(1e3 * evoked.times, evoked.data[grad_picks,:].T)
    ax2.set_xlabel('Time (ms)')
    #ax2.set_ylabel(units[ch_type])
    #ax2.set_ylim([vmin, vmax])
    ax2.axvline(0, color='m', linewidth=3, linestyle='--')
    plt.show()

###############################################################################
# View evoked response
# times = 1e3 * epochs.times  # time in miliseconds

# ch_max_name, latency = evoked.get_peak(mode='neg')

# import matplotlib.pyplot as plt
# evoked.plot()

# plt.xlim([times[0], times[-1]])
# plt.xlabel('time (ms)')
# plt.ylabel('Potential (uV)')
# plt.title('EEG evoked potential')

# plt.axvline(latency * 1e3, color='red',
#             label=ch_max_name, linewidth=2,
#             linestyle='--')
# plt.legend(loc='best')

# plt.show()

