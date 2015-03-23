# Since March 2015, ICA has been performed on the data, and we'll be using
# Savitzky-Golay for epochs filtering
#
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
CLOBBER=True
import mne
#try:
from mne.io import Raw
from mne import pick_types, compute_covariance, read_epochs
from mne import read_evokeds, write_evokeds
from mne import read_forward_solution, read_cov
from mne import read_source_estimate
from mne.forward import do_forward_solution
from mne.minimum_norm import (make_inverse_operator, write_inverse_operator,
        read_inverse_operator, apply_inverse)
from mne.viz import plot_cov
from mne.report import Report
import numpy as np
from operator import add # for stc reduction operation
#from viz_cjb import plot_evoked_topomap

import nibabel as nib

# get basic stuff like mkdir_p and some defaults
from VSC_utils import *

import matplotlib
matplotlib.use('agg') # force non-interactive plotting
#matplotlib.use('Qt4Agg') # force non-interactive plotting
import matplotlib.pyplot as plt

#import os, errno
import json, subprocess
from copy import deepcopy

machine_name = os.uname()[1].split('.')[0]

if 'isis' in machine_name:
    import sys
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
    import subprocess
    from access import Query
    from analysis_dict import Anadict

    db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
    ad=Anadict(db)
elif 'mba-cjb' in machine_name or 'hathor' in machine_name:
    class local_Anadict():
        def __init__(self):
            self._scratch_folder = '/Users/cjb/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scratch'

    class local_Query():
        def get_subjects(self):
            return ['030_WAH']
    
    db = local_Query()
    ad = local_Anadict()

fs_subjects_dir = ad._scratch_folder + '/fs_subjects_dir'
#from mne import read_evokeds
#except:
#    from mne.fiff import Raw, pick_types, read_evoked

do_evokeds = False
do_forward_solutions_evoked = False
do_inverse_operators_evoked = False

# localize the face vs blur (diff) condition
# also do just face to get a nice map
do_STC_FFA = False
plot_STC_FFA = True

# create an average brain from participants, not fsaverage!
do_make_average_subject = False
do_make_morph_maps_to_VSaverage = False
do_average_morphed_evokedSEs = False

do_morph_evokedSEs_to_fsaverage = False
do_grandaverage_CScontrasts = False

do_sourcelevel_rmanova_stclustering = False

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


###################################
epo_folder = ad._scratch_folder + '/epochs/ica/' + filter_params['input_files']
evo_folder = ad._scratch_folder + '/evoked/ica/' + filter_params['input_files']
opr_folder = ad._scratch_folder + '/operators/ica/' + filter_params['input_files']
stc_folder = ad._scratch_folder + '/estimates/ica/' + filter_params['input_files']
tra_folder = ad._scratch_folder + '/trans'
###################################

if do_evokeds: # do a couple of "main effects"

    rep_folder = evo_folder + '/report'
    mkdir_p(rep_folder)

    topo_times = np.concatenate((np.arange(0.05, 0.110,0.010), 
        np.arange(0.12, 0.210,0.020)))

    #for subj in ['007_SGF']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        report = Report(info_fname=None, subjects_dir=None, subject=subj,
                        title='Evoked responses',
                        verbose=None)

        epo_path = epo_folder + '/' + subj
        evo_path = evo_folder + '/' + subj
        rep_file = rep_folder + '/' + subj + '.html'
        mkdir_p(evo_path)

        session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        for trial_type in ['VS','FB','FFA']:
            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                evokeds = []
                # evoked_categories loaded from VSC_utils.py
                evocat_sorted = evoked_categories[trial_type].keys()[:]
                evocat_sorted.sort()
                for categ in evocat_sorted:
                    if len(evoked_categories[trial_type][categ]) == 2:
                        # remember to equalize counts! dropping info on which
                        # were dropped...
                        epo,_ = epochs.equalize_event_counts(
                            event_ids=evoked_categories[trial_type][categ],
                            method='mintime', copy=True)

                        evokeds.append(epo[evoked_categories[trial_type][categ][0]].average() - \
                                epo[evoked_categories[trial_type][categ][1]].average())
                        evokeds[-1].comment = categ
                    else:
                        evokeds.append(epochs[evoked_categories[trial_type][categ][0]].average())
                        evokeds[-1].comment = categ

                #print trial_type, session, ': Estimating (optimal) covariance matrix'
                #noise_cov = compute_covariance(epochs, 
                #        method='auto',return_estimators=False,
                #        tmin=baseline[0], tmax=baseline[1]) # take the BEST estimate!
                noise_cov = compute_covariance(epochs, method='shrunk',
                        return_estimators=True,
                        tmin=baseline[0], tmax=baseline[1])[0] #take best
                figs = plot_cov(noise_cov, epochs.info, show=False)
                captions = [trial_type+session+':COV',
                        trial_type+session+':SVD']
                report.add_figs_to_section(figs, captions=captions,
                        section='Covar', scale=None, image_format='png')
                #figs[0].savefig(img_path + '/' + trial_type + '_' + session + '_all_covMAT.png')
                #figs[1].savefig(img_path + '/' + trial_type + '_' + session + '_all_covSVD.png')
                plt.close(figs[0])
                plt.close(figs[1])

                for e in evokeds:
                    #if savgol_hf_evo is not None:
                    #    e.savgol_filter(savgol_hf_evo)
                    figs = []
                    figs.append(e.plot_white(noise_cov, show=False))
                    figs.append(e.plot_topomap(times=topo_times,ch_type='mag'))
                    figs.append(e.plot_topomap(times=topo_times,ch_type='grad'))
                    captions = [e.comment+'-butterfly',e.comment+'-MAGtopo',e.comment+'-GRADtopo'] 
                    report.add_figs_to_section(figs, captions=captions,
                        section=trial_type + session,
                        scale=None, image_format='png')
                    for fig in figs:
                        plt.close(fig)

                evo_out= evo_path + '/' + trial_type + session + '-avg.fif'
                write_evokeds(evo_out, evokeds)  # save evoked data to disk

                cov_out = evo_path + '/' + trial_type + session + '-cov.fif'
                print cov_out
                noise_cov.save(cov_out)  # save covariance data to disk

        report.save(fname=rep_file, open_browser=False, overwrite=CLOBBER)

if do_forward_solutions_evoked:
    # modified to use the mne-python wrapper instead of calling command line
    # directly. See pipeline.py for the bash-way, which might be interesting
    # for an OGE-aware implementation?

    # check that 'T1' is attached to subject first, assume then MR preproc OK
    #for subj in [x for x in ad.analysis_dict.keys() if 'T1' in ad.analysis_dict[x].keys()]:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        trans_fif = tra_folder + '/' + subj + '-trans.fif'
        evo_path = evo_folder + '/' + subj
        fwd_path = opr_folder + '/' + subj
        mkdir_p(fwd_path)

        # HAve to assume all require their own because of ICA. Is this so?
        session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        for trial_type in ['VS','FB','FFA']:
            for session in session_nos[trial_type]:
                evo_file = evo_path + '/' + trial_type + session + '-avg.fif'

                fwd_out = fwd_path + '/' + trial_type + session + \
                                    '-' + fwd_params['spacing'] + '-fwd.fif'
                if os.path.exists(fwd_out):
                    continue

                do_forward_solution(subj, evo_file,
                        fname=fwd_out, #destination name
                        src=None, # Because spacing-param gives standard name!
                        spacing=fwd_params['spacing'],
                        mindist=fwd_params['mindist'],
                        bem=subj + fwd_params['bem'],
                        trans=None, mri=trans_fif,
                        eeg=False, meg=True, fixed=False, grad=False,
                        mricoord=False, overwrite=CLOBBER, subjects_dir=None,
                        verbose=None)


if do_inverse_operators_evoked:

    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        evo_path = evo_folder + '/' + subj
        opr_path = opr_folder + '/' + subj

        session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        for trial_type in ['VS','FB','FFA']:
            for session in session_nos[trial_type]:

                evo_file = evo_path + '/' + trial_type + session + '-avg.fif'
                cov_file = evo_path + '/' + trial_type + session + '-cov.fif'
                fwd_file = opr_path + '/' + trial_type + session + \
                        '-' + fwd_params['spacing'] + '-fwd.fif'
                inv_file = opr_path + '/' + trial_type + session + \
                        '-' + fwd_params['spacing'] + '-inv.fif'
                if file_exists(inv_file) and not CLOBBER:
                    continue

                # Load data
                evokeds = read_evokeds(evo_file)
                fwd_opr = read_forward_solution(fwd_file, surf_ori=True)
                noise_cov = read_cov(cov_file) # assume regularized

                # Here assuming that the nave of evoked doesn't influence
                # the inverse operator. May want to check this via mailing
                # list? A quick perusal of the code doesn't make anything
                # stand out: nave isn't used. Plus, if en empty room
                # noise cov were used here, it would by definition be 
                # independent of nave, so scaling it (by nave) wouldn't
                # make sense anyway.
                # Thus: just using the info from evokeds[0]
                if isinstance(evokeds, mne.evoked.Evoked):
                    info = evokeds.info
                elif isinstance(evokeds, list):
                    info = evokeds[0].info

                inv_opr = make_inverse_operator(info,
                        fwd_opr, noise_cov,
                        limit_depth_chs=inv_params['limit_depth_params'],
                        loose=inv_params['loose'],
                        depth=inv_params['depth'],
                        fixed=inv_params['fixed'])

                write_inverse_operator(inv_file, inv_opr)

if do_STC_FFA:
    # looking at the evokeds, it seems there's plenty to
    # see even efter 200, probably even longer.
    time_range = (-0.100, 0.300)
    methods = ['MNE','dSPM']
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)

    trial_type = 'FFA'
    session = ''
    do_evoked_contrasts = {'face': True, 'diff': True}
    SNRs = {'face': 3., 'diff': 3.}

    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        evo_path = evo_folder + '/' + subj
        opr_path = opr_folder + '/' + subj
        stc_path = stc_folder + '/' + subj

        evo_file = evo_path + '/' + trial_type + session + '-avg.fif'
        inv_file = opr_path + '/' + trial_type + session + \
                '-' + fwd_params['spacing'] + '-inv.fif'

        print 'Loading inverse operator...'
        inv_opr = read_inverse_operator(inv_file, verbose=False)

        for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
            # Load data
            evoked = read_evokeds(evo_file, condition=cond, verbose=False)
        
            lambda2 = 1. / SNRs[cond] ** 2.
            for method in methods:
                # Save result in stc files
                stc_file = stc_path + '/' + trial_type + session + \
                        '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                if file_exists(stc_file) and not CLOBBER:
                    continue

                #print 'Applying inverse with method:', method
                stc = apply_inverse(evoked, inv_opr,
                        lambda2, method, pick_ori=ori_sel, verbose=False)
                stc.crop(tmin=time_range[0], tmax=time_range[1]) # CROP

                stc.save(stc_file, verbose=False)

if plot_STC_FFA:

    #methods = ['MNE','dSPM']
    # Don't run MNE, just use dSPM for visualization
    methods = ['dSPM',]
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)
    trial_type = 'FFA'
    session = ''
    do_evoked_contrasts = {'face': True, 'diff': True}

    # NB! This needs to run headless !
    # Note the need to specify the server number: default is 99 :)
    # xvfb-run -n 98 --server-args="-screen 0 1024x768x24" python evoked_pipeline.py
#    from mayavi import mlab
#    mlab.options.offscreen = True

    tmp_folder = ad._scratch_folder + '/tmp/'
    tmp_file_suffix = '.brain-tf_%02d.png'
    brain_times = np.array([60., 80., 100., 120., 140., 160., 180.])
    use_abs_idx = False # Just use increments
    montage = [['lat', 'med'],['cau','ven']]
    stcran = dict(MNE={'max': 0.9, 'min': 0.1},
            dSPM={'max': 0.8, 'min': 0.2})
    #views=dict(rh=[((-140,124),'med'), ((-33,92),'lat')],
    #        lh=[((-37,120),'med'),((-147,90),'lat')])
    #tmpfile=dict(lat=ad._scratch_folder + '/tmp/brain-lat.tiff',
    #        med=ad._scratch_folder + '/tmp/brain-med.tiff',
    #        both=ad._scratch_folder + '/tmp/brain.tiff')

    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        stc_path = stc_folder + '/' + subj
        rep_folder = stc_path + '/report' #here under subj to reduce clutter
        mkdir_p(rep_folder)
        rep_file = rep_folder + '/' + subj + '-FFA.html'

        #  cannot be loaded/appended :(
        report = Report(info_fname=None, 
                subjects_dir=fs_subjects_dir, subject=subj,
                title='FFA estimates', verbose=None)

        for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
            # Load data
            for method in methods:
                # Save result in stc files
                stc_file = stc_path + '/' + trial_type + session + \
                        '-' + fwd_params['spacing'] + '_' + cond + '_' + method

                stc = read_source_estimate(stc_file)

                fmax = stcran[method]['max']*np.ravel(stc.data).max()
                fmin = stcran[method]['min']*fmax
                fmid = (fmax - fmin) / 2.

                for hemi in ['lh','rh']:
                    print 'Hemi :', hemi
                    #fig = mlab.figure(size=(400, 400))
                    brain = stc.plot(surface='inflated', hemi=hemi,
                            subject=subj, alpha = 0.9,
                            subjects_dir=fs_subjects_dir,
                            fmin=fmin, fmid=fmid, fmax=fmax)
                                
                    #aparc_file= os.path.join(fs_subjects_dir,
                    #      subj, "label",
                    #      hemi + ".aparc.a2009s.annot")
                    #labels, ctab, names = nib.freesurfer.read_annot(aparc_file)
                    #FFA_idx = names.index('G_oc-temp_lat-fusifor')
                    #roi_data = np.zeros(labels.shape)
                    #roi_data[labels==FFA_idx] = 1.
                    #brain.add_data(roi_data)
                    brain.add_label("V1", color='springgreen',
                            borders=False, alpha=0.2)
                    brain.add_label("V1", color='springgreen',
                            borders=True, alpha=1.)
                    brain.add_label("fusiform", color='aquamarine',
                            borders=False, alpha=0.2)
                    brain.add_label("fusiform", color='aquamarine',
                            borders=True, alpha=1.)

                    time_idx = [brain.index_for_time(t) for t in brain_times]

                    tmp_pattern = tmp_folder + hemi + tmp_file_suffix
                    brain.save_image_sequence(time_idx, tmp_pattern,
                            use_abs_idx=use_abs_idx, montage=montage)

                for ii,tt in enumerate(brain_times):
                    cmd = 'montage -geometry +4+4 '
                    for hemi in ['lh', 'rh']:
                        cmd += tmp_folder + hemi + tmp_file_suffix % (ii) + ' '
                    tmpname = tmp_folder + 'both' + tmp_file_suffix
                    cmd +=  tmpname
                    caption = method + ' @ %.0fms' % (tt)
                    tmpname = tmp_pattern % (ii)
                    report.add_images_to_section(tmpname, captions=caption,
                            section=cond, scale=0.5)
        report.save(fname=rep_file, open_browser=False, overwrite=CLOBBER)



if do_make_average_subject:
    subj_list = db.get_subjects() #only included subjects
    subjects_str = ' '.join([s[1:] for s in subj_list]) #strip the first zero :(
    fs_cmd = 'make_average_subject --out VSaverage --subjects ' + subjects_str
    #print fs_cmd
    proc = subprocess.Popen([fs_cmd], shell=True)
    proc.communicate()
    print 'make_average_subject returned code', proc.returncode

if do_make_morph_maps_to_VSaverage:
    for subj in db.get_subjects():
        subj = subj[1:] #strip the first zero :(
        print 'Morphing', subj
        mne.surface.read_morph_map(subj, 'VSaverage')

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

    evoked = mne.read_evokeds(evo_in)
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

