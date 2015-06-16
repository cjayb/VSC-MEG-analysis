# Since March 2015, ICA has been performed on the data, and we'll be using
# Savitzky-Golay for epochs filtering
#
# After initial round of analysis in Matlab (Fieldtrip, J-R King), most sensor-
# level results were replicated in python.
#   * there were differences in the csXoddXsession interaction, even though
#     the csXodd's for pre and post seemed similar...?
#   * JRK used "robust averaging", whereas in python we reject
#   * ICA used here for removing eye and heart activity (as best we can)

#
CLOBBER=False

do_evokeds = False
do_forward_solutions_evoked = False
do_inverse_operators_evoked = False

# localize the face vs blur (diff) condition
# also do just face to get a nice map
do_STC_FFA = False
plot_STC_FFA = False
do_STC_FFA_groupavg = False

# Decoding
do_GAT_FFA = False
do_GAT_FFA_scaledLR = False
do_GAT_VS_N2pc = False
do_GAT_VS_anyTRG = False

do_GAT_FB_anyTRG = False
do_GAT_FB_AtoB = False

do_GAT_FB_identityCS = True

# Now all group stats done at once
do_GAT_groupstat = True

# Try to generate some N2pc plots
do_N2pc_evokeds = False
do_STC_N2pc = False
plot_STC_N2pc = False
do_STC_N2pc_groupavg = False

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

###############################

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
from mne.decoding import GeneralizationAcrossTime
from mne.stats import spatio_temporal_cluster_1samp_test
import numpy as np
from operator import add # for stc reduction operation
#from viz_cjb import plot_evoked_topomap

import nibabel as nib

# get basic stuff like mkdir_p and some defaults
from VSC_utils import *

import pickle
import copy

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
    #sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/stormdb')
    sys.path.append('/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/scripts/VSC-MEG-analysis')
    import subprocess
    from stormdb.access import Query
    from analysis_dict import Anadict

    db=Query('MINDLAB2013_01-MEG-AttentionEmotionVisualTracking')
    ad=Anadict(db)
#elif 'mba-cjb' in machine_name or 'hathor' in machine_name:
else:
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



###################################
epo_folder = ad._scratch_folder + '/epochs/ica/' + filter_params['input_files']
evo_folder = ad._scratch_folder + '/evoked/ica/' + filter_params['input_files']
opr_folder = ad._scratch_folder + '/operators/ica/' + filter_params['input_files']
stc_folder = ad._scratch_folder + '/estimates/ica/' + filter_params['input_files']
dec_folder = ad._scratch_folder + '/decoding/ica/' + filter_params['input_files']
rep_folder = ad._scratch_folder + '/reports/ica/' + filter_params['input_files']
tra_folder = ad._scratch_folder + '/trans'
###################################

if do_evokeds: # do a couple of "main effects"

    evo_rep_folder = evo_folder + '/report'
    mkdir_p(evo_rep_folder)

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
        rep_file = evo_rep_folder + '/' + subj + '.html'
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

        report.save(fname=rep_file, open_browser=False, overwrite=True)

if do_N2pc_evokeds: # 

    N2pc_rep_folder = rep_folder + '/N2pc'
    mkdir_p(N2pc_rep_folder)

    topo_times = np.arange(0.12, 0.37,0.020)

    #for subj in ['007_SGF']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        report = Report(info_fname=None, subjects_dir=None, subject=subj,
                        title='N2pc responses',
                        verbose=None)

        epo_path = epo_folder + '/' + subj
        evo_path = evo_folder + '/' + subj
        rep_file = N2pc_rep_folder + '/' + subj + '.html'
        mkdir_p(evo_path)

        session_nos = dict(VS=['1','2'])
        for trial_type in ['VS',]:
            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                evokeds = []
                # evoked_categories loaded from VSC_utils.py
                evocat_sorted = evoked_categories['N2pc'].keys()[:]
                evocat_sorted.sort()
                for categ in evocat_sorted:
                    if len(evoked_categories['N2pc'][categ]) == 2:
                        # remember to equalize counts! dropping info on which
                        # were dropped...
                        epo,_ = epochs.equalize_event_counts(
                            event_ids=evoked_categories['N2pc'][categ],
                            method='mintime', copy=True)

                        evokeds.append(epo[evoked_categories['N2pc'][categ][0]].average() - \
                                epo[evoked_categories['N2pc'][categ][1]].average())
                        evokeds[-1].comment = categ
                    else:
                        evokeds.append(epochs[evoked_categories['N2pc'][categ][0]].average())
                        evokeds[-1].comment = categ

                for e in evokeds:
                    #if savgol_hf_evo is not None:
                    #    e.savgol_filter(savgol_hf_evo)
                    figs = []
                    noise_cov = mne.read_cov(evo_path + '/VS' + session + '-cov.fif')
                    figs.append(e.plot_white(noise_cov, show=False))
                    figs.append(e.plot_topomap(times=topo_times,ch_type='mag'))
                    figs.append(e.plot_topomap(times=topo_times,ch_type='grad'))
                    captions = [e.comment+'-butterfly',e.comment+'-MAGtopo',e.comment+'-GRADtopo'] 
                    report.add_figs_to_section(figs, captions=captions,
                        section=trial_type + session,
                        scale=None, image_format='png')
                    for fig in figs:
                        plt.close(fig)

                evo_out= evo_path + '/N2pc' + session + '-avg.fif'
                write_evokeds(evo_out, evokeds)  # save evoked data to disk

        report.save(fname=rep_file, open_browser=False, overwrite=True)

if do_GAT_FFA: # Generalization across time

    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    rep_file = gat_rep_folder + '/' + 'GAT_FFA.html'

    report = Report(info_fname=None, subjects_dir=None, subject=None,
                    title='Generalization Across Time (FFA)',
                    verbose=None)


    #for subj in ['030_WAH']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        epo_path = epo_folder + '/' + subj
        gat_path = dec_folder + '/' + subj
        mkdir_p(gat_path)

        session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        #for trial_type in ['VS','FB','FFA']:
        for trial_type in ['FFA']:
            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)

                # equalize event counts when using SVM
                # do this by passing the diff key
                epochs.equalize_event_counts(
                    event_ids=evoked_categories[trial_type]['diff'],
                    method='mintime', copy=False)

                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                # Need to redifine the events to only include 2 classes!
                # Dirty hack: modify the events in the epochs object directly!
                eve=epochs.events
                eve[eve[:,2]==200,2] = 100
                epochs.event_id = {u'face': 100, u'blur': 150}

                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=2)
                # If (clf is) None the classifier will be a standard pipeline including 
                # StandardScaler and a linear SVM with default parameters.

                figs = []
                # fit and score
                gat.fit(epochs)
                gat.score(epochs)
                figs.append(gat.plot(vmin=0.2, vmax=0.8,
                         title="Generalization Across Time (faces vs. blurred)"))
                figs.append(gat.plot_diagonal())  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/FFA-GAT.pickle', 'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT','Class']

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)
                    

if do_GAT_FFA_scaledLR: # Generalization across time with scaled Log Reg


    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    rep_file = gat_rep_folder + '/' + 'GAT_FFA_scaledLR.html'

    report = Report(info_fname=None, subjects_dir=None, subject=None,
                    title='Generalization Across Time (FFA), ' + \
                            'scaled Logistic Regression model',
                    verbose=None)


    #for subj in ['030_WAH']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        epo_path = epo_folder + '/' + subj
        gat_path = dec_folder + '/' + subj
        mkdir_p(gat_path)

        session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        #for trial_type in ['VS','FB','FFA']:
        for trial_type in ['FFA']:
            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)

                # equalize event counts when using SVM
                # do this by passing the diff key
                epochs.equalize_event_counts(
                    event_ids=evoked_categories[trial_type]['diff'],
                    method='mintime', copy=False)

                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                triggers = epochs.events[:,2]
                # for two classes
                # The following three are identical:
                y = np.in1d(triggers, (100,200)).astype(int) # face is one, blur is zero
                #y = 1. * epochs.events[:,2] != 150 # face is True, blur is False
                #y = np.in1d(triggers, (100,200)).astype(int) + 1# face is two, blur is one

                # what about 3 classes (faceA, faceB & blur)?
                #y = triggers
                # Works like fuckall, not getting the 3-label classification...?

                scaler = StandardScaler()
                clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                    # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                pipeline = Pipeline([('scaler', scaler), ('clf', clf)])

                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer)

                figs = []
                # fit and score
                gat.fit(epochs, y)
                gat.score(epochs)
                figs.append(gat.plot(vmin=0.1, vmax=0.9,
                         title="Generalization Across Time (faces vs. blurred)"))
                figs.append(gat.plot_diagonal(chance=0.5))  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/FFA-GAT-scaledLR.pickle', 'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT','Class']

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)


if do_GAT_VS_N2pc: # Generalization across time for visual search, N2pc

    tmin, tmax = -0.1, 0.35

    # evoked_categories loaded from VSC_utils.py
    # we'll use the 'diff' defined therein
    trgs_RvsL = evoked_categories['N2pc']['diff']
    # trgs_RvsL[0] = ['A1',...]
    # trgs_RvsL[1] = ['A4',...]
    #anytarget = tuple([element for lst in evList for element in lst])

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    rep_file = gat_rep_folder + '/' + 'GAT_VS1_N2pc.html'

    report = Report(info_fname=None, subjects_dir=None, subject=None,
                    title='Generalization Across Time, visual search, N2pc',
                    verbose=None)


    #for subj in ['030_WAH']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        epo_path = epo_folder + '/' + subj
        gat_path = dec_folder + '/' + subj
        mkdir_p(gat_path)

        #session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        # Only do first session for now...
        session_nos = dict(VS=['1',], FB=['1',], FFA=['',])
        #for trial_type in ['VS','FB','FFA']:
        for trial_type in ['VS']:

            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                epochs.drop_channels(['EOG001','EOG003'])

                # equalize event counts when using SVM
                # do this by passing the diff key
#                epochs.equalize_event_counts(
#                    event_ids=trgs_LvsR, # targets left vs. right
#                    method='mintime', copy=False)
                # Don't bother equalising, as we are doing LR by default
                # However, we will have to drop some epochs corresponding
                # to the trials we don't want to model!

                evoked_categories[trial_type]['stdA']
                triggers = epochs.events[:,2]
                events_to_drop = np.in1d(triggers, (
                    epochs.event_id[evoked_categories[trial_type]['stdA'][0][0]], 
                    epochs.event_id[evoked_categories[trial_type]['stdB'][0][0]],
                    epochs.event_id[evoked_categories[trial_type]['devA'][0][0]], 
                    epochs.event_id[evoked_categories[trial_type]['devB'][0][0]])
                        )
                #print 'Epochs before drop:', epochs
                epochs.drop_epochs(events_to_drop, reason='only deviants')
                #print 'Epochs after drop:', epochs

                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                triggers = epochs.events[:,2]
                # right hemifield target is one, right is zero
                y = np.in1d(triggers, 
                    tuple(epochs.event_id[t] for t in trgs_RvsL[0])).astype(int)


                #### Use scaled Logistic Regression per default
                scaler = StandardScaler()
                clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                    # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                pipeline = Pipeline([('scaler', scaler), ('clf', clf)])
                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer)

                figs = []
                # fit and score
                gat.fit(epochs, y)
                gat.score(epochs)
                figs.append(gat.plot(vmin=0.2, vmax=0.8,
                         title="Generalization Across Time (N2pc)"))
                figs.append(gat.plot_diagonal())  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/N2pc%s-GAT.pickle' % session, 
                        'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT-ses%s' % session,'Class-ses%s' % session]

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)
                    

if do_GAT_VS_anyTRG: # Generalization across time for visual search, any target

    tmin, tmax = -0.1, 0.35

    # evoked_categories loaded from VSC_utils.py
    # we'll use the 'diff' defined therein
    evList = evoked_categories['N2pc']['diff']
    anytarget = tuple([element for lst in evList for element in lst])

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    rep_file = gat_rep_folder + '/' + 'GAT_VS1_anyTRG.html'

    report = Report(info_fname=None, subjects_dir=None, subject=None,
                    title='Generalization Across Time, visual search, any target',
                    verbose=None)


    #for subj in ['030_WAH']:
    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        epo_path = epo_folder + '/' + subj
        gat_path = dec_folder + '/' + subj
        mkdir_p(gat_path)

        #session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
        # Only do first session for now...
        session_nos = dict(VS=['1',], FB=['1',], FFA=['',])
        #for trial_type in ['VS','FB','FFA']:
        for trial_type in ['VS']:

            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                epochs.drop_channels(['EOG001','EOG003'])

                # equalize event counts when using SVM
                # do this by passing the diff key
#                epochs.equalize_event_counts(
#                    event_ids=trgs_LvsR, # targets left vs. right
#                    method='mintime', copy=False)
                # Don't bother equalising, as we are doing LR by default
                # However, we will have to drop some epochs corresponding
                # to the trials we don't want to model!

                triggers = epochs.events[:,2]
                events_to_drop = np.in1d(triggers, (
                    epochs.event_id[evoked_categories[trial_type]['devA'][0][0]], 
                    epochs.event_id[evoked_categories[trial_type]['devB'][0][0]])
                        )
                #print 'Epochs before drop:', epochs
                epochs.drop_epochs(events_to_drop, reason='duplicate deviants')
                #print 'Epochs after drop:', epochs

                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                triggers = epochs.events[:,2]
                # ANY target is one, NO TRG is   zero
                y = np.in1d(triggers, 
                    tuple(epochs.event_id[t] for t in anytarget)).astype(int)


                #### Use scaled Logistic Regression per default
                scaler = StandardScaler()
                clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                    # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                pipeline = Pipeline([('scaler', scaler), ('clf', clf)])
                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer)

                # Define decoder. The decision_function is employed to use AUC for scoring
                #gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=4)
                # If (clf is) None the classifier will be a standard pipeline including 
                # StandardScaler and a linear SVM with default parameters.
                ###########

                figs = []
                # fit and score
                gat.fit(epochs, y)
                gat.score(epochs)
                figs.append(gat.plot(vmin=0.2, vmax=0.8,
                         title="Generalization Across Time (any TarGeT)"))
                figs.append(gat.plot_diagonal(chance=0.5))  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/VS%s-anyTRG-GAT.pickle' % session, 
                                        'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT-ses%s' % session,'Class-ses%s' % session]

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)
                    

if do_GAT_FB_anyTRG: # Generalization across time for feedback, any target

    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    #session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
    # Only do first session for now...
    session_nos = dict(VS=['1',], FB=['1','2'], FFA=['',])
    #for trial_type in ['VS','FB','FFA']:
    anytarget_dict = dict(FB=('devA','devB'))
    for trial_type in ['FB']:
        rep_file = gat_rep_folder + '/' + 'GAT_%s_anyTRG.html' % trial_type

        report = Report(info_fname=None, subjects_dir=None, subject=None,
                        title='Generalization Across Time, %s, any target' % trial_type,
                        verbose=None)

        anytarget = anytarget_dict[trial_type]

        #for subj in ['030_WAH']:
        for subj in db.get_subjects():
            if len(subj) == 8:
                subj = subj[1:]

            epo_path = epo_folder + '/' + subj
            gat_path = dec_folder + '/' + subj
            mkdir_p(gat_path)

            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                epochs.drop_channels(['EOG001','EOG003'])


                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                triggers = epochs.events[:,2]
                # ANY target is one, NO TRG is   zero
                y = np.in1d(triggers, 
                    tuple(epochs.event_id[t] for t in anytarget)).astype(int)


                #### Use scaled Logistic Regression per default
                scaler = StandardScaler()
                clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                    # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                pipeline = Pipeline([('scaler', scaler), ('clf', clf)])
                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer)

                # Define decoder. The decision_function is employed to use AUC for scoring
                #gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=4)
                # If (clf is) None the classifier will be a standard pipeline including 
                # StandardScaler and a linear SVM with default parameters.
                ###########

                figs = []
                # fit and score
                gat.fit(epochs, y)
                gat.score(epochs)
                figs.append(gat.plot(vmin=0.2, vmax=0.8,
                         title="Generalization Across Time (any target)"))
                figs.append(gat.plot_diagonal(chance=0.5))  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/%s%s-anyTRG-GAT.pickle' % (trial_type, session), 
                                        'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT-ses%s' % session,'Class-ses%s' % session]

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)
                    
if do_GAT_FB_identityCS: # Generalization across time for feedback
    # Try to classify CS+ from CS- in std and dev condition separately

    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    session_nos = dict(VS=['1',], FB=['1','2'], FFA=['',])

    for trial_type in ['FB']:
        rep_file = gat_rep_folder + '/' + 'GAT_%s_identityCS.html' % trial_type

        report = Report(info_fname=None, subjects_dir=None, subject=None,
                        title='Generalization Across Time, %s, CS+/-' % trial_type,
                        verbose=None)


        #for subj in ['030_WAH']:
        for subj in db.get_subjects():
            if len(subj) == 8:
                subj = subj[1:]

            epo_path = epo_folder + '/' + subj
            gat_path = dec_folder + '/' + subj
            mkdir_p(gat_path)

            session_names = \
                [x for x in ad.analysis_dict[subj][filter_params['input_files']].keys()
                if ('FFA' not in x and 'empty' not in x)]

            if '1a' in session_names[0]:
                CScode = {'CS+': 'A', 'CS-': 'B'}
            elif '1b' in session_names[0]:
                CScode = {'CS+': 'B', 'CS-': 'A'}

            for facetype, dropme in zip(['std','dev'],['dev','std']):

                for session in session_nos[trial_type]:
                    fname = epo_path + '/' + trial_type + session + '-epo.fif'

                    epochs = read_epochs(fname)
                    epochs.drop_channels(['EOG001','EOG003'])

                    triggers = epochs.events[:,2]
                    events_to_drop = np.in1d(triggers, (
                        epochs.event_id[evoked_categories[trial_type][dropme+'A'][0][0]], 
                        epochs.event_id[evoked_categories[trial_type][dropme+'B'][0][0]])
                            )
                    #print 'Epochs before drop:', epochs
                    epochs.drop_epochs(events_to_drop, reason='drop the other types')


                    if epoch_params['savgol_hf'] is not None:
                        epochs.savgol_filter(epoch_params['savgol_hf'])

                    epochs.crop(tmin=tmin, tmax=tmax)

                    triggers = epochs.events[:,2]
                    # CS+  is one, CS- is   zero
                    y = np.in1d(triggers, 
                            (epochs.event_id[facetype + CScode['CS+']], )).astype(int) 


                    #### Use scaled Logistic Regression per default
                    scaler = StandardScaler()
                    clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                        # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                    pipeline = Pipeline([('scaler', scaler), ('clf', clf)])
                    # Define decoder. The decision_function is employed to use AUC for scoring
                    gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer)

                    # Define decoder. The decision_function is employed to use AUC for scoring
                    #gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=4)
                    # If (clf is) None the classifier will be a standard pipeline including 
                    # StandardScaler and a linear SVM with default parameters.
                    ###########

                    figs = []
                    # fit and score
                    gat.fit(epochs, y)
                    gat.score(epochs)
                    figs.append(gat.plot(vmin=0.2, vmax=0.8,
                             title="GAT, %s, CS+ vs CS-" % (facetype)))
                    figs.append(gat.plot_diagonal(chance=0.5))  # plot decoding across time (correspond to GAT diagonal)

                    with open(gat_path + '/%s%s-identityCS-%s-GAT.pickle' % \
                            (trial_type, session, facetype), 'wb') as f:
                        pickle.dump(gat, f, protocol=2) # use optimised binary format

                    captions = [subj+'-GAT',subj+'-Diag'] 
                    sections = ['ses%s-%s' % (session, facetype),
                            'ses%s-%s' % (session, facetype)]

                    print 'Generating plots for', subj
                    for ff,fig in enumerate(figs):
                        report.add_figs_to_section(fig, captions=captions[ff],
                            section=sections[ff],
                            scale=None, image_format='png')
                        plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)

if do_GAT_FB_AtoB: # Generalization across time for feedback, any target

    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    #session_nos = dict(VS=['1','2'], FB=['1','2'], FFA=['',])
    # Only do first session for now...
    session_nos = dict(VS=['1',], FB=['1','2'], FFA=['',])
    #for trial_type in ['VS','FB','FFA']:
    for trial_type in ['FB']:
        rep_file = gat_rep_folder + '/' + 'GAT_%s_AtoB.html' % trial_type

        report = Report(info_fname=None, subjects_dir=None, subject=None,
                        title='Generalization Across Time, %s, face A to B' % trial_type,
                        verbose=None)

        #for subj in ['030_WAH']:
        for subj in db.get_subjects():
            if len(subj) == 8:
                subj = subj[1:]

            epo_path = epo_folder + '/' + subj
            gat_path = dec_folder + '/' + subj
            mkdir_p(gat_path)

            for session in session_nos[trial_type]:
                fname = epo_path + '/' + trial_type + session + '-epo.fif'

                epochs = read_epochs(fname)
                epochs.drop_channels(['EOG001','EOG003'])


                if epoch_params['savgol_hf'] is not None:
                    epochs.savgol_filter(epoch_params['savgol_hf'])

                epochs.crop(tmin=tmin, tmax=tmax)

                triggers = epochs.events[:,2]
                y_devVSstd_A = (triggers[np.in1d(triggers, 
                    (epochs.event_id['devA'],epochs.event_id['stdA']))] == \
                            epochs.event_id['devA']).astype(int)
                y_devVSstd_B = (triggers[np.in1d(triggers, 
                    (epochs.event_id['devB'],epochs.event_id['stdB']))] == \
                            epochs.event_id['devB']).astype(int)

                #### Use scaled Logistic Regression per default
                scaler = StandardScaler()
                clf = force_predict(LogisticRegression(penalty='l2', C=1), axis=1)
                    # C=1, solver='lbfgs', multi_class='multinomial'), axis=1) # use this for 3-class
                pipeline = Pipeline([('scaler', scaler), ('clf', clf)])
                # Define decoder. The decision_function is employed to use AUC for scoring
                gat = GeneralizationAcrossTime(n_jobs=4, clf=pipeline, scorer=auc_scorer,
                        predict_mode='mean-prediction') # must use mean pred when X-scoring

                # Define decoder. The decision_function is employed to use AUC for scoring
                #gat = GeneralizationAcrossTime(predict_mode='cross-validation', n_jobs=4)
                # If (clf is) None the classifier will be a standard pipeline including 
                # StandardScaler and a linear SVM with default parameters.
                ###########

                figs = []
                # fit and score
                gat.fit(epochs['devA','stdA'], y_devVSstd_A)
                gat.score(epochs['devB','stdB'], y_devVSstd_B)

                figs.append(gat.plot(vmin=0.2, vmax=0.8,
                         title="Generalization Across Time (dev vs std, A to B)"))
                figs.append(gat.plot_diagonal(chance=0.5))  # plot decoding across time (correspond to GAT diagonal)

                with open(gat_path + '/%s%s-devVSstd-AtoB-GAT.pickle' % (trial_type, session), 
                                        'wb') as f:
                    pickle.dump(gat, f, protocol=2) # use optimised binary format

                captions = [subj,subj] 
                sections = ['GAT-ses%s' % session,'Class-ses%s' % session]

                print 'Generating plots for', subj
                for ff,fig in enumerate(figs):
                    report.add_figs_to_section(fig, captions=captions[ff],
                        section=sections[ff],
                        scale=None, image_format='png')
                    plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)
                    

if do_GAT_groupstat:

    contrast_list = [\
            dict(gat_name='FFA, scaled LR', fname='FFA-GAT-scaledLR'),
            dict(gat_name='N2pc in VS1', fname='N2pc1-GAT'),
            #dict(gat_name='VS1 any target', fname='VS1-anyTRG-GAT'),
            dict(gat_name='FB1 any target', fname='FB1-anyTRG-GAT'),
            dict(gat_name='FB2 any target', fname='FB2-anyTRG-GAT'),
            dict(gat_name='FB1-STD CS+ vs. CS-', fname='FB1-identityCS-std-GAT'),
            dict(gat_name='FB1-DEV CS+ vs. CS-', fname='FB1-identityCS-dev-GAT'),
            dict(gat_name='FB2-STD CS+ vs. CS-', fname='FB2-identityCS-std-GAT'),
            dict(gat_name='FB2-DEV CS+ vs. CS-', fname='FB2-identityCS-dev-GAT'),
            #dict(gat_name='FB1 (dev vs std), A-to-B', fname='FB1-devVSstd-AtoB-GAT'),
            #dict(gat_name='FB2 (dev vs std), A-to-B', fname='FB2-devVSstd-AtoB-GAT'),
            ]


    tmin, tmax = -0.1, 0.35

    gat_rep_folder = rep_folder
    mkdir_p(gat_rep_folder)

    rep_file = gat_rep_folder + '/' + 'GAT_groupstat.html'

    report = Report(info_fname=None, subjects_dir=None, subject=None,
                    title='Generalization Across Time cluster stats',
                    verbose=None)

    included_subjects = db.get_subjects()
    #for subj in ['030_WAH']:
    for cont in contrast_list:
        gat_scores_list = []
        gat_mean = None
        gat_sem = None
        for subj in included_subjects:
            if len(subj) == 8:
                subj = subj[1:]

            gat_path = dec_folder + '/' + subj

            print "Reading GAT pickle for", subj
            with open(os.path.join(gat_path, cont['fname']+'.pickle'), 'rb') as f:
                gat = pickle.load(f)
                gat_scores_list.append(gat.scores_)
                if gat_mean is None:
                    gat_mean = copy.deepcopy(gat)
                    gat_sem = copy.deepcopy(gat)

        # Gather all scores
        scores = np.array(gat_scores_list)
        gat_mean.scores_ = np.mean(scores, axis=0)

        gat_sem.scores_ = np.std(scores, axis=0) / np.sqrt(len(included_subjects))
         
        # STATS
        chance = 0.5  # chance level; if it's an AUC, it has to be .5
        alpha = 0.05

        T_obs_, clusters, p_values, _ = spatio_temporal_cluster_1samp_test(
                    scores - chance, out_type='mask', n_permutations=128,
                        threshold=dict(start=2, step=2.), n_jobs=4)
           
        p_values = p_values.reshape(scores.shape[1:])
           
        figs = []
        # PLOT
        fig = gat_mean.plot(show=False, vmin=0.2, vmax=0.8, title=cont['gat_name'])
        ax = fig.axes[0]
        xx, yy = np.meshgrid(gat_mean.train_times_['times'],
                                     gat_mean.test_times_['times'][0],
                                                  copy=False, indexing='xy')
        ax.contour(xx, yy, p_values < alpha, colors='black', levels=[0])

        figs.append(fig)
        
        figd = gat_mean.plot_diagonal(chance=0.5)
        ax = figd.axes[0]
        gat_mean.scores_ -= gat_sem.scores_ # draw negative first
        gat_mean.plot_diagonal(ax=ax, color='r', chance=0.5)
        gat_mean.scores_ += 2. * gat_sem.scores_ # then the positive!
        gat_mean.plot_diagonal(ax=ax, color='r', chance=0.5)

        figs.append(figd)

        captions = ['TFCE, alpha=%g' % (alpha), 'Diag, chance=%g' % (chance)]

        for ff,fig in enumerate(figs):
            report.add_figs_to_section(fig, captions=captions[ff],
                section=cont['gat_name'],
                scale=None, image_format='png')
            plt.close(fig)

    report.save(fname=rep_file, open_browser=False, overwrite=True)




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

    from mayavi import mlab
    # need to run offscreen on isis (VNC)
    mlab.options.offscreen = True
    # This also works for running from terminal, but turns out not to be necc.
    # Note the need to specify the server number: default is 99 :)
    # xvfb-run -n 98 --server-args="-screen 0 1024x768x24" python evoked_pipeline.py

    #methods = ['MNE','dSPM']
    # Don't run MNE, just use dSPM for visualization
    methods = ['dSPM',]
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)
    trial_type = 'FFA'
    session = ''
    do_evoked_contrasts = {'face': True, 'diff': True}

    rep_folder = rep_folder + '/plot_STC_FFA/'
    mkdir_p(rep_folder)

    tmp_folder = ad._scratch_folder + '/tmp/'
    tmp_file_suffix = '.brain-tf_%02d.png'
    brain_times = np.array([60., 80., 100., 120., 140., 160., 180.,200.])
    use_abs_idx = False # Just use increments
    # found lh, then rh = 180 - az(lh)
    views = dict(
            lh={ # NB: swapping lat and med to make prettier plots!
                'lat': dict(azimuth=-40.,  elevation=130.),
                'med': dict(azimuth=-123., elevation=100.)},
            rh={
                'med': dict(azimuth=220., elevation=130.),
                'lat': dict(azimuth=303., elevation=100.)})

    stcran = dict(MNE={'max': 0.9, 'min': 0.1},
            dSPM={'max': 0.8, 'min': 0.2})

    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        stc_path = stc_folder + '/' + subj
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
                    fig = mlab.figure(size=(400,350))
                    #fig = mlab.figure(size=(400, 400))
                    brain = stc.plot(surface='inflated', hemi=hemi,
                            subject=subj, alpha = 0.9,
                            subjects_dir=fs_subjects_dir,
                            fmin=fmin, fmid=fmid, fmax=fmax,
                            figure=fig)
                                
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
                    #montage = [['lat', 'med'],['cau','ven']]
                    montage = [views[hemi]['med'], views[hemi]['lat']]

                    brain.save_image_sequence(time_idx, tmp_pattern,
                            use_abs_idx=use_abs_idx, montage=montage)

                    mlab.close(fig)

                for ii,tt in enumerate(brain_times):
                    cmd = 'montage -geometry 640x480+4+4 '
                    #cmd = 'montage -geometry +4+4 '
                    for hemi in ['lh', 'rh']:
                        cmd += tmp_folder + hemi + tmp_file_suffix % (ii) + ' '
                    tmpname = tmp_folder + 'both' + tmp_file_suffix % (ii)
                    cmd +=  tmpname

                    proc = subprocess.Popen([cmd], shell=True)
                    proc.communicate()

                    caption = method + ' @ %.0fms' % (tt)
                    report.add_images_to_section(tmpname, captions=caption,
                            section=cond, scale=None)
        report.save(fname=rep_file, open_browser=False, overwrite=True)

if do_STC_N2pc:
    # looking at the evokeds, it seems there's plenty to
    # see even efter 200, probably even longer.
    time_range = (-0.100, 0.360)
    methods = ['MNE','dSPM']
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)

    trial_type = 'VS'
    sessions = ['1','2']
    contrast_name = 'N2pc' # = trial_type? inverse operator taken for trial_type!
    do_evoked_contrasts = {'diff': True, 'diffA': True, 'diffB': True,
                            'devLH': True, 'devRH': True} 

    # NB, assuming really poor SNR!
    SNRs = {'diff': 1., 'diffA': 1., 'diffB': 1., 'devLH': 1., 'devRH': 1.}

    for subj in db.get_subjects():
        if len(subj) == 8:
            subj = subj[1:]

        evo_path = evo_folder + '/' + subj
        opr_path = opr_folder + '/' + subj
        stc_path = stc_folder + '/' + subj

        for session in sessions:
            # NB: contrast vs trial type
            evo_file = evo_path + '/' + contrast_name + session + '-avg.fif'
            inv_file = opr_path + '/' + trial_type + session + \
                    '-' + fwd_params['spacing'] + '-inv.fif'

            print 'Loading inverse operator...'
            inv_opr = read_inverse_operator(inv_file, verbose=False)

            for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
                # Load data
                evoked = read_evokeds(evo_file, condition=cond, verbose=False)
            
                lambda2 = 1. / SNRs[cond] ** 2.
                stc_path_SNR = stc_path + '/SNR%.0f' % (SNRs[cond])
                mkdir_p(stc_path_SNR)

                for method in methods:
                    # Save result in stc files
                    stc_file = stc_path_SNR + '/' + contrast_name + session + \
                            '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                    if file_exists(stc_file) and not CLOBBER:
                        continue

                    #print 'Applying inverse with method:', method
                    stc = apply_inverse(evoked, inv_opr,
                            lambda2, method, pick_ori=ori_sel, verbose=False)
                    stc.crop(tmin=time_range[0], tmax=time_range[1]) # CROP

                    stc.save(stc_file, verbose=False)

############################
if plot_STC_N2pc:

    from mayavi import mlab
    # need to run offscreen on isis (VNC)
    mlab.options.offscreen = True
    # This also works for running from terminal, but turns out not to be necc.
    # Note the need to specify the server number: default is 99 :)
    # xvfb-run -n 98 --server-args="-screen 0 1024x768x24" python evoked_pipeline.py

    #methods = ['MNE','dSPM']
    # Don't run MNE, just use dSPM for visualization
    methods = ['dSPM',]
    ori_sel = None # 'normal' leads to the SIGN of the estimates remaining (not good!)
    trial_type = 'VS'
    contrast_name = 'N2pc'
    sessions = ['1','2']
    do_evoked_contrasts = {'diff': True, 'devLH': True, 'devRH': True}
    # NB, these must be the same as when generated!!
    SNRs = {'diff': 1., 'diffA': 1., 'diffB': 1., 'devLH': 1., 'devRH': 1.}

    rep_folder = rep_folder + '/plot_STC_N2pc/'
    mkdir_p(rep_folder)

    tmp_folder = ad._scratch_folder + '/tmp/'
    tmp_file_suffix = '.brain-tf_%02d.png'
    brain_times = np.arange(180., 350., 20.)
    use_abs_idx = False # Just use increments
    # found lh, then rh = 180 - az(lh)
    views = dict(
            lh={ # NB: swapping lat and med to make prettier plots!
                'caulo': dict(azimuth=-80., elevation=120.),
                'lat': dict(azimuth=-40.,  elevation=130.),
                'med': dict(azimuth=-123., elevation=100.)},
            rh={
                'caulo': dict(azimuth=-100., elevation=120.),
                'med': dict(azimuth=220., elevation=130.),
                'lat': dict(azimuth=303., elevation=100.)},
            both={
                'caulo': dict(azimuth=-90., elevation=120.),
                })

    #stcran = dict(MNE={'max': 0.9, 'min': 0.1},
    #        dSPM={'max': 0.8, 'min': 0.2})
    stc_clim = dict(kind='percent', lims=(90.,98.,100.))

    for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
        # Load data
        for method in methods:
            for session in sessions:

                rep_file = rep_folder + '/' + contrast_name + session + '-allsubs_' + \
                                cond + '-' + method + '-SNR%.0f'%(SNRs[cond]) + '.html'

                #  cannot be loaded/appended :(
                report = Report(info_fname=None, 
                        subjects_dir=fs_subjects_dir, subject='VSaverage',
                        title='N2pc estimates', verbose=None)

                for subj in db.get_subjects():
                    if len(subj) == 8:
                        subj = subj[1:]

                    stc_path = stc_folder + '/' + subj+ '/SNR%.0f' % (SNRs[cond])
                    # Load results from stc files
                    stc_file = stc_path + '/' + contrast_name + session + \
                            '-' + fwd_params['spacing'] + '_' + cond + '_' + method

                    stc = read_source_estimate(stc_file)

                    #fmax = stcran[method]['max']*np.ravel(stc.data).max()
                    #fmin = stcran[method]['min']*fmax
                    #fmid = (fmax - fmin) / 2.

                    for hemi in ['lh','rh']:
                        print 'Hemi :', hemi
                        fig = mlab.figure(size=(400,350))
                        #fig = mlab.figure(size=(400, 400))

                        brain = stc.plot(surface='inflated', hemi=hemi,
                                subject=subj, alpha = 0.9,
                                subjects_dir=fs_subjects_dir,
                                clim=stc_clim,
                                figure=fig)
                                    
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
                        #montage = [['lat', 'med'],['cau','ven']]
                        #montage = [views[hemi]['caulo'],]
			montage='current'
			mlab.view(views[hemi]['caulo']['azimuth'],
				views[hemi]['caulo']['elevation'])

                        brain.save_image_sequence(time_idx, tmp_pattern,
                                use_abs_idx=use_abs_idx, montage=montage)

                        mlab.close(fig)

                    for ii,tt in enumerate(brain_times):
                        cmd = 'montage -geometry 640x480+4+4 '
                        #cmd = 'montage -geometry +4+4 '
                        for hemi in ['lh', 'rh']:
                            cmd += tmp_folder + hemi + tmp_file_suffix % (ii) + ' '
                        tmpname = tmp_folder + 'both' + tmp_file_suffix % (ii)
                        cmd +=  tmpname

                        proc = subprocess.Popen([cmd], shell=True)
                        proc.communicate()

                        caption = subj + '-' + method
                        secname = '%.0fms' % (tt)

                        report.add_images_to_section(tmpname, captions=caption,
                                section=secname, scale=None)

            # at level with sessions
            report.save(fname=rep_file, open_browser=False, overwrite=True)

if do_STC_FFA_groupavg:

    from mayavi import mlab
    # need to run offscreen on isis (VNC)
    mlab.options.offscreen = True

    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_to = 'VSaverage'
    methods = ['MNE','dSPM']

    trial_type = 'FFA'
    sessions = ['',]
    contrast_name = 'FFA' # = trial_type? inverse operator taken for trial_type!
    do_evoked_contrasts = {'diff': True, 'face': True}
    SNRs = {'face': 3., 'diff': 3.}

    rep_folder = rep_folder + '/plot_STC_FFA_groupave/'
    mkdir_p(rep_folder)

    tmp_folder = ad._scratch_folder + '/tmp/'
    tmp_file_suffix = '.brain-tf_%02d.png'
    brain_times = np.arange(80., 220., 20.)
    use_abs_idx = False # Just use increments
    # found lh, then rh = 180 - az(lh)
    views = dict(
            lh={ # NB: swapping lat and med to make prettier plots!
                'caulo': dict(azimuth=-70., elevation=110.),
                'lat': dict(azimuth=-40.,  elevation=130.),
                'med': dict(azimuth=-123., elevation=100.)},
            rh={
                'caulo': dict(azimuth=-110., elevation=110.),
                'med': dict(azimuth=220., elevation=130.),
                'lat': dict(azimuth=303., elevation=100.)},
            both={
                'caulo': dict(azimuth=-90., elevation=120.),
                })

    #stcran = dict(MNE={'max': 0.9, 'min': 0.1},
    #        dSPM={'max': 0.8, 'min': 0.2})
    stc_clim = dict(kind='percent', lims=(90.,98.,100.))

    included_subjects = db.get_subjects()

    for session in sessions:
        for method in methods:
            rep_file = rep_folder + '/' + contrast_name + session + \
                    '-VSaverage-' + method + '.html'

            #  cannot be loaded/appended :(
            report = Report(info_fname=None, 
                    subjects_dir=fs_subjects_dir, subject='VSaverage',
                    title='FFA estimates', verbose=None)

            for cond in \
                    [k for k in do_evoked_contrasts.keys() if \
                    do_evoked_contrasts[k]]:

                stc_list = [] # for holding the stc before averaging

                ave_stc_path = stc_folder + '/VSaverage'
                ave_stc_file = ave_stc_path + '/' + contrast_name + session + \
                        '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                if file_exists(ave_stc_file+'-lh.stc') and not CLOBBER:
                    print "Average already exists & CLOBBER is False"
                    print 'Reading average %s -> %s -> %s' % (contrast_name, cond, method)
                    stc_ave = mne.read_source_estimate(ave_stc_file)
                else:
                    for subj in included_subjects:
                        if len(subj) == 8:
                            subj = subj[1:]

                        opr_path = opr_folder + '/' + subj
                        stc_path = stc_folder + '/' + subj + '/SNR%.0f' % (SNRs[cond])

                        # Load stc file
                        stc_file = stc_path + '/' + contrast_name + session + \
                                '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                        stc_from = mne.read_source_estimate(stc_file)

                        print 'Morphing', subj, 'to', subject_to
                        stc_to = mne.morph_data(subj, subject_to,
                                stc_from, grade=vertices_to, n_jobs=4, verbose=False)

                        stc_list.append(stc_to)

                    print 'Computing average...'
                    stc_ave = reduce(add, stc_list)
                    stc_ave /= len(stc_list)

                    print 'Saving average %s -> %s -> %s' % (contrast_name, cond, method)
                    stc_ave.save(ave_stc_file, verbose=False)

                for hemi in ['lh','rh']:
                    print 'Hemi :', hemi
                    fig = mlab.figure(size=(400,350))
                    #fig = mlab.figure(size=(400, 400))

                    brain = stc_ave.plot(surface='inflated', hemi=hemi,
                            subject='VSaverage', alpha = 0.9,
                            subjects_dir=fs_subjects_dir,
                            clim=stc_clim,
                            views=[views[hemi]['caulo']],
                            figure=fig)
                                
                    brain.add_label("Pole_occipital", color='springgreen',
                            borders=False, alpha=0.2)
                    brain.add_label("S_calcarine", color='aquamarine',
                            borders=False, alpha=0.2)
                    brain.add_label("G_oc-temp_lat-fusifor", color='aquamarine',
                            borders=False, alpha=0.2)
                    brain.add_label("Pole_occipital", color='springgreen',
                            borders=True, alpha=1.)
                    brain.add_label("S_calcarine", color='aquamarine',
                            borders=True, alpha=1.)
                    brain.add_label("G_oc-temp_lat-fusifor", color='aquamarine',
                            borders=True, alpha=1.)

                    time_idx = [brain.index_for_time(t) for t in brain_times]

                    tmp_pattern = tmp_folder + hemi + tmp_file_suffix
                    #montage = [['lat', 'med'],['cau','ven']]
                    montage = [views[hemi]['med'], views[hemi]['lat']]
                    montage='current'

                    brain.save_image_sequence(time_idx, tmp_pattern,
                            use_abs_idx=use_abs_idx, montage=montage)

                    mlab.close(fig)

                for ii,tt in enumerate(brain_times):
                    cmd = 'montage -geometry 640x480+4+4 '
                    #cmd = 'montage -geometry +4+4 '
                    for hemi in ['lh', 'rh']:
                        cmd += tmp_folder + hemi + tmp_file_suffix % (ii) + ' '
                    tmpname = tmp_folder + 'both' + tmp_file_suffix % (ii)
                    cmd +=  tmpname

                    proc = subprocess.Popen([cmd], shell=True)
                    proc.communicate()

                    secname = cond
                    caption = '%.0fms' % (tt)

                    report.add_images_to_section(tmpname, captions=caption,
                            section=secname, scale=None)

            # at level with condition
            report.save(fname=rep_file, open_browser=False, \
                    overwrite=True)


if do_STC_N2pc_groupavg:
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_to = 'VSaverage'
    methods = ['MNE','dSPM']

    trial_type = 'VS'
    sessions = ['1','2']
    contrast_name = 'N2pc' # = trial_type? inverse operator taken for trial_type!
    do_evoked_contrasts = {'diff': True, 'diffA': True, 'diffB': True,
                            'devLH': True, 'devRH': True} 

    included_subjects = db.get_subjects()

    for session in sessions:
        for cond in [k for k in do_evoked_contrasts.keys() if do_evoked_contrasts[k]]:
            for method in methods:

                stc_list = [] # for holding the stc before averaging

                ave_stc_path = stc_folder + '/VSaverage'
                ave_stc_file = ave_stc_path + '/' + contrast_name + session + \
                        '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                if file_exists(ave_stc_file) and not CLOBBER:
                    print "Average already exists, skipping..."
                    continue

                for subj in included_subjects:
                    if len(subj) == 8:
                        subj = subj[1:]

                    opr_path = opr_folder + '/' + subj
                    stc_path = stc_folder + '/' + subj

                    # Load stc file
                    stc_file = stc_path + '/' + contrast_name + session + \
                            '-' + fwd_params['spacing'] + '_' + cond + '_' + method
                    stc_from = mne.read_source_estimate(stc_file)

                    print 'Morphing', subj, 'to', subject_to
                    stc_to = mne.morph_data(subj, subject_to,
                            stc_from, grade=vertices_to, n_jobs=4, verbose=False)

                    stc_list.append(stc_to)

                print 'Computing average...'
                stc_ave = reduce(add, stc_list)
                stc_ave /= len(stc_list)

                print 'Saving average %s -> %s -> %s' % (contrast_name, cond, method)
                stc_ave.save(ave_stc_file, verbose=False)


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

