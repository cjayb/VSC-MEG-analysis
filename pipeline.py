# 
#
# License: BSD (3-clause)
import matplotlib
matplotlib.use('agg') # force non-interactive plotting
import numpy as np
import os, errno
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
do_simple_contrasts_univar = True

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
                
                    
if do_simple_contrasts_univar: # do a couple of "main effects"

    import matplotlib.pyplot as plt
    clim_all = dict(mag=[-250, 250], grad=[0, 50])
    clim_con = dict(mag=[-125, 125], grad=[0, 25])
    topo_times = np.arange(0.0, 0.210,0.020)
    #for subj in ['007_SGF']:
    for subj in ad.analysis_dict.keys():

        epo_path = ad._scratch_folder + '/epochs/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj
        img_path = ad._scratch_folder + '/epochs/' + filt_dir + '/' + filter_params['input_files'] + '/' + subj + '/img'
        mkdir_p(img_path)

        for session in ['pre','post']:
            for trial_type in ['VS','FB']:
                fname = epo_path + '/' + trial_type + '_' + session + '-epo.fif' 
                epochs = mne.read_epochs(fname)

                evoked_all = epochs[['stdB','devB','stdA','devA']].average()
                evoked_face = epochs[['stdB','devB']].average() - epochs[['stdA','devA']].average()
                evoked_odd = epochs[['devA','devB']].average() - epochs[['stdA','stdB']].average()
                cov = mne.compute_covariance(epochs, tmin=baseline[0], tmax=baseline[1]) # same covariance for all contrasts

                evoked_all.plot_image(clim=clim_all, show=False)
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_allERF_time.png')
                plot_evoked_topomap(evoked_all,topo_times, show=False, vmin=[clim_all['grad'][0],clim_all['mag'][0]], vmax=[clim_all['grad'][1],clim_all['mag'][1]])
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_allERF_topo.png')
                evoked_face.plot_image(clim=clim_con, show=False)
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_faceERF_time.png')
                plot_evoked_topomap(evoked_face,topo_times, show=False, vmin=[clim_con['grad'][0],clim_con['mag'][0]], vmax=[clim_con['grad'][1],clim_con['mag'][1]])
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_faceERF_topo.png')
                evoked_odd.plot_image(clim=clim_con, show=False)
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_oddERF_time.png')
                plot_evoked_topomap(evoked_odd,topo_times, show=False, vmin=[clim_con['grad'][0],clim_con['mag'][0]], vmax=[clim_con['grad'][1],clim_con['mag'][1]])
                plt.savefig(img_path + '/' + trial_type + '_' + session + '_oddERF_topo.png')
                #eve_dict, con_dict, con_names = events_logic(events, contrast) # not efficient but will do

                # average epochs and get an Evoked dataset.
                #evoked = epochs[con_names[0]].average() -  \
                #            epochs[con_names[1]].average() 
                        
                #evo_out = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-ave.fif'
                #cov_out = evo_path + '/' + trial_type + '_' + contrast + '_' + session + '-cov.fif'
                #print evo_out
                #evoked.save(evo_out)  # save evoked data to disk
                #cov.save(cov_out)  # save evoked data to disk

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

