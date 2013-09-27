# -*- coding: utf-8 -*-
"""
Due to delay and JITTER of the presentation of visual frames, event files must
be created and corrected for. This simultaneously gives us corrected reaction
times. 

Output: 
scratch/events.fif/SUB_ID/task-correve.fif
scratch/RTs.csv/SUB_ID/task-corrRT.csv

Created on Thu Sep 26 10:25:28 2013

@author: cjb
"""

from database import Query
from analysis_dict import Anadict

from mne import find_events, write_events
from mne.fiff import Raw, pick_channels

import numpy as np

import csv

from sys import exit as sysexit
import os
import errno

def check_path_exists(chkpath):
    
    try: 
        os.makedirs(chkpath)        
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(chkpath):
            pass
        else:
            raise

def write_RTs(fname, data):
    
    fp = open(fname, 'w')
    writ = csv.writer(fp)  
    writ.writerows(data)
    fp.close()
    
def write_frame_delays(fname, data):
    
    fp = open(fname, 'w')
    for d in data:    
        fp.write("%.0f\n" % d)    
    fp.close()

def _next_crossing(a,lowlim):
    for a_ind,val in enumerate(a):
        if val > lowlim:
            return a_ind
    
    print 'ERROR: No analogue trigger found within %d samples of the digital trigger' % a_ind
    print 'Cannot continue, aborting...'
    sysexit(-1)

#def _find_next_analogue_trigger(raw, ind, lim, ana_channel='MISC001', 
#                                maxdelay_samps=200):
                                
    #print 'Nothing here yet'
    #pick = pick_channels(raw.info['ch_names'], include=ana_channel)
    #ana_data, _ = raw[pick,ind:ind+maxdelay_samps]
def _find_next_analogue_trigger(ana_data, ind, lim, maxdelay_samps=100):

    return _next_crossing(ana_data[ind:ind+maxdelay_samps].squeeze(),lim)


def _find_analogue_trigger_limit(ana_data):
    
    return ana_data.mean()

proj_code = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'

stimchan = 'STI101'
anachan = 'MISC001'

# These are the triggers that are followed by an image
# FB03_neutral, FB03_angry,FB10_neutral, FB10_angry, 
# no_target_S03, blurred face (FFA only) no_target_S10, targets at 6 positions
imtriggers = np.r_[10,11,20,21,100,150,200,np.arange(111,117), np.arange(211,217)]

db = Query(proj_code=proj_code)
anadict = Anadict(db)    
# I sent an empty room file from another subject in to 007, but it screwed up the
# series numbering! On 26 Sep I /really/ messed it up and have to get help from
# Jesper to fix it! This script ran with the below line (before major cockup),
# so assume it is OK and that 007 will be fixed.
## anadict.analysis_dict['007_SGF']['001.VS_1a_1']['raw'] = '/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/raw/0007/20130910_000000/MEG/001.VS_1a_1/files/PROJ0103_SUBJ0007_SER001_FILESNO001.fif'

for subj in anadict.analysis_dict.keys():
    for task in anadict.analysis_dict[subj].keys():

        
        if 'empty' in task:
            pass # this correction only applies to the actual paradigm
        else:
        #elif '003' in task:
            raw_name = anadict.analysis_dict[subj][task]['raw']

            # Reading events
            raw = Raw(raw_name, preload=False)

            events = find_events(raw, stim_channel=stimchan, min_duration=0.002)
            # This returns, oddly the sample indices to which
            # raw.first_samp has been ADDED! So to get the actual index
            # into raw[:,__], we need to subtract the first_samp!

            corr_events = events.copy()

            frame_delays = np.zeros(events.shape[0])
            
            if not 'FFA' in task:
                # Get total number of responses (should be 192!)
                codes = events[:,2]
                lefts = codes==2
                rights = codes==3
                resps = lefts + rights
                responses = np.zeros((resps.sum(),4))
                if resps.sum() != 192:
                    print 'Warning: Total number of responses is not 192! (%d)' % resps.sum()

            print "Loading analogue channel data and deciding on threshold..."
            pick = pick_channels(raw.info['ch_names'], include=anachan)
            ana_data, _ = raw[pick,:]
            ana_data = ana_data[0]
            triglimit = _find_analogue_trigger_limit(ana_data)
            print "Analogue data trigger limit set to %.2f" % triglimit


            # Scan through all events, even though not terribly pretty
            # The list isn't that huge to warrant any coolness here
            row=0
            resp=0
            prev_target = -1
            prev_target_time = -1
            prev_trigger = -1
            for ind, before, after in events:
            
                # Now check whether to search for the frame trigger
                # Also: since we'll want to fix the reaction times, search also for 
                # the response strigger
            
                corr_ind = ind  
                
                if after == 2 or after == 3:
                    RT = raw.index_as_time(ind) - prev_target_time
                    if prev_target == after:
                        responses[resp,:] = np.r_[resp+1,True,prev_trigger,RT]
                    else:
                        responses[resp,:] = np.r_[resp+1,False,prev_trigger,RT]
                        
                    resp += 1
                
                elif after in imtriggers:
                    raw_ind = ind - raw.first_samp    # now these really are indices into raw!
#                    anatrig_ind = _find_next_analogue_trigger(raw, raw_ind,triglimit 
#                                                             ana_channel='MISC001', 
#                                                             maxdelay_samps=200)
                    anatrig_ind = _find_next_analogue_trigger(ana_data, raw_ind,triglimit ,
                                                             maxdelay_samps=100)
                    corr_ind += anatrig_ind
                    corr_events[row,0] = corr_ind                                         
                    frame_delays[row] = anatrig_ind
                    
                    if not 'FFA' in task:    
                        # Check that this is not a feedback pic!
                        if after != 10 and after != 11 and after != 20 and after != 21:
                            if after == 100 or after == 200:
                                prev_target = 2
                            else:
                                prev_target = 3
                            prev_target_time = raw.index_as_time(corr_ind)
                    
                prev_trigger = after                                             
                row += 1
            

            corr_events_path = db._scratch_folder + '/events.fif/' + subj
            check_path_exists(corr_events_path)
            corr_events_fname = corr_events_path + '/' + task[4:] + '-correve.fif'
            write_events(corr_events_fname, corr_events)

            frame_delays_fname = corr_events_path + '/' + task[4:] + '-frame_delays.txt'
            write_frame_delays(frame_delays_fname, frame_delays)
            
            if not 'FFA' in task: 
                responses_path = db._scratch_folder + '/RTs.csv/' + subj
                check_path_exists(responses_path)
                responses_fname = responses_path + '/' + task[4:] + '-corrRTs.csv'
                write_RTs(responses_fname, responses)


            dict_entry = anadict.analysis_dict[subj][task]
            dict_entry.update({'correve': corr_events_fname})
            dict_entry.update({'frame_delays': frame_delays_fname})
            if not 'FFA' in task: 
                dict_entry.update({'corrRT': responses_fname})
            

anadict.save('Added corrected event files, reaction times and frame delays')    