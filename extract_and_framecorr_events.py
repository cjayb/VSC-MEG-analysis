"""
=========================
Find visual events in a raw file and correct for the frame delay 
measured using an optoprobe on MISC001
=========================

Find events from the stimulation/trigger channel in the raw data.
"""
# Author: Chris Bailey <cjb@cfin.au.dk>
#
# License: BSD (3-clause)
# Based on code from mne-python, A.Gramfort et al.

print __doc__

import mne
from mne.fiff import Raw, pick_channels
import numpy as np

def next_crossing(a,lowlim):
    for a_ind,val in enumerate(a):
        if val > lowlim:
            return a_ind
    
    print 'ERROR: Shouldn''t get this far!'
    return -1

def find_next_analogue_trigger(raw, ind, ana_channel='MISC001', 
                                lowlim=0.1, hilim=0.2):
                                
    #print 'Nothing here yet'
    pick = pick_channels(raw.info['ch_names'], include=ana_channel)
    ana_data, _ = raw[pick,ind:ind+200]
    return next_crossing(ana_data[0,:].squeeze(),lowlim)
    
    

projname = 'MINDLAB2013_01-MEG-AttentionEmotionVisualTracking'
subjid_num = 10
subjid_str = '%04d' % subjid_num
data_path = '/projects/' + projname + '/raw/' + subjid_str

hard_coded_file = '/20130910_000000/MEG/002.VS_1b_1/files/PROJ0103_SUBJ0010_SER002_FILESNO001.fif'
fname = data_path + hard_coded_file
eve_name = data_path + hard_coded_file + '-correve.fif'

# Reading events
raw = Raw(fname)

events = mne.find_events(raw, stim_channel='STI101')
# This returns, oddly the sample indices to which
# raw.first_samp has been ADDED! So to get the actual index
# into raw[:,__], we need to subtract the first_samp!

corr_events = events.copy()

frame_delays = np.zeros(events.shape[0])
responses = np.zeros((192,4))

# These are the triggers that are followed by an image
# FB_neutral, FB_angry, no_target_S03, no_target_S10, targets at 6 positions
imtriggers = np.r_[10,11,100,200,np.arange(111,117), np.arange(211,217)]

# Scan through all events, even though not terribly pretty
# The list isn't that huge to warrant any coolness here
row=0
resp=0
prev_target = -1
prev_target_time = -1
prev_trigger = -1
for ind, before, after in events:

    
#    print "At sample %d stim channel went from %d to %d" % (
#                                            ind, before, after)

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
        ind -= raw.first_samp    # now these really are indices into raw!
        anatrig_ind = find_next_analogue_trigger(raw, ind, 
                                                 ana_channel='MISC001', 
                                                 lowlim=0.1, hilim=0.2)
        corr_ind += anatrig_ind
        corr_events[row,0] = corr_ind                                         
        frame_delays[row] = anatrig_ind
        
        if after != 10 and after != 11:
            if after == 100 or after == 200:
                prev_target = 2
            else:
                prev_target = 3
            prev_target_time = raw.index_as_time(corr_ind)
    
    prev_trigger = after                                             
    row += 1

#frame_delays[frame_delays < 0.1] = np.NaN
#masked_frame_delays = np.ma.masked_equal(frame_delays,np.NaN)
hist,bins = np.histogram(frame_delays,bins=[24.6, 41.3,58.3,74.0])
center=(bins[:-1]+bins[1:])/2
bar(center,hist)
#bar(bins,hist)

# Writing events
#mne.write_events(eve_name, events)

#for ind, before, after in events[:5]:
#    print "At sample %d stim channel went from %d to %d" % (
#                                                    ind, before, after)
