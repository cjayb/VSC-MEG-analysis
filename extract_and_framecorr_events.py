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
from mne.fiff import Raw

def find_next_analogue_trigger(raw, ind, ana_channel='MISC001', 
                                lowlim=1, hilim=3):
                                
    


data_path = '/projects/'
hard_coded_file = ''
fname = data_path + hard_coded_file
eve_name = data_path + hard_coded_file + '-correve.fif'

# Reading events
raw = Raw(fname)

events = mne.find_events(raw, stim_channel='STI101')

# These are the triggers that are followed by an image
# FB_neutral, FB_angry, no_target_S03, no_target_S10, targets at 6 positions
imtriggers = (10,11,100,200,range(111,117), range(211,217))

# Scan through all events, even though not terribly pretty
# The list isn't that huge to warrant any coolness here
for ind, before, after in events:
    print "At sample %d stim channel went from %d to %d" % (
                                            ind, before, after)
                                            
    # Now check whether to search for the frame trigger
    # Also: since we'll want to fix the reaction times, search also for 
    # the response trigger
    
    find_next_analogue_trigger(raw, ind, ana_channel='MISC001', 
                                lowlim=1, hilim=3)


# Writing events
mne.write_events(eve_name, events)

for ind, before, after in events[:5]:
    print "At sample %d stim channel went from %d to %d" % (
                                                    ind, before, after)
