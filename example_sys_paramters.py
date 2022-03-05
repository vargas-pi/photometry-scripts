"""
an example of what to put in the sys_parameters.py file
these are settings for our systems. our systems have different naming
schemes for the gizmos and channels so the parameters for each scheme
are shown here. we have essentially this same file on all systems but just
uncomment the necessary parameters on each.

the important thing to note here is the structure of the dictionaries and their names
the number of the line is the first level of keys in all dictionaries and if relevant, 
the channels are the next level. 

NOTE: our systems are currently set up with 2 lines, each
with a 490nm and a 405nm fluorescence channel. As a result all of the code using these
parameters is written with this assumption for the keys in the dictionary. this should be
fixed in future versions of the code
"""

#old system parameters
DEFAULT_SAVE_CHANNELS = {1: {490: 'x19A',
                             405: 'x15A'}, 
                         2: {490: 'x29B', 
                             405: 'x25B'}
                        }
GIZMOS = {1: 'FibPho1', 
          2: 'FibPho1'}

CHANNELS = {1: {490: 'Response-1A', 
                405: 'Response-2A'},
            2: {490: 'Response-1B', 
                405: 'Response-2B'}
            }

#new system parameters
# DEFAULT_SAVE_CHANNELS = {1: {490: 'x19A',
#                              405: 'x15A'}, 
#                          2: {490: 'x29B', 
#                              405: 'x25B'}
#                         }
# GIZMOS = {1: 'FibPho1', 2: 'FibPho2'}
# CHANNELS = {1: {490: 'Response-2A', 
#                 405: 'Response-1A'},
#             2: {490: 'Response-2A', 
#                 405: 'Response-1A'}
#             }