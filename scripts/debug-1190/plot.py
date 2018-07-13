#!/usr/bin/python3.6

import sys, pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if (sys.argv[1] != 'hms' and sys.argv[1] != 'shms') :
    print ('Usage: python debug1190.py hms (shms)')
    sys.exit(1)

if (sys.argv[1] == 'hms')  : 
    spec   = 'h'
    slno   = [2, 4, 5, 7,  8, 10, 13, 14, 16, 17]
if (sys.argv[1] == 'shms') : 
    spec   = 'p'
    slno   = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

edf = pd.read_pickle(spec+'edf.pkl')
bdf = pd.read_pickle(spec+'bdf.pkl')

plt.figure('Bunch ID for All Slots')
for slot in slno :
    #if (slot == 4 or slot == 16) :
    plt.plot(edf['slot-%d' % slot]['chip-0'], bdf['slot-%d' % slot]['chip-0'], label = 'Slot %d' %slot)
    plt.legend(loc = 'best')
    plt.title('Chip-0 Bunch ID')
    plt.xlabel('Event Number')
    plt.ylabel('Bunch ID')
# plt.show()

if (sys.argv[1] == 'hms')  : 
    for index in range(len(slno)) :
        #if (index == 0) : continue
        plt.figure('S%d - All Slots' % slno[index])
        for slot in slno :
            #if (slot == 2) : continue
            #if (slot == 4 or slot == 8 or slot == 13 or slot == 16) : 
            plt.plot(edf['slot-%d' % slot]['chip-0'], bdf['slot-%d' % slno[index]]['chip-0'] - bdf['slot-%d' % slot]['chip-0'], 
                     label = 'S%d - S%d' % (slno[index], slot))
            plt.legend(loc = 'best')
            plt.title('Chip-0 Bunch ID (S%d - All Other Slots)' % slno[index])
            plt.xlabel('Event Number')
            plt.ylabel('Bunch ID')
            #plt.savefig('/home/pooser/Desktop/slot-%d-bid-diff.pdf' % slno[index])

# if (sys.argv[1] == 'shms')  : 
#     plt.figure()
#     for slot in slno :
#         plt.plot(edf['slot-%d' % slot]['chip-0'], bdf['slot-6']['chip-0'] - bdf['slot-%d' % slot]['chip-0'], 
#                  label = 'S6 - S%d' %slot)
#         plt.legend(loc = 'best')

#plt.show()
