#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#repFile = np.empty(5, dtype = str)
repFile = np.array([])
repFilePrefix = '/Users/pooser/hallc_replay/REPORT_OUTPUT/SHMS/PRODUCTION/replay_shms_all_production_'
repFileSuffix = '_-1.root'

lh2_21d_5p1 = { 'runRange' : np.arange(2484, 2488+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

index = 0
for run in lh2_21d_5p1['runRange']:
    print repFile
    print run
    #lh2_21d_5p1.update( { 'repFiles' : repFile.append(repFilePrefix + str(run) + repFileSuffix) } )
    lh2_21d_5p1.update( { 'repFiles' : np.append(repFile, repFilePrefix + str(run) + repFileSuffix) } )
    #lh2_21d_5p1.update( { 'repFiles' : repFile.insert(index, repFilePrefix + str(run) + repFileSuffix) } )
    print lh2_21d_5p1['repFiles']
    index += 1

print lh2_21d_5p1['repFiles']

ald_21d_5p1 = { 'runRange' : np.arange(2489, 2489+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

ld2_21d_5p1 = { 'runRange' : np.arange(2490, 2494+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

b10_21d_5p1 = { 'runRange' : np.arange(2495, 2496+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

b11_21d_5p1 = { 'runRange' : np.arange(2497, 2498+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

be_21d_5p1  = { 'runRange' : np.arange(2502, 2505+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }

c_21d_5p1   = { 'runRange' : np.arange(2507, 2508+1),
                'pCentral' : -5.1,
                'theta'    : 21.035 }




        
