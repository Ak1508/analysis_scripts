#!/usr/bin/python

# Import various modules
import sys
import string
import numpy as np
import matplotlib.pyplot as plt
import ROOT as r
import rootpy as rp

from glob import *
from rootpy import *

# Acquire 
if (sys.argv[1] != 'hms' and sys.argv[1] != 'shms'):
    print 'Usage: python getCharge.py hms (or shms)'
    sys.exit(1)

# Define run list and report files
if (sys.argv[1] == 'hms') : 
    rf = glob('hms-reports/replay_hms_production_*_-1.report')
    df = glob('hms-data/hms_replay_production_*_-1.root')
if (sys.argv[1] == 'shms') : 
    rf = glob('shms-reports/replay_shms_production_*_-1.report')
    df = glob('shms-data/shms_replay_production_*_-1.root')

# Define dictionaries
# xem data dictionary
xem = {}
# Target dictionary
xem_tar = { 'lh2' : 1.01,
            'ld2' : 2.01,
            'be9' : 9.01,
            'b10' : 10.01,
            'b11' : 11.01,
            'c12' : 12.01,
            'ald' : 26.98 }
# Report file dictionary
xem_rf = {  'data'     : [],  # data file
            'rn'       : [],  # run number
            'pcent'    : [],  # central momentum
            'tamu'     : [],  # target amu
            'theta'    : [],  # spectrometer theta
            'ebeam'    : [],  # beam energy
            'i4a'      : [],  # bcm4a current (uA)
            'i4a_cut'  : [],  # bcm4a current (uA, cut > 5 uA)
            'q4a'      : [],  # bcm4a charge (mC)
            'q4a_cut'  : [],  # bcm4a charge (mC, cut > 5 uA)
            'clt'      : [],  # computer live time
            'elt'      : [],  # electronic live time
            'tr_eff'   : [],  # tracking efficiency
            'etr_eff'  : [],  # electron tracking efficiency
            'htr_eff'  : [],  # hadron tracking efficiency
            'scin_eff' : [] } # 3/4 trigger efficiency           

# Store values of interest in arrays
for index, run in enumerate(rf):
    xem_rf['data'].append(df[index])
    with open(rf[index]) as fobj:
        for line in fobj:
            data = line.split(':')
            #if ('' in data[0]) : xem_rf[''].append(data[1].strip())
            # Kinematic configurations
            if ('Run Num'     in data[0]) : xem_rf['rn'].append(data[1].strip())
            if ('Momentum'    in data[0]) : xem_rf['pcent'].append(data[1].strip())
            if ('Target AMU'  in data[0]) : xem_rf['tamu'].append(data[1].strip())
            if ('Spec Theta'  in data[0]) : xem_rf['theta'].append(data[1].strip())
            if ('Beam Energy' in data[0]) : xem_rf['ebeam'].append(data[1].strip())
            # Charge and current
            # if ('' in data[0]) : xem_rf[''].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Current' in data[0])          : xem_rf['i4a'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Beam Cut Current' in data[0]) : xem_rf['i4a_cut'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Charge' in data[0])           : xem_rf['q4a'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Beam Cut Charge' in data[0])  : xem_rf['q4a_cut'].append(filter(lambda x: x in string.digits + '.', data[1]))
            # Live times (must be multiplied by 0.01 -> done later)
            if (sys.argv[1] == 'hms')  : 
                if ('Pre-Scaled Ps2 HMS Computer Live Time' in data[0])  : xem_rf['clt'].append(data[1][:8].strip())
            if (sys.argv[1] == 'shms') : 
                if ('Pre-Scaled Ps2 SHMS Computer Live Time' in data[0]) : xem_rf['clt'].append(data[1][:8].strip())
            if ('OG 6 GeV Electronic Live Time (100, 150)' in data[0])   : xem_rf['elt'].append(data[1][:8].strip())
            # Tracking efficiencies
            if (sys.argv[1] == 'hms')  : 
                if ('SING FID TRACK EFFIC' in data[0])        : xem_rf['tr_eff'].append(data[1][:11].strip())
                if ('E SING FID TRACK EFFIC' in data[0])      : xem_rf['etr_eff'].append(data[1][:11].strip())
                if ('HADRON SING FID TRACK EFFIC' in data[0]) : xem_rf['htr_eff'].append(data[1][:11].strip())
            if (sys.argv[1] == 'shms')  : 
                if ('SING FID TRACK EFFIC' in data[0])        : xem_rf['tr_eff'].append(data[1][:8].strip())
                if ('E SING FID TRACK EFFIC' in data[0])      : xem_rf['etr_eff'].append(data[1][:8].strip())
                if ('HADRON SING FID TRACK EFFIC' in data[0]) : xem_rf['htr_eff'].append(data[1][:8].strip())
            # Trigger efficiency
            if ('3_of_4 EFF' in data[0]) : xem_rf['scin_eff'].append(data[1].strip())

# Enumerate xem targets
for tar_str, tar_amu in xem_tar.items():
    # Initialize dictionary
    xem[tar_str] = {}
    # Enumerate variables
    for var_str, var in xem_rf.items():
        # Initialize lists in xem dictionary
        xem[tar_str][var_str] = []
        # Enumerate target list from report files
        for index, target in enumerate(xem_rf['tamu']):
            # Append lists when enumerated targets are identical
            if (float(xem_tar[tar_str]) == float(target)):
                xem[tar_str][var_str].append(xem_rf[var_str][index])

# Convert lists to arrays, kludgy as hell but oh well
for tar, tar_dict in xem.items():
    for rf_vars, rf_list in xem[tar].items():
        if (rf_vars == 'data') : continue
        rf_array = np.asarray(rf_list, dtype = float)
        del xem[tar][rf_vars]
        if   (rf_vars == 'clt') : xem[tar][rf_vars] = rf_array*0.01
        elif (rf_vars == 'elt') : xem[tar][rf_vars] = rf_array*0.01
        else : xem[tar][rf_vars] = rf_array

# Parse root files into list corresponding the central momentum
for tar, tar_dict in xem.items():
    # Creats sorted array of unique central momentum settings
    xem[tar]['pcent_list'] = np.unique(xem[tar]['pcent'])
    # Initialize root file list containers
    tmp_rf_list = []
    rf_list = []
    xem[tar]['chain_list'] = []
    # Enumerate condensed central momentum list
    for index, pcent_val in enumerate(xem[tar]['pcent_list']):
        # Make shallow copy of list so that when the temporary list is deleted an instance remains
        rf_list = list(tmp_rf_list)
        del tmp_rf_list[:]
        # Cleanup vacancy as a result of deleting the instance of the temporary list
        if (len(rf_list) != 0) : 
            if (len(xem[tar]['chain_list'][index-1]) == 0) : xem[tar]['chain_list'].pop(index-1)
            xem[tar]['chain_list'].append(rf_list)
        # Enumerate the full central momentum list
        for iindex, ppcent_val in enumerate(xem[tar]['pcent']):
            # If the central momenta from the two lists then fill root file containers
            if (xem[tar]['pcent_list'][index] == xem[tar]['pcent'][iindex]) :
                tmp_rf_list.append(xem[tar]['data'][iindex])
        # Populate the root file list corresponding to the respective momenta
        xem[tar]['chain_list'].append(tmp_rf_list)
