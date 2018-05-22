#!/usr/bin/python

# Import various modules
import os, sys, time
import string
import numpy as np
import matplotlib.pyplot as plt
import ROOT as r

from glob import *

# Define the system clock
startTime = time.time()

# Acquire user input for spectrometer
if (sys.argv[1] != 'hms' and sys.argv[1] != 'shms'):
    print 'Usage: python getCharge.py hms (or shms)'
    sys.exit(1)

# Define run list files if using text file
#if (sys.argv[1] == 'hms') : 
#    rlf = 'hms-xem-list.txt'
#    rl = np.genfromtxt(hrlf, dtype = int, skip_header=1)
#    rfp = 'hms-reports/replay_hms_production_'
#    rfs = '_-1.report'
#    rf  = []
#    dfp = 'hms-data/hms_replay_production_'
#    dfs = '_-1.root'
#    df  = []
#if (sys.argv[1] == 'shms') : 
#    rlf = 'shms-xem-list.txt'
#    rl = np.genfromtxt(hrlf, dtype = int, skip_header=1)
#    rfp = 'shms-reports/replay_hms_production_'
#    rfs = '_-1.report'
#    rf  = []
#    dfp = 'shms-data/shms_replay_production_'
#    dfs = '_-1.root'
#    df  = []

# Define run list and report files
if (sys.argv[1] == 'hms') :
    spec = 'h'
    rf   = glob('hms-reports/replay_hms_production_*_-1.report')
    df   = glob('hms-data/hms_replay_production_*_-1.root')
if (sys.argv[1] == 'shms') : 
    spec = 'p'
    rf   = glob('shms-reports/replay_shms_production_*_-1.report')
    df   = glob('shms-data/shms_replay_production_*_-1.root')

# Define constants
al_den = 2.699     # density (g/cm^3) of AL7075
# Define target properties
l2_entr  = 0.0150 # thickness (cm) of loop 2 (lh2) entrance window
l2_exit  = 0.0205 # thickness (cm) of loop 2 (lh2) exit window
l3_entr  = 0.0130 # thickness (cm) of loop 3 (ld2) entrance window
l3_exit  = 0.0186 # thickness (cm) of loop 3 (ld2) exit window
ald_entr = 0.1816 # thickness (g/cm^2) of 10 cm aluminum (AL7075) dummy entrance
ald_exit = 0.1815 # thickness (g/cm^2) of 10 cm aluminum (AL7075) dummy exit
# Calculate target scale factors
l2w = l2_entr+l2_exit
l3w = l3_entr+l3_exit
ald = (ald_entr+ald_exit)/al_den
l2_sf = l2w/ald
l3_sf = l3w/ald

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
xem_rpf = { 'data'     : [],  # data file
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
            'scin_eff' : [] }  # 3/4 trigger efficiency

# Store values of interest in arrays
for index, run in enumerate(rf):
    xem_rpf['data'].append(df[index])
    with open(rf[index]) as fobj:
        for line in fobj:
            data = line.split(':')
            #if ('' in data[0]) : xem_rpf[''].append(data[1].strip())
            # Kinematic configurations
            if ('Run Num'     in data[0]) : xem_rpf['rn'].append(data[1].strip())
            if ('Momentum'    in data[0]) : xem_rpf['pcent'].append(data[1].strip())
            if ('Target AMU'  in data[0]) : xem_rpf['tamu'].append(data[1].strip())
            if ('Spec Theta'  in data[0]) : xem_rpf['theta'].append(data[1].strip())
            if ('Beam Energy' in data[0]) : xem_rpf['ebeam'].append(data[1].strip())
            # Charge and current
            # if ('' in data[0]) : xem_rpf[''].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Current' in data[0])          : xem_rpf['i4a'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Beam Cut Current' in data[0]) : xem_rpf['i4a_cut'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Charge' in data[0])           : xem_rpf['q4a'].append(filter(lambda x: x in string.digits + '.', data[1]))
            if ('BCM4a Beam Cut Charge' in data[0])  : xem_rpf['q4a_cut'].append(filter(lambda x: x in string.digits + '.', data[1]))
            # Live times (must be multiplied by 0.01 -> done later)
            if (sys.argv[1] == 'hms')  : 
                if ('Pre-Scaled Ps2 HMS Computer Live Time' in data[0])  : xem_rpf['clt'].append(data[1][:8].strip())
            if (sys.argv[1] == 'shms') : 
                if ('Pre-Scaled Ps2 SHMS Computer Live Time' in data[0]) : xem_rpf['clt'].append(data[1][:8].strip())
            if ('OG 6 GeV Electronic Live Time (100, 150)' in data[0])   : xem_rpf['elt'].append(data[1][:8].strip())
            # Tracking efficiencies
            if (sys.argv[1] == 'hms')  : 
                if ('SING FID TRACK EFFIC' in data[0])        : xem_rpf['tr_eff'].append(data[1][:11].strip())
                if ('E SING FID TRACK EFFIC' in data[0])      : xem_rpf['etr_eff'].append(data[1][:11].strip())
                if ('HADRON SING FID TRACK EFFIC' in data[0]) : xem_rpf['htr_eff'].append(data[1][:11].strip())
            if (sys.argv[1] == 'shms')  : 
                if ('SING FID TRACK EFFIC' in data[0])        : xem_rpf['tr_eff'].append(data[1][:8].strip())
                if ('E SING FID TRACK EFFIC' in data[0])      : xem_rpf['etr_eff'].append(data[1][:8].strip())
                if ('HADRON SING FID TRACK EFFIC' in data[0]) : xem_rpf['htr_eff'].append(data[1][:8].strip())
            # Trigger efficiency
            if ('3_of_4 EFF' in data[0]) : xem_rpf['scin_eff'].append(data[1].strip())

# Enumerate xem targets
for tar_str, tar_amu in xem_tar.items():
    # Initialize dictionary
    xem[tar_str] = {}
    # Enumerate variables
    for var_str, var in xem_rpf.items():
        # Initialize lists in xem dictionary
        xem[tar_str][var_str] = []
        # Enumerate target list from report files
        for index, target in enumerate(xem_rpf['tamu']):
            # Append lists when enumerated targets are identical
            if (float(xem_tar[tar_str]) == float(target)):
                xem[tar_str][var_str].append(xem_rpf[var_str][index])

# Convert lists to arrays, kludgy as hell but oh well
for tar, tar_dict in xem.items():
    for rpf_vars, rpf_list in xem[tar].items():
        if (rpf_vars == 'data') : continue
        rpf_array = np.asarray(rpf_list, dtype = float)
        del xem[tar][rpf_vars]
        if   (rpf_vars == 'clt') : xem[tar][rpf_vars] = rpf_array*0.01
        elif (rpf_vars == 'elt') : xem[tar][rpf_vars] = rpf_array*0.01
        else : xem[tar][rpf_vars] = rpf_array
    # Calculate the per run efficiency
    xem[tar]['tot_eff']  = xem[tar]['tr_eff']*xem[tar]['scin_eff']*xem[tar]['clt']*xem[tar]['elt']
    xem[tar]['etot_eff'] = xem[tar]['etr_eff']*xem[tar]['scin_eff']*xem[tar]['clt']*xem[tar]['elt']
    xem[tar]['htot_eff'] = xem[tar]['htr_eff']*xem[tar]['scin_eff']*xem[tar]['clt']*xem[tar]['elt']
    # Calculate the efficiency corrected charge (electrons)
    xem[tar]['eff_corr_q4a']     = xem[tar]['etot_eff']*xem[tar]['q4a']
    xem[tar]['eff_corr_q4a_cut'] = xem[tar]['etot_eff']*xem[tar]['q4a_cut']

# Parse root files into list corresponding the central momentum
for tar, tar_dict in xem.items():
    # Sorted array of unique central momentum settings
    xem[tar]['pcent_list'] = np.unique(xem[tar]['pcent'])
    # Initialize root file list containers
    rof_list     = []
    tmp_rof_list = []
    xem[tar]['chain_list'] = []
    # Initialize the efficiency corrected charge containers
    ecq_list     = []
    tmp_ecq_list = []
    xem[tar]['ecq_list'] = []
    # Enumerate condensed central momentum list
    for index, pcent_val in enumerate(xem[tar]['pcent_list']):
        # Make shallow copy of list so that when the temporary list is deleted an instance remains
        rof_list = list(tmp_rof_list)
        del tmp_rof_list[:]
        ecq_list = list(tmp_ecq_list)
        del tmp_ecq_list[:]
        # Cleanup vacancy as a result of deleting the instance of the temporary list
        if (len(rof_list) != 0) : 
            if (len(xem[tar]['chain_list'][index-1]) == 0) : xem[tar]['chain_list'].pop(index-1)
            xem[tar]['chain_list'].append(rof_list)
        # Enumerate the full central momentum list
        for iindex, ppcent_val in enumerate(xem[tar]['pcent']):
            # If the central momenta from the two lists then fill root file containers
            if (xem[tar]['pcent_list'][index] == xem[tar]['pcent'][iindex]) :
                tmp_rof_list.append(xem[tar]['data'][iindex])
                tmp_ecq_list.append(xem[tar]['eff_corr_q4a_cut'][iindex])
        # Populate the root file list corresponding to the respective momenta
        xem[tar]['chain_list'].append(tmp_rof_list)
        xem[tar]['ecq_list'].append(np.asarray(tmp_ecq_list))
    # Calculate the efficienct corrected charge for each target and momentum setting
    xem[tar]['ecq'] = []
    for index, pcent_val in enumerate(xem[tar]['pcent_list']):
        xem[tar]['ecq'].append(np.sum(xem[tar]['ecq_list'][index]))
            
# Chain ROOT files together per momentum setting
for tar, tar_dict in xem.items():
    # Initialize the tree chain lists
    xem[tar]['tree_chain'] = []
    tree_chain = []
    # Enumerate the individual momentum settings
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        # Initialze the TChain object for each momentum setting
        tree_chain = r.TChain('T')
        # Enumerate the list of ROOT files to be chained together
        for df_index, df_list in enumerate(xem[tar]['chain_list'][index]):
            tree_chain.Add(xem[tar]['chain_list'][index][df_index])
        # Populate the list of TChain objects for each momentum setting
        xem[tar]['tree_chain'].append(tree_chain)
        
# Create ROOT file with histograms
if (sys.argv[1] == 'hms') :  xem_rof = r.TFile('xem_hms.root', 'recreate')
if (sys.argv[1] == 'shms') : xem_rof = r.TFile('xem_shms.root', 'recreate')
for tar, tar_dict in xem.items():
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        xem_rof.mkdir('%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        xem_rof.cd('%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        #nentries = xem[tar]['tree_chain'][index].GetEntries() 
        nentries = 10000
        # Define histograms
        hxbj            = r.TH1F('hxbj_%s_%s' % (tar, xem[tar]['pcent_list'][index]),            'x_{Bj} for %s, %s GeV; x_{Bj}; Number of Entries / 0.025' % (tar.upper(), xem[tar]['pcent_list'][index]), 60, 0, 1.5)
        hytar           = r.TH1F('hytar_%s_%s' % (tar, xem[tar]['pcent_list'][index]),           'y_{tar} for %s, %s GeV; y_{tar}; Number of Entries / 0.1 cm' % (tar.upper(), xem[tar]['pcent_list'][index]), 100, -5.0, 5.0)
        hw2_vs_xbj      = r.TH2F('hw2_vs_xbj_%s_%s' % (tar, xem[tar]['pcent_list'][index]),      'W^{2} vs. x_{Bj} for %s, %s GeV; x_{Bj} / 0.025; W^{2} / 0.010 GeV^{2}' % (tar.upper(), xem[tar]['pcent_list'][index]), 60, 0, 1.5, 1500, 0, 15.0)
        hdp_vs_theta    = r.TH2F('hdp_vs_theta_%s_%s' % (tar, xem[tar]['pcent_list'][index]),    '#deltap vs. (#theta_{c}-#theta) for %s, %s GeV; #theta_{c}-#theta / 0.01 deg; #deltap / 0.1' % (tar.upper(), xem[tar]['pcent_list'][index]), 100, -5.0, 5.0, 340, -12.0, 22.0)
        hxptar_vs_yptar = r.TH2F('hxptar_vs_yptar_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'y\'_{tar} vs. x\'_{tar} for %s, %s GeV; x\'_{tar} / 1 mrad; y\'_{tar} / 1 mrad' % (tar.upper(), xem[tar]['pcent_list'][index]), 200, -100, 100, 200, -100, 100.0)
        # Loop over the entries in the trees
        print '\nAnalyzing the %s target at %s GeV.  There are %d events to be analyzed.\n' % (tar.upper(), xem[tar]['pcent_list'][index], nentries)
        for entry in range(nentries+1):
            xem[tar]['tree_chain'][index].GetEntry(entry)
            if ((entry % 100000) == 0 and entry != 0) : print 'Analyzed %d events...' % entry
            # Acquire the leaves of interest
            # = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.')
            # PID variables
            if (sys.argv[1] == 'hms') :
                lnpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.cer.npeSum'); npeSum = lnpeSum.GetValue(0)
            if (sys.argv[1] == 'shms') :
                lhgcNpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.hgcer.npeSum'); hgcNpeSum  = lhgcNpeSum.GetValue(0)
                lngcNpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.ngcer.npeSum'); ngcNpeSum  = lngcNpeSum.GetValue(0)                
            letracknorm = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.cal.etracknorm');  etracknorm = letracknorm.GetValue(0)
            # Phase space variables
            ldelta = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.dp'); delta  = ldelta.GetValue(0)
            lxtar  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.x');  xtar  = lxtar.GetValue(0) 
            lytar  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.y');  ytar  = lytar.GetValue(0) 
            lxptar = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.th'); xptar = 1000.0*lxptar.GetValue(0) # convert to mrad
            lyptar = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.ph'); yptar = 1000.0*lyptar.GetValue(0) # convert to mrad
            # Kinematic variables
            lw2    = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.W2'); w2 = lw2.GetValue(0)
            lq2    = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.Q2'); q2 = lq2.GetValue(0)
            lxbj   = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.x_bj'); xbj = lxbj.GetValue(0)
            ltheta = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.scat_ang_deg'); theta = ltheta.GetValue(0)
            # Fill histograms prior to fiducial cuts
            hw2_vs_xbj.Fill(xbj, w2)
            if (sys.argv[1] == 'hms') :  hdp_vs_theta.Fill(xem[tar]['theta'][index] + theta, delta)
            if (sys.argv[1] == 'shms') : hdp_vs_theta.Fill(xem[tar]['theta'][index] - theta, delta)
            hxptar_vs_yptar.Fill(xptar, yptar)
            # Define the fiducial cuts
            if (sys.argv[1] == 'hms') :
                npeCut   = bool(npeSum < 1.5)
                deltaCut = bool(abs(delta) > 9.0)
                xptarCut = bool(abs(xptar) > 90.0)
            if (sys.argv[1] == 'shms') :
                hgcNpeCut = bool(hgcNpeSum < 1.5)
                ngcNpeCut = bool(ngcNpeSum < 1.5)
                npeCut    = bool(hgcNpeCut or ngcNpeCut)
                deltaCut  = bool(delta < -10.0 or delta > 20.0)
                xptarCut  = bool(abs(xptar) > 70.0)
            w2Cut         = bool(w2 < 2.0) # select the DIS regime
            yptarCut      = bool(abs(yptar) > 50.0)
            etracknormCut = bool(etracknorm < 0.7)
            # Impose fiducial cuts
            if (npeCut or deltaCut or etracknormCut or w2Cut or xptarCut or yptarCut) : continue
            # Fill the histograms
            hxbj.Fill(xbj)
            hytar.Fill(ytar)
        # Populate efficency corrected charge histograms
        # xbj
        hxbj_qNorm = hxbj.Clone()
        hxbj_qNorm.SetNameTitle('hxbj_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized x_{Bj} for %s, %s GeV; x_{Bj} / 0.025; Y / #epsilonQ (Counts / mC)' % (tar.upper(), xem[tar]['pcent_list'][index]))
        hxbj_qNorm.Sumw2()
        hxbj_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # ytar
        hytar_qNorm = hytar.Clone()
        hytar_qNorm.SetNameTitle('hytar_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized y_{tar} for %s, %s GeV; y_{tar} / 0.1 cm; Y / #epsilonQ (Counts / mC)' % (tar.upper(), xem[tar]['pcent_list'][index]))
        hytar_qNorm.Sumw2()
        hytar_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # Write the histograms to tape
        xem_rof.Write()
        xem_rof.cd('../')
# Close the ROOT file
xem_rof.Close()

print '\nThe analysis took %.3f minutes\n' % ((time.time() - startTime) / (60.))
