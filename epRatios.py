#!/usr/bin/python

# Import various modules
import os, sys, time
import string
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import ROOT as R

from glob import *
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')
mpl.rcParams.update({'errorbar.capsize': 2})

# Define the system clock
startTime = time.time()

# Acquire user input for spectrometer
if (sys.argv[1] != 'hms' and sys.argv[1] != 'shms') :
    print 'Usage: python getCharge.py hms (shms)'
    sys.exit(1)

# Define run list files if using text file
if (sys.argv[1] == 'hms') :
    spec = 'h'
    rlf = 'hms-xem-list.txt'
    rl = np.genfromtxt(rlf, dtype = int, skip_header=1)
    rfp = 'hms-reports/replay_hms_production_'
    rfs = '_-1.report'
    rf  = []
    dfp = 'hms-data/hms_replay_production_'
    dfs = '_-1.root'
    df  = []
if (sys.argv[1] == 'shms') :
    spec = 'p'
    rlf = 'shms-xem-list.txt'
    rl = np.genfromtxt(rlf, dtype = int, skip_header=1)
    rfp = 'shms-reports/replay_shms_production_'
    rfs = '_-1.report'
    rf  = []
    dfp = 'shms-data/shms_replay_production_'
    dfs = '_-1.root'
    df  = []

# Define run list and report files
#if (sys.argv[1] == 'hms') :
#    spec = 'h'
#    rf   = glob('hms-reports/replay_hms_production_*_-1.report')
#    df   = glob('hms-data/hms_replay_production_*_-1.root')
#if (sys.argv[1] == 'shms') : 
#    spec = 'p'
#    rf   = glob('shms-reports/replay_shms_production_*_-1.report')
#    df   = glob('shms-data/shms_replay_production_*_-1.root')

# Define constants
mp      = 0.93827231 # (GeV) mass of proton
avn     = 6.0221409e+23 # Avogadro's number
al_den  = 2.699   # density (g/cm^3) of AL7075
lh2_den = 0.07080 # density (g/cm^3) of LH2
ld2_den = 0.1638  # density (g/cm^3) of LD2 
be9_ath = 1.3140  # areal thickness (g/cm^2) of 9Be
b10_ath = 0.5722  # areal thickness (g/cm^2) of 10B4C
b11_ath = 0.6348  # areal thickness (g/cm^2) of 11B4C
c12_ath = 0.5244  # areal thickness (g/cm^2) of 1.5% C12
# Atomic masses 
lh2_am  = 1.008   # (g/mole) of LH2
ld2_am  = 2.01410 # (g/mole) of LD2
be9_am  = 9.01218 # (g/mole) of 9Be
b10_am  = 10.0129 # (g/mole) of 10B
b11_am  = 11.0093 # (g/mole) of 11B
c12_am  = 12.0107 # (g/mole) of C12
ald_am  = 26.9815 # (g/mole) of AL7075
# Define target properties
l2_entr  = 0.0150 # thickness (cm) of loop 2 (lh2) entrance window
l2_exit  = 0.0191 # thickness (cm) of loop 2 (lh2) exit window
l3_entr  = 0.0130 # thickness (cm) of loop 3 (ld2) entrance window
l3_exit  = 0.0188 # thickness (cm) of loop 3 (ld2) exit window
ald_entr = 0.1816 # areal thickness (g/cm^2) of 10 cm aluminum (AL7075) dummy entrance
ald_exit = 0.1815 # areal thickness (g/cm^2) of 10 cm aluminum (AL7075) dummy exit
# Calculate target scale factors
l2_len  = 10.0    # length (cm) of loop 2 (lh2) target cell
l3_len  = 10.0    # length (cm) of loop 3 (ld2) target cell
ald     = (ald_entr+ald_exit)/al_den # total thickness (g/cm^2)
l2_dsf  = (l2_entr+l2_exit)/ald
l3_dsf  = (l3_entr+l3_exit)/ald
num_lh2 = lh2_den*l2_len*avn/lh2_am
num_ld2 = ld2_den*l3_len*avn/ld2_am
num_be9 = be9_ath*avn/be9_am
num_b10 = b10_ath*avn/b10_am
num_b11 = b11_ath*avn/b11_am
num_c12 = c12_ath*avn/c12_am
num_ald = (ald_entr+ald_exit)*avn/ald_am

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
# Target dictionaries
#xem_tar = { 'ld2' : 2.01,
#            'c12' : 12.01,
#            'ald' : 26.98 }
# Number of nuclei in target
#xem_tar_num_nucl = { 'ld2' : num_ld2,
#                     'c12' : num_c12,
#                     'ald' : num_ald }
#xem_tar_atmc_num = { 'ld2' : 1.0,
#                     'c12' : 6.0,
#                     'ald' : 13.0 }
xem_tar_num_nucl = { 'lh2' : num_lh2,
                     'ld2' : num_ld2,
                     'be9' : num_be9,
                     'b10' : num_b10,
                     'b11' : num_b11,
                     'c12' : num_c12,
                     'ald' : num_ald }
xem_tar_atmc_num = { 'lh2' : 1.0,
                     'ld2' : 1.0,
                     'be9' : 4.0,
                     'b10' : 5.0,
                     'b11' : 5.0,
                     'c12' : 6.0,
                     'ald' : 13.0 }

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
            'scin_eff' : [],  # 3/4 trigger efficiency
            'psfactor' : [] } # el_real (ptrig2) pre-scale factor

# Store values of interest in arrays
#for index, run in enumerate(rf):
for index, run in enumerate(rl):
    rf.append(rfp + str(rl[index]) + rfs)
    df.append(dfp + str(rl[index]) + dfs)
    xem_rpf['data'].append(df[index])
    with open(rf[index]) as fobj:
        for line in fobj:
            data = line.split(':')
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
            psdata = data[0].split('=')
            if ('Ps2_factor' in psdata[0]) : xem_rpf['psfactor'].append(psdata[1].strip())

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
    # Calculate the efficiency (and pre-scale) corrected charge (electrons)
    xem[tar]['eff_corr_q4a']        = xem[tar]['etot_eff']*xem[tar]['q4a']
    xem[tar]['eff_corr_q4a_cut']    = xem[tar]['etot_eff']*xem[tar]['q4a_cut']
    xem[tar]['eff_ps_corr_q4a_cut'] = xem[tar]['etot_eff']*xem[tar]['q4a_cut'] / xem[tar]['psfactor']
    
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
                #tmp_ecq_list.append(xem[tar]['eff_corr_q4a_cut'][iindex])
                tmp_ecq_list.append(xem[tar]['eff_ps_corr_q4a_cut'][iindex])
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
        tree_chain = R.TChain('T')
        # Enumerate the list of ROOT files to be chained together
        for df_index, df_list in enumerate(xem[tar]['chain_list'][index]):
            tree_chain.Add(xem[tar]['chain_list'][index][df_index])
        # Populate the list of TChain objects for each momentum setting
        xem[tar]['tree_chain'].append(tree_chain)
        
# Create ROOT file with histograms
if (sys.argv[1] == 'hms') :  xem_rof = R.TFile('xem_hms_test.root', 'recreate')
if (sys.argv[1] == 'shms') : xem_rof = R.TFile('xem_shms_test.root', 'recreate')
#if (sys.argv[1] == 'hms') :  xem_rof = R.TFile('xem_hms_eprime_full.root', 'recreate')
#if (sys.argv[1] == 'shms') : xem_rof = R.TFile('xem_shms_eprime_full.root', 'recreate')
for tar, tar_dict in xem.items():
    # Add LaTeX format for target strings
    if (tar == 'ald') : tarStr = 'Al Dummy'
    if (tar == 'ld2') : tarStr = 'LD_{2}'
    if (tar == 'lh2') : tarStr = 'LH_{2}'
    if (tar == 'be9') : tarStr = '{}^{9}Be'
    if (tar == 'b10') : tarStr = '{}^{10}B'
    if (tar == 'b11') : tarStr = '{}^{11}B'
    if (tar == 'c12') : tarStr = '{}^{12}C'
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        xem_rof.mkdir('%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        xem_rof.cd('%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        #nentries = xem[tar]['tree_chain'][index].GetEntries() 
        nentries = 0
        # Define histograms
        hxbj            = R.TH1F('hxbj_%s_%s' % (tar, xem[tar]['pcent_list'][index]),            'x_{Bj} for %s, %s GeV; x_{Bj}; Number of Entries / 0.025' % (tarStr, xem[tar]['pcent_list'][index]), 60, 0, 1.5)
        hytar           = R.TH1F('hytar_%s_%s' % (tar, xem[tar]['pcent_list'][index]),           'y_{tar} for %s, %s GeV; y_{tar} (cm); Number of Entries / 0.1 cm' % (tarStr, xem[tar]['pcent_list'][index]), 100, -5.0, 5.0)
        heprime         = R.TH1F('heprime_%s_%s' % (tar, xem[tar]['pcent_list'][index]),         'E\' for %s, %s GeV; E\' (GeV); Number of Entries / 0.100 GeV' % (tarStr, xem[tar]['pcent_list'][index]), 120, 0.0, 12.0)
        hw2_vs_xbj      = R.TH2F('hw2_vs_xbj_%s_%s' % (tar, xem[tar]['pcent_list'][index]),      'W^{2} vs. x_{Bj} for %s, %s GeV; x_{Bj} / 0.025; W^{2} / 0.010 GeV^{2}' % (tarStr, xem[tar]['pcent_list'][index]), 60, 0, 1.5, 1500, 0, 15.0)
        hdp_vs_theta    = R.TH2F('hdp_vs_theta_%s_%s' % (tar, xem[tar]['pcent_list'][index]),    '#deltap vs. (#theta_{c}-#theta) for %s, %s GeV; #theta_{c}-#theta / 0.01 deg; #deltap / 0.5%%' % (tarStr, xem[tar]['pcent_list'][index]), 100, -5.0, 5.0, 68, -12.0, 22.0)
        hxptar_vs_yptar = R.TH2F('hxptar_vs_yptar_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'y\'_{tar} vs. x\'_{tar} for %s, %s GeV; x\'_{tar} / 1 mrad; y\'_{tar} / 1 mrad' % (tarStr, xem[tar]['pcent_list'][index]), 200, -100, 100, 200, -100, 100.0)
        # Loop over the entries in the trees
        print '\nAnalyzing the %s target at %s GeV.  There are %d events to be analyzed.\n' % (tar.upper(), xem[tar]['pcent_list'][index], nentries)
        for entry in range(nentries):
            xem[tar]['tree_chain'][index].GetEntry(entry)
            if ((entry % 100000) == 0 and entry != 0) : print 'Analyzed %d events...' % entry
            # Acquire the leaves of interest
            # PID variables
            if (sys.argv[1] == 'hms') :
                lnpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.cer.npeSum'); npeSum = lnpeSum.GetValue(0)
            if (sys.argv[1] == 'shms') :
                lhgcNpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.hgcer.npeSum'); hgcNpeSum  = lhgcNpeSum.GetValue(0)
                lngcNpeSum = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.ngcer.npeSum'); ngcNpeSum  = lngcNpeSum.GetValue(0)                
            letracknorm = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.cal.etracknorm');  etracknorm = letracknorm.GetValue(0)
            # Phase space & acceptance variables
            ldelta  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.dp'); delta  = ldelta.GetValue(0)
            lxtar   = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.x');  xtar   = lxtar.GetValue(0) 
            lytar   = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.y');  ytar   = lytar.GetValue(0) 
            lxptar  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.th'); xptar  = 1000.0*lxptar.GetValue(0) # convert to mrad
            lyptar  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.ph'); yptar  = 1000.0*lyptar.GetValue(0) # convert to mrad
            leprime = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.gtr.p');  eprime = leprime.GetValue(0) # convert to mrad
            # Kinematic variables
            lw2     = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.W2'); w2 = lw2.GetValue(0)
            lq2     = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.Q2'); q2 = lq2.GetValue(0)
            lxbj    = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.x_bj'); xbj = lxbj.GetValue(0)
            ltheta  = xem[tar]['tree_chain'][index].GetLeaf(spec.upper() + '.kin.scat_ang_deg'); theta = ltheta.GetValue(0)
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
                hgcNpeCut = bool(hgcNpeSum < 0.5)
                ngcNpeCut = bool(ngcNpeSum < 7.5)
                npeCut    = bool(hgcNpeCut or ngcNpeCut)
                deltaCut  = bool(delta < -10.0 or delta > 20.0)
                xptarCut  = bool(abs(xptar) > 70.0)
            w2Cut         = bool(w2 < 2.0) # select the DIS regime
            yptarCut      = bool(abs(yptar) > 50.0)
            etracknormCut = bool(etracknorm < 0.85)
            # Impose fiducial cuts
            if (npeCut or deltaCut or etracknormCut or w2Cut or xptarCut or yptarCut) : continue
            # Fill the histograms
            hxbj.Fill(xbj)
            heprime.Fill(eprime)
            hytar.Fill(ytar)
        # Populate efficency corrected charge histograms
        # xbj
        hxbj_qNorm = hxbj.Clone()
        hxbj_qNorm.SetNameTitle('hxbj_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized x_{Bj} for %s, %s GeV; x_{Bj} / 0.025; Y / #epsilonQ (Counts / mC)' % (tarStr, xem[tar]['pcent_list'][index]))
        hxbj_qNorm.Sumw2()
        hxbj_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # eprime
        heprime_qNorm = heprime.Clone()
        heprime_qNorm.SetNameTitle('heprime_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized E\' for %s, %s GeV; E\' / 0.100 GeV; Y / #epsilonQ (Counts / mC)' % (tarStr, xem[tar]['pcent_list'][index]))
        heprime_qNorm.Sumw2()
        heprime_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # ytar
        hytar_qNorm = hytar.Clone()
        hytar_qNorm.SetNameTitle('hytar_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized y_{tar} for %s, %s GeV; y_{tar} / 0.1 cm; Y / #epsilonQ (Counts / mC)' % (tarStr, xem[tar]['pcent_list'][index]))
        hytar_qNorm.Sumw2()
        hytar_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # delta
        hdp_qNorm = hdp_vs_theta.ProjectionY()
        hdp_qNorm.SetNameTitle('hdp_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized #deltap for %s, %s GeV; #deltap / 0.5%%; Y / #epsilonQ (Counts / mC)' % (tarStr, xem[tar]['pcent_list'][index]))
        hdp_qNorm.Sumw2()
        hdp_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # w2
        hw2_qNorm = hw2_vs_xbj.ProjectionY()
        hw2_qNorm.SetNameTitle('hw2_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]), 'Charge Normalized W^{2} for %s, %s GeV; W^{2} / 0.010 GeV^{2}; Y / #epsilonQ (Counts / mC)' % (tarStr, xem[tar]['pcent_list'][index]))
        hw2_qNorm.Sumw2()
        hw2_qNorm.Scale(1. / xem[tar]['ecq'][index])
        # Write the histograms to tape
        xem_rof.Write()
        hdp_qNorm.Delete() # address a behavior with projection I do not understand
        hw2_qNorm.Delete() # address a behavior with projection I do not understand
        xem_rof.cd('../')
# Close the ROOT file
xem_rof.Close()

print '\nThe analysis took %.3f minutes\n' % ((time.time() - startTime) / (60.))

# Open ROOT files produced above so that ratios can be calculated
if (sys.argv[1] == 'hms')  : xem_rof = R.TFile('xem_hms_eprime_full.root',  'read')
if (sys.argv[1] == 'shms') : xem_rof = R.TFile('xem_shms_eprime_full.root', 'read')
#if (sys.argv[1] == 'hms')  : xem_rof = R.TFile('xem_hms_test.root',  'read')
#if (sys.argv[1] == 'shms') : xem_rof = R.TFile('xem_shms_test.root', 'read')

# Convert histos in numpy arrays for easier manipulation
for tar, tar_dict in xem.items():
    # Create containers to store yields and bin centered values
    eprime_raw_yield_list = []
    eprime_val_list = []
    eprime_yield_list = []
    xem[tar]['eprime_raw_yield'] = []
    xem[tar]['eprime_val'] = []
    xem[tar]['eprime_yield'] = []
    xem[tar]['eprime_yield_err'] = []
    xem[tar]['eprime_yield_max'] = []
    xem[tar]['eprime_yield_min'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        # Descend into directory with histos of interest
        xem_rof.cd('%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        # Define temporary place holders for histo and array objects
        # Get raw histo and place contents in array
        tmp_raw_heprime = xem_rof.FindObjectAny('heprime_%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        tmp_raw_aeprime = tmp_raw_heprime.GetArray() # returns number of x bins +2 (over&underflow)
        tmp_raw_aeprime.SetSize(tmp_raw_heprime.GetNbinsX()) # returns number of x bins +2 (over&underflow)
        # Get charge normalized histo and place contents in array
        tmp_heprime     = xem_rof.FindObjectAny('heprime_qNorm_%s_%s' % (tar, xem[tar]['pcent_list'][index]))
        tmp_aeprime     = tmp_heprime.GetArray()
        tmp_aeprime.SetSize(tmp_heprime.GetNbinsX()) # returns number of x bins +2 (over&underflow)
        # Define bin centering arrays
        eprime_xval    = np.linspace(tmp_heprime.GetXaxis().GetXmin(), tmp_heprime.GetXaxis().GetXmax() - tmp_heprime.GetXaxis().GetBinWidth(1), num = tmp_heprime.GetNbinsX())
        eprime_offset  = tmp_heprime.GetXaxis().GetBinWidth(1) / 2.
        tmp_eprime_arr = eprime_xval + eprime_offset
        # Get maximum and conditional value of each histogram
        max_yield_val = tmp_heprime.GetMaximum()
        min_yield_val = 0.30*max_yield_val
        # Fill arrays with histo content
        tmp_raw_heprime_arr = np.array(tmp_raw_heprime)[:-2] # delete the last two over&underflow elements
        tmp_heprime_arr     = np.array(tmp_heprime)[:-2]     # delete the last two over&underflow elements
        # Fill conditional arrays if desired
        cond_raw_heprime_arr = tmp_raw_heprime_arr
        #cond_heprime_arr = tmp_heprime_arr[tmp_heprime_arr > min_yield_val]
        #cond_xval_arr = tmp_eprime_arr[tmp_heprime_arr > min_yield_val]
        cond_heprime_arr     = tmp_heprime_arr
        cond_xval_arr     = tmp_eprime_arr
        # Store yields and bin centers in lists for each momentum
        eprime_raw_yield_list.append(cond_raw_heprime_arr)
        eprime_val_list.append(cond_xval_arr)
        eprime_yield_list.append(cond_heprime_arr)
        xem[tar]['eprime_yield_max'].append(max_yield_val)
        xem[tar]['eprime_yield_min'].append(min_yield_val)
    # Store yields and bin center lists in xem dictionary indexed on xem['pcent_list']
    xem[tar]['eprime_val']       = eprime_val_list
    xem[tar]['eprime_raw_yield'] = eprime_raw_yield_list
    xem[tar]['eprime_yield']     = eprime_yield_list
    xem[tar]['eprime_yield_err'] = np.sqrt(eprime_raw_yield_list)*(1. / xem[tar]['ecq'][index])

# Calculate the dummy corrected yields for cryo targets
for tar, tar_dict in xem.items():
    # Create containers to store dummy corrected yields
    eprime_dc_yield_list = []
    xem[tar]['eprime_dc_yield'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        if (tar == 'lh2') :
            eprime_dc_yield_list.append(xem['lh2']['eprime_yield'][index] - xem['ald']['eprime_yield'][index]*l2_dsf)
        if (tar == 'ld2') :
            eprime_dc_yield_list.append(xem['ld2']['eprime_yield'][index] - xem['ald']['eprime_yield'][index]*l3_dsf)
        else :
            eprime_dc_yield_list.append(xem[tar]['eprime_yield'][index])
    xem[tar]['eprime_dc_yield'] = eprime_dc_yield_list
    xem[tar]['eprime_dc_yield']
# Calculate the ratios of yields
for tar, tar_dict in xem.items():
    # Create containers to store ratios and error on ratios
    eprime_ratio_list = []
    xem[tar]['eprime_ratio'] = []
    eprime_ratio_err_list = []
    xem[tar]['eprime_ratio_err'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        eprime_ratio_list.append(np.divide(xem[tar]['eprime_dc_yield'][index]*(1./(xem_tar_num_nucl[tar]*xem_tar_atmc_num[tar])), 
                                        xem['ld2']['eprime_dc_yield'][index]*(1./(xem_tar_num_nucl['ld2']*xem_tar_atmc_num['ld2'])), 
                                        where = xem['ld2']['eprime_dc_yield'][index] > 0.0))
    xem[tar]['eprime_ratio'] = eprime_ratio_list
    # Calculate error on ratios
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        eprime_ratio_err_list.append(xem[tar]['eprime_ratio'][index]*np.sqrt(np.divide(xem[tar]['eprime_yield_err'][index],   xem[tar]['eprime_yield'][index],   where = xem[tar]['eprime_yield'][index] > 0.0)**2.0 + 
                                                                       np.divide(xem['ld2']['eprime_yield_err'][index], xem['ld2']['eprime_yield'][index], where = xem['ld2']['eprime_yield'][index] > 0.0)**2.0))
    xem[tar]['eprime_ratio_err'] = eprime_ratio_err_list

# Truncate all non-zero ratios for plotting
for tar, tar_dict in xem.items():
    # Initialize lists
    eprime_nz_val_list = []
    eprime_nz_ratio_list = []
    eprime_nz_ratio_err_list = []
    xem[tar]['eprime_nz_val'] = []
    xem[tar]['eprime_nz_ratio'] = []
    xem[tar]['eprime_nz_ratio_err'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        # Truncate all non-zero ratios for plotting
        eprime_nz_val_list.append(xem[tar]['eprime_val'][index][xem[tar]['eprime_ratio'][index]>1.0e-6])
        eprime_nz_ratio_list.append(xem[tar]['eprime_ratio'][index][xem[tar]['eprime_ratio'][index]>1.0e-6])
        eprime_nz_ratio_err_list.append(xem[tar]['eprime_ratio_err'][index][xem[tar]['eprime_ratio'][index]>1.0e-6])
    xem[tar]['eprime_nz_val'] = eprime_nz_val_list
    xem[tar]['eprime_nz_ratio'] = eprime_nz_ratio_list
    xem[tar]['eprime_nz_ratio_err'] = eprime_nz_ratio_err_list

# Plot the ratios in bins of E'
hmkr = ['bo', 'g^', 'rs', 'kd', 'm*']
pmkr = ['bo', 'g^', 'rs', 'kd']
for tar, tar_dict in xem.items():
    # Add LaTeX format for target strings
    if (tar == 'c12') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['eprime_nz_val'][index], xem[tar]['eprime_nz_ratio'][index], yerr = xem[tar]['eprime_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') : 
                plt.errorbar(xem[tar]['eprime_nz_val'][index], xem[tar]['eprime_nz_ratio'][index], yerr = xem[tar]['eprime_nz_ratio_err'][index],
                             fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(2.0, 6.5)
        plt.ylim(0.85, 1.2)
        plt.legend(loc = 2)
        plt.show()
        plt.clf()

# Define function to calculate xbj from bins in E'
def calc_xbj(ep, eb, theta) :
    return ((eb*ep*(1.0 - np.cos(np.deg2rad(theta))))/(mp*(eb - ep)))

# Calculate xbj from E' and store in dictionary
for tar, tar_dict in xem.items():
    xbj_calc_nz_val_list = []
    xbj_calc_nz_ratio_list = []
    xem[tar]['xbj_calc_nz_val'] = []
    xem[tar]['xbj_calc_nz_ratio'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        xbj_calc_nz_val_list.append(calc_xbj(xem[tar]['eprime_nz_val'][index], xem[tar]['ebeam'][index], xem[tar]['theta'][index]))
        xbj_calc_nz_ratio_list.append(xem[tar]['eprime_nz_ratio'][index])
    xem[tar]['xbj_calc_nz_val'] = xbj_calc_nz_val_list
    xem[tar]['xbj_calc_nz_ratio'] = xbj_calc_nz_ratio_list

# Plot the ratios in bins of xbj from bins of E'
for tar, tar_dict in xem.items():
    if (tar == 'c12') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio'][index], yerr = xem[tar]['eprime_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio'][index], yerr = xem[tar]['eprime_nz_ratio_err'][index],
                             fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(0.0, 1.15)
        plt.ylim(0.85, 1.2)
        plt.legend(loc = 2)
        plt.show()
        plt.clf()
    
# Import the radiative corrections table and parse the columns into arrays
ebeam, eprime, theta, xbj, q2, w2, rcd, rcc, ratio, cc = np.loadtxt('rc_f1f2_ineft.txt', skiprows = 1, unpack=True)
# Define xbj range to inerpolate over and interpolate the data
xr     = np.linspace(0.2, 1.25, 10000)
interp = np.interp(xr, xbj, ratio)
# Define function to return the same indexed value of the interpolated data
def get_rcf(x):
    index = (np.abs(xr - x)).argmin()
    return interp[index]
# Create array with same bin centering as histogram defined above
for tar, tar_dict in xem.items():
    xbj_full_list = []
    rcf_list = []
    tmp_rcf_arr = []
    tmp_cp_rcf_arr = []
    xem[tar]['rcf_list'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        xbj_full_list.append(calc_xbj(xem[tar]['eprime_val'][index], xem[tar]['ebeam'][index], xem[tar]['theta'][index]))
        tmp_cp_rcf_arr = list(tmp_rcf_arr)
        del tmp_rcf_arr[:]
        for each_xbj in np.nditer(xbj_full_list[index]):
            tmp_rcf_arr.append(get_rcf(each_xbj))
        rcf_list.append(np.asarray(tmp_rcf_arr))
    xem[tar]['rcf_list'] = rcf_list
# Plot the data
fig, (ax0, ax1) = plt.subplots(nrows = 2, sharex = True, sharey=True)
ax0.plot(xbj, ratio, 'bo', markersize = 7)
ax0.plot(xr, interp, 'r-', linewidth = 2)
ax1.plot(xbj_full_list[0], rcf_list[0], 'dk', markersize = 7)
plt.xlim(0.2, 1.15)
plt.ylim(0.875, 1.05)
plt.show()

# Apply the radiative corrections to the data
for tar, tar_dict in xem.items():
    xbj_calc_nz_val_list = []
    xbj_calc_nz_ratio_list = []
    xbj_calc_nz_ratio_err_list = []
    xbj_nz_rcf_list = []
    xbj_calc_nz_ratio_corr_list = []
    xem[tar]['xbj_calc_nz_val'] = []
    xem[tar]['xbj_calc_nz_ratio'] = []
    xem[tar]['xbj_calc_nz_ratio_err'] = []
    xem[tar]['xbj_nz_rcf'] = []
    xem[tar]['xbj_calc_nz_ratio_corr'] = []
    for index, mom_list in enumerate(xem[tar]['pcent_list']):
        xbj_calc_nz_val_list.append(calc_xbj(xem[tar]['eprime_nz_val'][index], xem[tar]['ebeam'][index], xem[tar]['theta'][index]))
        xbj_calc_nz_ratio_list.append(calc_xbj(xem[tar]['eprime_nz_val'][index], xem[tar]['ebeam'][index], xem[tar]['theta'][index]))
        xbj_calc_nz_ratio_err_list.append(calc_xbj(xem[tar]['eprime_nz_val'][index], xem[tar]['ebeam'][index], xem[tar]['theta'][index])*np.sqrt(2.)*
                                          np.divide(xem[tar]['eprime_nz_ratio_err'][index], xem[tar]['eprime_nz_ratio'][index]))
        xbj_nz_rcf_list.append(xem[tar]['rcf_list'][index][xem[tar]['eprime_ratio'][index]>1.0e-6])
        xbj_calc_nz_ratio_corr_list.append(xem[tar]['eprime_nz_ratio'][index]*xem[tar]['rcf_list'][index][xem[tar]['eprime_ratio'][index]>1.0e-6])
    xem[tar]['xbj_calc_nz_val'] = xbj_calc_nz_val_list
    xem[tar]['xbj_calc_nz_ratio'] = xbj_calc_nz_ratio_list
    xem[tar]['xbj_calc_nz_ratio_err'] = xbj_calc_nz_ratio_err_list
    xem[tar]['xbj_nz_rcf'] = xbj_nz_rcf_list
    xem[tar]['xbj_calc_nz_ratio_corr'] = xbj_calc_nz_ratio_corr_list

# Plot the radiative corrected ratios in bins of xbj from bins of E'
for tar, tar_dict in xem.items():
    if (tar == 'c12') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio_corr'][index], yerr = xem[tar]['xbj_calc_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio_corr'][index], yerr = xem[tar]['xbj_calc_nz_ratio_err'][index],
                             fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(0.0, 1.15)
        plt.ylim(0.85, 1.2)
        plt.legend(loc = 2)
        plt.show()
        plt.clf()


# Plot the radiative corrected ratios in bins of xbj from bins of E'
for tar, tar_dict in xem.items():
    if (tar == 'c12') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio_corr'][index], yerr = xem[tar]['xbj_calc_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') : 
                if (xem[tar]['pcent_list'][index] == 3.3 or xem[tar]['pcent_list'][index] == 4.0) :
                    tv  = np.delete(xem[tar]['xbj_calc_nz_val'][index], -1)
                    tr  = np.delete(xem[tar]['xbj_calc_nz_ratio_corr'][index], -1)
                    tre = np.delete(xem[tar]['xbj_calc_nz_ratio_err'][index], -1)
                    plt.errorbar(tv, np.zeros(np.size(tr))+1.0, yerr = tre,
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
                else :
                    plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], np.zeros(np.size(xem[tar]['xbj_calc_nz_ratio_corr'][index]))+1.0, yerr = xem[tar]['xbj_calc_nz_ratio_err'][index],
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(0.0, 1.025)
        plt.ylim(0.9, 1.1)
        plt.legend(loc = 2)
        plt.xlabel(r'$x_{Bj}$', size = 20)
        plt.ylabel(r'$\sigma_{{}^{12}C} / \sigma_{D}$', size = 20)
        plt.show()
        plt.clf()
        #plt.savefig('c12_emc_ratio_stat_errors.png', format='png', figsize = (15, 15), dpi=1000)
    if (tar == 'b10') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio_corr'][index], yerr = xem[tar]['xbj_calc_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') :
                if (xem[tar]['pcent_list'][index] == 3.3 or xem[tar]['pcent_list'][index] == 4.0) :
                    tv  = np.delete(xem[tar]['xbj_calc_nz_val'][index], -1)
                    tr  = np.delete(xem[tar]['xbj_calc_nz_ratio_corr'][index], -1)
                    tre = np.delete(xem[tar]['xbj_calc_nz_ratio_err'][index], -1)
                    plt.errorbar(tv, np.zeros(np.size(tr))+1.0, yerr = tre,
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
                else :
                    plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], np.zeros(np.size(xem[tar]['xbj_calc_nz_ratio_corr'][index]))+1.0, yerr = xem[tar]['xbj_calc_nz_ratio_err'][index],
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(0.0, 1.025)
        plt.ylim(0.9, 1.1)
        plt.legend(loc = 2)
        plt.xlabel(r'$x_{Bj}$', size = 20)
        plt.ylabel(r'$\sigma_{{}^{10}B} / \sigma_{D}$', size = 20)
        plt.show()
        plt.clf()
        #plt.savefig('b10_emc_ratio_stat_errors.png', format='png', figsize = (15, 15), dpi=1000)
    if (tar == 'b11') :
        for index, mom_list in enumerate(xem[tar]['pcent_list']):
            # Truncate all non-zero ratios for plotting     
            if (sys.argv[1] == 'hms' and xem[tar]['pcent_list'][index] != 3.3) : 
                plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], xem[tar]['xbj_calc_nz_ratio_corr'][index], yerr = xem[tar]['xbj_calc_nz_ratio_err'][index], 
                             fmt = '%s' % hmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
            elif (sys.argv[1] == 'shms') :
                if (xem[tar]['pcent_list'][index] == 3.3 or xem[tar]['pcent_list'][index] == 4.0) :
                    tv  = np.delete(xem[tar]['xbj_calc_nz_val'][index], -1)
                    tr  = np.delete(xem[tar]['xbj_calc_nz_ratio_corr'][index], -1)
                    tre = np.delete(xem[tar]['xbj_calc_nz_ratio_err'][index], -1)
                    plt.errorbar(tv, np.zeros(np.size(tr))+1.0, yerr = tre,
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
                else :
                    plt.errorbar(xem[tar]['xbj_calc_nz_val'][index], np.zeros(np.size(xem[tar]['xbj_calc_nz_ratio_corr'][index]))+1.0, yerr = xem[tar]['xbj_calc_nz_ratio_err'][index],
                                 fmt = '%s' % pmkr[index], label = '%s GeV' % xem[tar]['pcent_list'][index], markersize=7, alpha=0.75)
        plt.xlim(0.0, 1.025)
        plt.ylim(0.9, 1.1)
        plt.legend(loc = 2)
        plt.xlabel(r'$x_{Bj}$', size = 20)
        plt.ylabel(r'$\sigma_{{}^{11}B} / \sigma_{D}$', size = 20)
        plt.show()
        plt.clf()
        #plt.savefig('b11_emc_ratio_stat_errors.png', format='png', figsize = (15, 15), dpi=1000)
