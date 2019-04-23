#!/usr/bin/python3.6

import glob, re, json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib import rc

rc('text', usetex = True)
rc('font', family = 'serif')
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'errorbar.capsize': 6})

numFecs = 5; numLinks = 2

dfStats = glob.glob('data/run*/*stats.json')
dfAdc   = glob.glob('data/run*/*adc.txt')

dfStats.sort(); dfAdc.sort()

dataDict = {}

for fecId in range(numFecs) :
    dataDict['fec%d' % fecId] = {}
    for linkId in range(numLinks) :
        dataDict['fec%d' % fecId]['link%d' % linkId] = {}
        dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles'] = []
        dataDict['fec%d' % fecId]['link%d' % linkId]['adcFiles']  = []
        dataDict['fec%d' % fecId]['link%d' % linkId]['statData']  = []
        dataDict['fec%d' % fecId]['link%d' % linkId]['adcData']   = []
        dataDict['fec%d' % fecId]['link%d' % linkId]['chanNums']  = []
        for index, run in enumerate(dfStats) :
            if (fecId+1 == int(re.findall(r'\d+', dfStats[index])[0][3:4]) and 
                linkId  == int(re.findall(r'\d+', dfStats[index])[3][1:2])) :
                dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles'].append(dfStats[index])
                dataDict['fec%d' % fecId]['link%d' % linkId]['adcFiles'].append(dfAdc[index])
        for index, run in enumerate(dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles']) :   
            with open(dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles'][index]) as jdf :
                dataDict['fec%d' % fecId]['link%d' % linkId]['statData'].append(json.load(jdf))
        for index, run in enumerate(dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles']) :
            if   (index == 0) :
                dataDict['fec%d' % fecId]['link%d' % linkId]['chanNums'].append(np.arange(1, len(np.asarray(dataDict['fec%d' % fecId]['link%d' % linkId]['statData'][index]['mean']))+1, 1))
            elif (index == 1) :
                dataDict['fec%d' % fecId]['link%d' % linkId]['chanNums'].append(np.arange(81, len(np.asarray(dataDict['fec%d' % fecId]['link%d' % linkId]['statData'][index]['mean']))+81, 1))
      
for fecId in range(numFecs) :
    plt.figure(fecId)
    for linkId in range(numLinks) :
        for index, run in enumerate(dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles']) :
            runNum = int(re.findall(r'\d+', dataDict['fec%d' % fecId]['link%d' % linkId]['statFiles'][index])[0])
            plt.subplot(1, 2, index+1)
            dpc = np.where(dataDict['fec%d' % fecId]['link%d' % linkId]['chanNums'][linkId]%2 == 0, 'b', 'r')
            plt.scatter(dataDict['fec%d' % fecId]['link%d' % linkId]['chanNums'][linkId], 
                     dataDict['fec%d' % fecId]['link%d' % linkId]['statData'][index]['mean'], 
                     marker = 'o', color = dpc)
            plt.title('FEC %d, Run %d' % (fecId, runNum))
            plt.xlabel('Channel Number')
            plt.ylabel('$\mathrm{\sigma_{pedestal}\ (channels)}$') 
            plt.tight_layout()
    plt.show()
    plt.savefig('plots/fec%d.pdf' % fecId)
