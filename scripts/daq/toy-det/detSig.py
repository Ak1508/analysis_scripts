#!/usr/bin/python3.6

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import gumbel_r

class dataStream :

    def __init__(self, sampleRate, numEvents) :
        # sampling rate (s^-1), number of events to simulate
        self.sampleRate, self.numEvents = sampleRate, numEvents
        # number of samples in read-out window, up-to 1024 samples for SAMPA chips
        self.numSamples = 1024
        # read-out window length (s)
        self.windowWidth = (1.0 / self.sampleRate)*self.numSamples
        # length of run (s)
        self.runTime = self.numEvents * self.windowWidth
        # inital time sample
        self.initTime = 0.0

    # method to simulate time samples
    def simTimeSamples(self, windowStart, windowEnd) :
        self.timeSamples   = np.linspace(windowStart, windowEnd, int (self.numSamples))
        self.initTime      = self.timeSamples[-1]
        self.randTimeIndex = np.random.randint(int (len(self.timeSamples)*0.25), 
                                               high = int (len(self.timeSamples)*0.75))
    # method to simulate adc sample
    # def simAdcSignals(self, numSamples, randTimeIndex) :
    def simAdcSignals(self) :
        self.baseLine   = np.random.rand(int (self.numSamples))*1.0e-2
        self.hitJudge   = np.random.random()
        if (self.hitJudge < 0.5) :    # simulate no hit
            self.adcSamples = self.baseLine
        elif (self.hitJudge >= 0.5) : # simulate hit
            self.gumbelMean = np.random.random()
            self.gumbelBeta = np.random.random()
            self.gumbelPpf  = np.linspace(gumbel_r.ppf(0.001, loc = self.gumbelMean, scale = self.gumbelBeta), 
                                      gumbel_r.ppf(0.999, loc = self.gumbelMean, scale = self.gumbelBeta), 100)
            self.gumbelPdf  = gumbel_r.pdf(self.gumbelPpf, loc = self.gumbelMean, scale = self.gumbelBeta)
            self.adcSamples = np.insert(self.baseLine, self.randTimeIndex, self.gumbelPdf)[:int (self.numSamples)]

    def __iter__(self) :
        return self

    def __next__(self) :
        if (self.initTime == 0.0) :
            # simulate time and adc samples
            self.simTimeSamples(self.initTime, self.windowWidth)
            self.simAdcSignals()
        elif (self.initTime > 0.0 and self.initTime < self.runTime) :
            # simulate the time and adc signals
            self.simTimeSamples(self.initTime, self.initTime + self.windowWidth)
            self.simAdcSignals()
        elif (self.initTime >= self.runTime) : 
            raise StopIteration

sampleRate = 10  # MHz
numEvents  = 10
numChans   = 10

# data file for streaming analysis
# datFile = open('run-%2.0d-mhz-%d-chan-%d-ev.dat' % (sampleRate, numChans, numEvents), 'w+')
# for chan in range(1, numChans+1) : 
#     datFile.write('# channel = %d\n' % chan)
#     dataObj = dataStream(sampleRate*1.0e+6, numEvents)
#     for event in dataObj : 
#         np.savetxt(datFile, (dataObj.timeSamples, dataObj.adcSamples), fmt='%.9f')
#     #     plt.plot(dataObj.timeSamples, dataObj.adcSamples) 
#     # plt.show()
# datFile.close()

# data file for event based analysis
datFile = open('run-%2.0d-mhz-%d-chan-%d-ev.dat' % (sampleRate, numChans, numEvents), 'w+')
for ievent in range(1, numEvents+1) :
	datFile.write('@ event = %d\n' % ievent)
	for chan in range(1, numChans+1) : 
		datFile.write('# channel = %d\n' % chan)
		dataObj = dataStream(sampleRate*1.0e+6, 1)
		for event in dataObj : 
			np.savetxt(datFile, (dataObj.timeSamples, dataObj.adcSamples), fmt='%.9f')
		# 	plt.plot(dataObj.timeSamples, dataObj.adcSamples) 
		# plt.show()
datFile.close()
