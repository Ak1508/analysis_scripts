#!/usr/bin/python3.6

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import gumbel_r

class dataStream :

    def __init__(self, sampleRate, windowWidth, runTime) :
        # sampling rate (s^-1), read-out window length (s), lenth of run (s)
        self.sampleRate, self.windowWidth, self.runTime = sampleRate, windowWidth, runTime
        # number of samples in read-out window
        self.numSamples = sampleRate*windowWidth
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
            # self.simAdcSignals(self.numSamples, self.randTimeIndex)
            self.simAdcSignals()
        elif (self.initTime > 0.0 and self.initTime < self.runTime) :
            # simulate the time and adc signals
            self.simTimeSamples(self.initTime, self.initTime + self.windowWidth)
            # self.simAdcSignals(self.numSamples, self.randTimeIndex)
            self.simAdcSignals()
        elif (self.initTime >= self.runTime) : 
            raise StopIteration

sampleRate  = 5.0e+6  # s^-1
windowWidth = 1.0e-4  # s
runLength   = 50.0e-4 # s
numChans    = 3
 
for chan in range(1, numChans) : 
    dataObj = dataStream(sampleRate, windowWidth, runLength)
    for event in dataObj : 
        plt.plot(dataObj.timeSamples, dataObj.adcSamples) 
    plt.show()
