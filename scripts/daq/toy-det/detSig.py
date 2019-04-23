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
    def simAdcSignals(self, numSamples, randTimeIndex) :
        self.baseLine   = np.random.rand(int (self.numSamples))*1.0e-2
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
            self.simAdcSignals(self.numSamples, self.randTimeIndex)           
        elif (self.initTime > 0.0 and self.initTime < self.runTime) :
            # simulate the time and adc signals
            self.simTimeSamples(self.initTime, self.initTime + self.windowWidth)
            self.simAdcSignals(self.numSamples, self.randTimeIndex)          
        elif (self.initTime >= self.runTime) : 
            raise StopIteration                                           
