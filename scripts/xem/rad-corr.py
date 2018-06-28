#!/usr/bin/python

# Import various modules
import os, sys, time
import string
import numpy as np
import matplotlib.pyplot as plt

import numpy.polynomial.polynomial as plyf
import numpy.polynomial.legendre as lf

ebeam, eprime, theta, xbj, q2, w2, rcd, rcc, ratio, cc = np.loadtxt('rc_f1f2_ineft.txt', skiprows = 1, unpack=True)

xfr    = np.linspace(0.2, 1.25, 10000)
interp = np.interp(xfr, xbj, ratio)

def get_rcf(x):
    index = (np.abs(xfr - x)).argmin()
    return interp[index]

offset = 0.0125
xbj_bin = np.linspace(offset, 1.5-offset, 60)

rcf_list = []
for xval in xbj_bin:
    rcf_list.append(get_rcf(xval))
rcf = np.asarray(rcf_list)

plt.plot(xbj, ratio, 'bo', xfr, interp, 'r-')
plt.show()
