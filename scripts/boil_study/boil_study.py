
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

current = []
hQBcm1  = []
hQBcm2  = []
                   
hRunTime      = []
hScinTrigs    = []
hElRealTrigs  = []
hElCleanTrigs = []

#with open('hms_LH2_boil.data') as fobj:
#with open('shms_LH2_boil.data') as fobj:
#with open('hms_LD2_boil.data') as fobj:
#with open('shms_LD2_boil.data') as fobj:
#with open('hms_0p5C_boil.data') as fobj:
#with open('shms_0p5C_boil.data') as fobj:
#with open('hms_10cm_al_dummy_boil.data') as fobj:
with open('shms_10cm_al_dummy_boil.data') as fobj:
    for line in fobj:
        row = line.split()
        current.append(row[1])
        hQBcm1.append(row[2])
        hQBcm2.append(row[3])
        hRunTime.append(row[4])
        hScinTrigs.append(row[5])
        hElRealTrigs.append(row[6])
        hElCleanTrigs.append(row[7])

current.pop(0)
hQBcm1.pop(0)
hQBcm2.pop(0)
hRunTime.pop(0)
hScinTrigs.pop(0)
hElRealTrigs.pop(0)
hElCleanTrigs.pop(0)

fcurrent  = [float(i) for i in current]
fhQBcm1   = [float(i) for i in hQBcm1]
fhQBcm2   = [float(i) for i in hQBcm2]
fhRunTime = [float(i) for i in hRunTime]

fhScinTrigs    = [float(i) for i in hScinTrigs]
fhElRealTrigs  = [float(i) for i in hElRealTrigs]
fhElCleanTrigs = [float(i) for i in hElCleanTrigs]

hQSum = np.add(fhQBcm1, fhQBcm2)
hAvgQ = np.divide(hQSum, 2.)

# Only utilize BCM1
hScinYield    = np.divide(fhScinTrigs, fhQBcm1)
hElRealYield  = np.divide(fhElRealTrigs, fhQBcm1)
hElCleanYield = np.divide(fhElCleanTrigs, fhQBcm1)

# Normalize to 10 uA current
hScinNormYield    = hScinYield/hScinYield[0]
hElRealNormYield  = hElRealYield/hElRealYield[0]
hElCleanNormYield = hElCleanYield/hElCleanYield[0]

plt.plot(current, hScinNormYield, 'ro', label='3/4', markersize=8)
plt.plot(current, hElRealNormYield, 'bs', label='EL-REAL', markersize=8)
plt.plot(current, hElCleanNormYield, 'g^', label='EL-CLEAN', markersize=8)

plt.ylabel('Normalized Yield', fontsize=20)
plt.xlabel('Beam Current uA', fontsize=20)
plt.title('SHMS 10 cm Al Dummy Boiling Studies', fontsize=20)

plt.ylim(0.95, 1.25)
plt.xlim(0, 70)

plt.legend(loc=2)

plt.show()


