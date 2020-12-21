import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import scipy.io
import struct
from scipy import stats

FILE_NAME = 'data.mat'    # change for personal file system
proj = scipy.io.loadmat(FILE_NAME)

rRT = proj['data']['responseRT']

rRTlist = list()
rRTlist = rRT[0][0][0]

half=len(rRTlist)//2
mean1=float(round(np.mean(rRTlist[:half]),2))
mean2=float(round(np.mean(rRTlist[half:]),2))
SE1=float(stats.sem(rRTlist[:half]))
SE2=float(stats.sem(rRTlist[half:]))

print mean1
print mean2
print SE1
print SE2

groupGraph = {'First Half': mean1,'Second Half':mean2}
groupnames = list(groupGraph.keys())
groupdata=list(groupGraph.values())
fig, ax = plt.subplots()
plt.bar(groupnames,groupdata, color='maroon', width=0.4)
plt.show()

#rate = proj['data']['rating']

#t2, p2 = stats.ttest_ind(a,b)