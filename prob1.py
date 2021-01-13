import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import scipy.io
import struct
from scipy import stats
import pandas as pd
from scipy.stats import norm

from meta_d import MetaD 

# Problem 1 Question 1

FILE_PREFIX = '/Users/joonhwakim/Desktop/git/research/'    # change for personal file system
proj = scipy.io.loadmat(FILE_PREFIX + 'grating2AFC S11.mat')


rRT = proj['data']['responseRT']
rRTlist = list(rRT[0][0][0])
half=len(rRTlist)//2
rRTfirst = rRTlist[:half]
rRTsecond = rRTlist[half:]

meanFirst=np.mean(rRTfirst)
meanSecond=np.mean(rRTsecond)
SE1=stats.sem(rRTfirst)
SE2=stats.sem(rRTsecond)

expGraph = {'First Half': meanFirst,'Second Half':meanSecond, }
expNames = list(expGraph.keys())
expData=list(expGraph.values())
expXvalues=[1,0]
expError = [SE1,SE2]

plt.bar(expXvalues,expData, yerr=expError, ecolor = 'black', color='black', width=0.2, align='center')
plt.xticks(expXvalues, expNames)
plt.ylabel('Mean Reaction Time')
plt.xlabel('Experiment Halves')
plt.title('Mean Reaction Time of Each Half of the Experiment')
plt.show()


#Problem 1 Question 2

t2, p2 = stats.ttest_ind(rRTfirst,rRTsecond)
print "t = ", t2
print "p = ", p2
print "The RTs for the first and the second halves of the experiment differed significantly."

#Problem 1 Question 3

rate = proj['data']['rating']
rateList = list(rate[0][0][0])
con1 = list()
con2 = list()
con3 = list()
con4 = list()


for x,y in enumerate(rateList):
    if y == 1:
        con1.append(rRTlist[x])
    elif y == 2:
        con2.append(rRTlist[x])
    elif y == 3:
        con3.append(rRTlist[x])
    elif y == 4:
        con4.append(rRTlist[x])
    else:
        pass

mean1=np.mean(con1)
mean2=np.mean(con2)
mean3=np.mean(con3)
mean4=np.mean(con4)

conGraph = {'1': mean1,'2':mean2,'3':mean3,'4':mean4}
conNames = list(conGraph.keys())
conData=list(conGraph.values())
conXvalues=[0,2,1,3]

plt.bar(conXvalues,conData, color='black', width=0.2, align='center')
plt.xticks(conXvalues, conNames)
plt.ylabel('Mean Reaction Time')
plt.xlabel('Confidence Ratings')
plt.title('Mean Reaction Time of Each Confidence Rating')
plt.show()

#Problem 1 Question 4

    # Ratios
    #hit_rate = n_Hit/(n_Hit + n_Miss)
    #fa_rate = n_FA/(n_FA + n_CR)

    # Adjusted ratios
    #hit_rate_adjusted = (n_Hit+ 0.5)/((n_Hit+ 0.5) + n_Miss + 1)
    #fa_rate_adjusted = (n_FA+ 0.5)/((n_FA+ 0.5) + n_CR + 1)

    # dprime
    #dprime = scipy.stats.norm.ppf(hit_rate_adjusted) - scipy.stats.norm.ppf(fa_rate_adjusted)

sID = proj['data']['stimID']
sidList = list(sID[0][0][0])
res = proj['data']['response']
resList = list(res[0][0][0])

print len(sidList)
print len(resList)

def d_prime(conLevel): #<------ This is where I'm confused on how to actually look at the data for d' because of the two conditions for hits 1,1 and 0,0 
    cHit=0
    stim1Count=0
    cFA=0
    stim0Count=0
    for x, y in enumerate(rateList):
        if y == conLevel:
                if resList[x]<0:
                    pass
                elif sidList[x]==1:
                    stim1Count=stim1Count+1
                    if resList[x]==1:
                        cHit=cHit+1
                else:
                    stim0Count=stim0Count+1
                    if resList[x]==1:
                        cFA=cFA+1
                    
    

    hitRate = float(cHit)/float(stim1Count)
    faRate = float(cFA)/float(stim0Count)

    dPrime = scipy.stats.norm.ppf(hitRate) - scipy.stats.norm.ppf(faRate)
    return dPrime

d_prime(1)
d_prime(2)
d_prime(3)
d_prime(4)

dGraph = {'1': d_prime(1),'2':d_prime(2),'3':d_prime(3),'4':d_prime(4)}
dNames = list(dGraph.keys())
dData=list(dGraph.values())
dXvalues=[0,2,1,3]

plt.bar(dXvalues,dData, color='black', width=0.2, align='center')
plt.xticks(dXvalues, dNames)
plt.ylabel('d-prime')
plt.xlabel('Confidence Ratings')
plt.title('d-prime of Each Confidence Rating')
plt.show()

#Problem 1 Question 5


sIDl = proj['data']['stimID']
sIDlist = list(sIDl[0][0][0])
half=len(sIDlist)//2
sIDfirst = sIDlist[:half]
sIDsecond = sIDlist[half:]

resL = proj['data']['response']
resList = list(resL[0][0][0])
half=len(resList)//2
resLfirst = resList[:half]
resLsecond = resList[half:]

ratL = proj['data']['rating']
ratList = list(ratL[0][0][0])
half=len(sIDlist)//2
ratFirst = ratList[:half]
ratSecond = ratList[half:]

meta_d_classFirst=MetaD(sIDfirst, resLfirst,ratFirst, 4, 1)
meta_d_classFirst.type2_sdt_sse()
metaDfirst= meta_d_classFirst.meta_d_a

meta_d_classSecond=MetaD(sIDsecond, resLsecond,ratSecond, 4, 1)
meta_d_classSecond.type2_sdt_sse()
metaDsecond= meta_d_classSecond.meta_d_a

mdGraph = {'First Half': metaDfirst,'Second Half':metaDsecond}
mdNames = list(mdGraph.keys())
mdData=list(mdGraph.values())
mdXvalues=[1,0]

plt.bar(mdXvalues,mdData, color='black', width=0.2, align='center')
plt.xticks(mdXvalues, mdNames)
plt.ylabel('meta-d prime')
plt.ylim([0,.7])
plt.xlabel('Experiment Halves')
plt.title('meta-d prime of Each Half of the Experiment')
plt.show()
