import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import scipy.io
import struct
from scipy import stats
import pandas as pd
from scipy.stats import norm


data = pd.read_csv("results_4.csv")

wsList = []
gList = []
wsError =0
gError =0

for index, row in data.iterrows():
    if row["stimType"] == 0:
        wsList.append(row["responseRT"])
        wsError = wsError+row["error"]
    elif row["stimType"] ==1:
        gList.append(row["responseRT"])
        gError = gError+row["error"]

wsMeanRT = np.mean(wsList)
gMeanRT = np.mean(gList)
wsSE=stats.sem(wsList)
gSE=stats.sem(gList)

eGraph = {'White Square': wsMeanRT,'Gabor':gMeanRT, }
gNames = list(eGraph.keys())
gData=list(eGraph.values())
gXvalues=[1,0]
gabError = [wsSE,gSE]

plt.bar(gXvalues,gData, yerr=gabError, ecolor = 'black', color='black', width=0.2, align='center')
plt.xticks(gXvalues, gNames)
plt.ylabel('Mean Reaction Time')
plt.xlabel('Stimulus Type')
plt.title('Mean Reaction Time of Each Condition')
plt.show()


wsErrorRate = wsError/len(data.index)
gErrorRate = gError/len(data.index)


eGraph = {'White Square': wsErrorRate,'Gabor':gErrorRate, }
eNames = list(eGraph.keys())
eData=list(eGraph.values())
eXvalues=[1,0]

plt.bar(eXvalues,eData, ecolor = 'black', color='black', width=0.2, align='center')
plt.xticks(eXvalues, eNames)
plt.ylabel('Error Rates')
plt.xlabel('Stimulus Type')
plt.title('Error Rates of Each Condition')
plt.show()
