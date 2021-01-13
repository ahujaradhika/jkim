import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import scipy.io
import struct
from scipy import stats
import pandas as pd
from scipy.stats import norm

# trial parameters
size_of_distribution = 5
mu_non_target, sigma_non_target = 0, 1 # mean and standard deviation of nTarget distribution
mu_target, sigma_target = 1.5, 1.5 # mean and standard deviation of non-nTarget distribution
num_of_trials = 200
# c is the criteria: midpoint of the two means
crit = ((mu_target-mu_non_target)/2)


def run_trial(mu_non_target, sigma_non_target, mu_target, sigma_target, num_of_trials, size, crit):

    stim = np.random.rand(num_of_trials) > .5

    # Create distributions
    nTarget = np.random.normal(mu_non_target, sigma_non_target, size=stim.size)
    nTarget[stim] = np.random.normal(mu_target, sigma_target, size=stim.sum())
    df=pd.DataFrame({"trial": range(len(nTarget)), "stim": stim, "nTarget": nTarget})

    df['response'] = df.nTarget > crit

    hit = df.response[df.stim]
    miss = ~df.response[df.stim]
    fa = df.response[~df.stim]
    cr = ~df.response[~df.stim]

    # print("Hit rate: {:.2f}".format(hit.mean()))
    # print("Miss rate: {:.2f}".format(miss.mean()))
    # print("False alarm rate: {:.2f}".format(fa.mean()))
    # print("Correct rejection rate: {:.2f}".format(cr.mean()))

    dprime = stats.norm.ppf(hit.mean()) - stats.norm.ppf(fa.mean())
    # print("d prime: {:.2f}".format(dprime))

    c = -(stats.norm.ppf(hit.mean()) + stats.norm.ppf(fa.mean()))/2.0
    # print("c: {:.2f}".format(c))

    return (stim, nTarget, c)

result_200 = run_trial(mu_non_target, sigma_non_target, mu_target, sigma_target, 200, size_of_distribution, crit)
stim = result_200[0]
nTarget = result_200[1]
c = result_200[2]

criterions = np.linspace(crit+(2*sigma_target), crit-(2*sigma_target), 5)
hit_rates = [(nTarget[stim] > c).mean() for c in criterions]
fa_rates = [(nTarget > c).mean() for c in criterions]

plt.plot(fa_rates, hit_rates, color='black')
plt.ylim([0,1])
plt.xlim([0,1])
plt.ylabel('Hit Rate')
plt.xlabel('FA Rate')
plt.title('ROC Plot w/ 200 trials')
plt.show()


result_30 = run_trial(mu_non_target, sigma_non_target, mu_target, sigma_target, 30, size_of_distribution, crit)
stim = result_30[0]
nTarget = result_30[1]
c = result_30[2]

criterions = np.linspace(crit+(2*sigma_target), crit-(2*sigma_target), 5)
hit_rates = [(nTarget[stim] > c).mean() for c in criterions]
fa_rates = [(nTarget > c).mean() for c in criterions]

plt.plot(fa_rates, hit_rates, color='black')
plt.ylim([0,1])
plt.xlim([0,1])
plt.ylabel('Hit Rate')
plt.xlabel('FA Rate')
plt.title('ROC Plot w/ 30 trials')
plt.show()