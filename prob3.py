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
mu_target, sigma_target = 2, 1 # mean and standard deviation of target distribution
mu_non_target, sigma_non_target = 0.5, 1 # mean and standard deviation of non-target distribution
num_of_trials = 200


def run_trial(mu_target, sigma_target, mu_non_target, sigma_non_target, num_of_trials, size):

    # Create distributions
    target_dist = np.random.normal(mu_target, sigma_target, size_of_distribution)
    non_target_dist = np.random.normal(mu_non_target, sigma_non_target, size_of_distribution)

    # c is the criteria: midpoint of the two means
    c = (mu_target - mu_non_target)/2

    yes_count = 0
    no_count = 0
    hitCount = 0
    faCount = 0
    stimCount = 0
    noStimCount = 0
    orderList=[i for i in range(100)]


    for i in range(num_of_trials):
        chance = np.random.choice(orderList)
        if chance %2 == 0:
            val = np.random.choice(target_dist, 1, replace=False)
            stimCount += 1
            if val > c:
                yes_count += 1
                hitCount +=1
            else:
                no_count += 1
        else:
            val = np.random.choice(non_target_dist, 1, replace=False)
            noStimCount += 1
            if val > c:
                yes_count += 1
                faCount +=1
            else:
                no_count += 1
    
    hitRate = float(hitCount)/float(stimCount)
    faRate = float(faCount)/float(noStimCount)
    print (hitCount)
    print (faCount)
    print (stimCount)
    print (noStimCount)
    print (hitRate)
    print (faRate)

    dPrime = scipy.stats.norm.ppf(hitRate) - scipy.stats.norm.ppf(faRate)

    print(yes_count)
    print(no_count)

    return (yes_count/no_count, dPrime)

print(run_trial(mu_target, sigma_target, mu_non_target, sigma_non_target, num_of_trials, size_of_distribution))
