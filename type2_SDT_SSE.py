#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 19:18:36 2017

out = type2_SDT(nR_S1, nR_S2)

@author: Vincent Xiao

Given data from an experiment where an observer discriminates between two
stimulus alternatives on every trial and provides confidence ratings,
provides a type 2 SDT analysis of the data.

This is the function that does a standard type 1 SDT analysis on the
raw behavioral data and then does a type 2 SDT analysis using the function
fit_meta_d with d_min = -5, d_grain = .01, d_max = 5

INPUTS

nR_S1, nR_S2 (both are python lists that should have agreeable, even length)

nR_S1 and nR_S2 are vectors containing the total number of responses in
each response category, conditional on presentation of S1 and S2.
size of each array is 2*nRatings, where each element corresponds to a
count of responses in each response category. Response categories are
ordered as follows:
highest conf "S1" ... lowest conf "S1", lowest conf "S2", ... highest conf "S2"

e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
presented, the subject had the following response counts:
responded S1, rating=3 : 100 times
responded S1, rating=2 : 50 times
responded S1, rating=1 : 20 times
responded S2, rating=1 : 10 times
responded S2, rating=2 : 5 times
responded S2, rating=3 : 1 time

The ordering of response / rating counts for S2 should be the same as it
is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
presented, the subject had the following response counts:
responded S1, rating=3 : 3 times
responded S1, rating=2 : 7 times
responded S1, rating=1 : 8 times
responded S2, rating=1 : 12 times
responded S2, rating=2 : 27 times
responded S2, rating=3 : 89 times


OUTPUTS

out['d_a']       : d_a for input data. If s=1, d_a = d'
out['meta_d_a'] : meta_d_a for input data
out['M_ratio']  : meta_d_a / d_a; measure of metacognitive efficiency
out['M_diff']   : meta_d_a - d_a; measure of metacognitive efficiency
out['s']        : ratio of evidence distribution standard deviations assumed for the analysis. 
out['type2_fit']: output of fit_meta_d_SSE for the type 2 SDT fit.

"""
import numpy as np
import sys
import scipy.stats as sstats
from math import sqrt
import fit_meta_d_SSE

def type2_SDT_SSE(nR_S1, nR_S2):
    if not (len(nR_S1) == len(nR_S2) and len(nR_S2) % 2 == 0):
        print "nR_S1 and nR_S2 must be the same length and have an even number of elements."
        sys.exit()

    nRatings = len(nR_S1) / 2

    # Standard SDT analysis
    HR1 = sum(nR_S2[nRatings:]) / float(sum(nR_S2))
    FAR1 = sum(nR_S1[nRatings:]) / float(sum(nR_S1))

    ratingHRs = []
    ratingFARs = []
    for i in range(1, 2*nRatings, 1):
        ratingHRs.append(sum(nR_S2[i:]) / float(sum(nR_S2)))
        ratingFARs.append(sum(nR_S1[i:]) / float(sum(nR_S1)))

    s = 1.0

    # d' and c in terms of S1 distribution standard deviation units
    d_1 = (1/s) * sstats.norm.ppf(HR1) - sstats.norm.ppf(FAR1)
    c_1 = (-1/(1 + s)) * (sstats.norm.ppf(HR1) + sstats.norm.ppf(FAR1))
    c_prime = c_1 / d_1

    # Type 2 SDT analysis

    # Get type 2 HR and FAR for S1 responses
    HR2_rS1 = []
    FAR2_rS1 = []
    for i in range(nRatings - 1):
        HR2_rS1.append(sum(nR_S1[0:(i + 1)]) / float(sum(nR_S1[:nRatings])))
        FAR2_rS1.append(sum(nR_S2[0:(i + 1)]) / float(sum(nR_S2[:nRatings])))

    HR2_rS2 = []
    FAR2_rS2 = []
    for i in range(nRatings + 2, 2*nRatings + 1, 1):
        HR2_rS2.append(sum(nR_S2[(i - 1):]) / float(sum(nR_S2[nRatings:])))
        FAR2_rS2.append(sum(nR_S1[(i - 1):]) / float(sum(nR_S1[nRatings:])))

    d_min = -5
    d_grain = 0.01
    d_max = 5

    fit = fit_meta_d_SSE.fit_meta_d_SSE(HR2_rS1, FAR2_rS1, HR2_rS2, FAR2_rS2, c_prime, s, d_min, d_max, d_grain)

    # Package output as a dictionary
    out = {
        'd_a': d_1 * s * sqrt((2.0 / (1 + s**2))),
        'meta_d_a': fit['meta_d1'] * s * sqrt((2.0 / (1 + s**2))),
        'M_ratio': (fit['meta_d1'] * s * sqrt((2.0 / (1 + s**2)))) / (d_1 * s * sqrt((2.0 / (1 + s**2)))),
        'M_diff': fit['meta_d1'] * s * sqrt((2.0 / (1 + s**2))) - d_1 * s * sqrt((2.0 / (1 + s**2))),
        's': s,
        'type2_fit': fit,
        }

    return out
