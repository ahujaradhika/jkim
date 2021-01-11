#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:38:56 2017

meta d' calculation helper function: trials2counts

@author: Vincent

This is a python script that contains the trials2counts function.

Given data from an experiment where an observer discriminates between
two stimulus alternatives on every trial and provides confidence ratings,
convert trial by trial experimental information for N trials into response
counts.

INPUTS:
stimID    : 1 * N vector in the form of python list. 
            stimID[i] = 0 --> stimulus on i'th trial was S1.
            stimID[i] = 1 --> stimulus on i'th trial was S2.

response  : 1 * N vector in the form of python list.
            response[i] = 0 --> response on i'th trial was "S1".
            response[i] = 1 --> response on i'th trial was "S2".

rating    : 1 * N vector in the form of python list.
            rating[i] = X --> rating on i'th trial was X.
            X must be in the range 1 <= X <= nRatings.
            
N.B. All trials where stimID is not 0 or 1, response is not 0 or 1, or
rating is not in the range [1, nRatings], are omitted from the response
count.

nRatings  : total # of available subjective ratings available for the
            subject. e.g. if subject can rate confidence on a scale of 1-4,
            then nRatings = 4

optional inputs

padCells : if set to 1, each response count in the output has the value of
           padAmount added to it. Padding cells is desirable if trial counts 
           of 0 interfere with model fitting.
           if set to 0, trial counts are not manipulated and 0s may be
           present in the response count output.
           default value for padCells is 0.

padAmount: the value to add to each response count if padCells is set to 1.
           default value is 1/(2*nRatings)


OUTPUTS
nR_S1, nR_S2
these are vectors containing the total number of responses in
each response category, conditional on presentation of S1 and S2.

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
"""
import sys

def trials2counts(stimID, response, rating, nRatings, padCells=False, padAmount=None):
    if not(len(stimID) == len(response) and len(stimID) == len(rating)):
        print "stimID, response, and rating input vectors must have the same lengths!"
        sys.exit()

    # Filter out bad trials
    f = [True] * len(stimID)
    for i in range(len(stimID)):
        if not (((stimID[i] == 0 or stimID[i] == 1) and (response[i] == 0
            or response[i] == 1)) and (rating[i] >= 1 and rating[i] <= nRatings)):
            f[i] == False

    stimID = [item for item in stimID if f[stimID.index(item)] == True]
    response = [item for item in response if f[response.index(item)] == True]
    rating = [item for item in rating if f[rating.index(item)] == True]

    # Set default inputs
    if padAmount is None:
        padAmount = 1 / (2*nRatings)

    # Compute response counts
    nR_S1 = []
    nR_S2 = []

    # S1 responses
    for i in range(nRatings, 0, -1):
        temp = 0
        for j in range(len(stimID)):
            if stimID[j] == 0 and response[j] == 0 and rating[j] == i:
                temp += 1
        nR_S1.append(temp)
        temp = 0
        for j in range(len(stimID)):
            if stimID[j] == 1 and response[j] == 0 and rating[j] == i:
                temp += 1
        nR_S2.append(temp)

    # S2 responses
    for i in range(1, nRatings + 1, 1):
        temp = 0
        for j in range(len(stimID)):
            if stimID[j] == 0 and response[j] == 1 and rating[j] == i:
                temp += 1
        nR_S1.append(temp)
        temp = 0
        for j in range(len(stimID)):
            if stimID[j] == 1 and response[j] == 1 and rating[j] == i:
                temp += 1
        nR_S2.append(temp)

    # Pad response counts to avoid zeros
    if padCells:
        nR_S1 = [elem + padAmount for elem in nR_S1]
        nR_S2 = [elem + padAmount for elem in nR_S2]

    return [nR_S1, nR_S2]