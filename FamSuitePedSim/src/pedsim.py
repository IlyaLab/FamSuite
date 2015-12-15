__author__ = 'davidgibbs'

import numpy as np
from numpy.random import randint, ranf
from ped import Ped
from simPheno import genPheno, callPheno

# each family is assigned to one of a number of populations
def choosePop(param):
    n = range(0,param["numPopulations"])
    x = ranf()
    if x < param["popProportions"]:
        return(0)
    else:
        return(1)


def oneProb():
    x = np.random.gamma(shape=2, scale=2)/35.0
    if x > 0.5:
        x = 0.4999
    if x < 0.0:
        x = 0.0001
    return(x)


# each population has a set of marker freqencies.
def markerProbs(param):
    popMarkerProbs = dict()
    for i in range(0, int(param["numPopulations"])):
        probname = "P"+str(i)
        popMarkerProbs[probname] = [oneProb() for j in range(0, param["numMarkers"])]
    return(popMarkerProbs)


# simulate a number of Ped objects, one for each family.
def simped(param):

    # Number of children for each family
    pedNs = randint(low=param["minChildren"], high=param["maxChildren"], size=param["numPeds"])


    # which population they belong to
    pedPops = [choosePop(param) for i in range(0,len(pedNs))]

    # marker proportions by population
    markerFreqs = markerProbs(param)
    param["markerFreqs"] = markerFreqs

    # the ped objects.
    peds = [Ped(pedPops[i], pedNs[i], i, param) for i in range(0,len(pedNs))]

    # compute phenotypes
    (peds, varParams) = genPheno(peds, param)
    callPheno(peds)

    return((peds, varParams))
