__author__ = 'davidgibbs'

import numpy as np
from numpy import random, median
from math import log, exp
from numpy.random import ranf


def prepOneVar(params, pop, k):
    varParams = dict()
    varParams["ID"]      = k
    varParams["VarID"]   = "RS"+str(abs(random.random()))
    varParams["lambdaS"] = params["riskratios"][pop][0]
    varParams["lambdaG"] = params["riskratios"][pop][1]
    varParams["backMu"]  = params["backgroundPenetrance"][pop][0]
    varParams["backSd"]  = params["backgroundPenetrance"][pop][1]
    varParams["penMu"]   = params["penetranceParams"][pop][0]
    varParams["penSd"]   = params["penetranceParams"][pop][1]
    varParams["f0"]      = abs(random.normal(varParams["backMu"], varParams["backSd"]))
    varParams["f2"]      = varParams["f0"] + abs(random.normal(varParams["penMu"], varParams["penSd"]))
    varParams["f1"]      = (varParams["f0"] + varParams["f2"]) / 2.0
    varParams["beta1"]   = log(varParams["f1"]/varParams["f0"])
    varParams["beta2"]   = log(varParams["f2"]/varParams["f0"])
    varParams["sigma2"]  = 2.0*log(varParams["lambdaS"]/varParams["lambdaG"])
    varParams["mu"]      = log(varParams["f0"]) - varParams["sigma2"]/2.0
    return(varParams)


def getAlphas(ped, params):
    # for each family
    #    for each var
    #          pull out an alpha.
    # so we have a matrix of sibs X vars
    # where the sibs are correlated within the family for a particular var
    # but independent compared to other families.
    pop     = ped.getPopInt()
    nSibs   = ped.getSibs()
    falpha  = abs(random.normal(params["backgroundPenetrance"][pop][0], params["backgroundPenetrance"][pop][1]))
    sigma2  = 2.0*log(params["riskratios"][pop][0]/params["riskratios"][pop][1])
    mu      = log(falpha) - (sigma2/2.0)
    muList   = np.array([mu for i in range(0,nSibs)])
    sigmaMat = np.ndarray(shape=(nSibs,nSibs), dtype=float, order='F')
    sigmaMat.fill(0.5)
    np.fill_diagonal(sigmaMat, sigma2)
    alphaMat = np.random.multivariate_normal(muList,sigmaMat,1)[0]
    return(alphaMat)


def toDummy(ped, j, k):
    d = [0.0,0.0,0.0]
    if ped.getSibGeno(j,k) == (0,0):
        return(np.array([1.0, 0.0, 0.0]))
    if ped.getSibGeno(j,k) == (1,0) or ped.getSibGeno(j,k) == (0,1):
        return(np.array([0.0, 1.0, 0.0]))
    if ped.getSibGeno(j,k) == (1,1):
        return(np.array([0.0, 0.0, 1.0]))
    return(np.array([0.0,0.0,0.0]))  # probably should error our here.


def computePheno(pedi, alphaBlock, causalIdx, varParams):
    ys = np.zeros(pedi.getSibs())
    vp = varParams[pedi.getPop()]  # get the proper var params for this population
    for j in range(0,pedi.getSibs()):
        ys[j] += alphaBlock[j]     # the starting place
        for k in range(0,len(causalIdx)):
            betas = [0.0, vp[k]["beta1"], vp[k]["beta2"]]
            XofG = toDummy(pedi, j, causalIdx[k])
            ys[j] += np.dot(betas, XofG)
    return(ys)


def genPheno(peds, param):
    # each causal variant has some parameters associated with it
    # and each population has a different set of parameters
    varParams = dict()
    for i in range(0,param["numPopulations"]):
        varParamID = "P"+str(i)  # same as what's used by the families.
        varParams[varParamID] = [prepOneVar(param,i, k) for k in range(0,param["numCausal"])]

    # for each ped get the starting place ... alphas[family][sib]
    alphas = [getAlphas(p_i, param) for p_i in peds]

    # index of causal SNPs #
    causalIdx = np.sort(np.random.choice(range(0,param["numMarkers"]),size=param["numCausal"], replace=False))

    for i in range(0,param["numPeds"]):
        popk = peds[i].getPop()
        for j in range(0,peds[i].getSibs()):
            logpi = computePheno(peds[i], alphas[i], causalIdx, varParams)
            peds[i].addNBPheno(logpi)
            peds[i].saveCausalIdx(causalIdx)

    return((peds, varParams))


def callPheno(peds, params):
    for ped in peds:
        for j in range(0, ped.getSibNum()):
            thisPheno = ped.getSibPheno(j)
            prob = min(1.0,exp(thisPheno))  # cap the probability
            if params['modelType'] == 'random':
                ped.setSibPhenoCall(j, float(random.randint(1,3,1)[0]))
            else:
                if ranf() < prob:
                    ped.setSibPhenoCall(j, 2.0)
                else:
                    ped.setSibPhenoCall(j, 1.0)
    return()
