__author__ = 'davidgibbs'


# Depending on the configuration,
# we use different methods of inducing
# missing data.

from numpy.random import ranf


# MCAR - missing completely at random.

def mcar(param, pedlist):
    percMissing = float(param["missingParam1"])
    for ped in pedlist:
        for i in range(0,2+ped.getSibNum()):
            for j in range(0, param["numMarkers"]):
                # if we roll the dice and it comes up missing, choose a missing var.
                if ranf() < percMissing:
                    ped.makeMissing(i,j)
    return(pedlist)


def marMF(param, pedlist, target, tp2A, tp2B):
    for ped in pedlist:                                   #for each ped
        for i in range(0,2+ped.getSibNum()):              # for each person
            for j in range(0, param["numMarkers"]):       #   for each marker
                if ped.getGender(i) == target:            #      if they are in the target:
                    if ranf() < tp2A:                     #         use the target specific parameter
                        ped.makeMissing(i,j)              #         one gender will have more missing data.
                else:
                    if ranf() < tp2B:
                        ped.makeMissing(i,j)
    return(pedlist)


def marPop(param, pedlist, target, tp2A, tp2B):
    for ped in pedlist:                                   #for each ped
        for i in range(0,2+ped.getSibNum()):              # for each person
            for j in range(0, param["numMarkers"]):       #   for each marker
                if ped.getPop() == target:            #      if they are in the target:
                    if ranf() < tp2A:                     #         use the target specific parameter
                        ped.makeMissing(i,j)              #         one gender will have more missing data.
                else:
                    if ranf() < tp2B:
                        ped.makeMissing(i,j)
    return(pedlist)


def mar(param, pedlist):
    percMissing = float(param["missingParam1"])
    target = param["missingParam2"]
    targetParam1 = param["missingParam3"]
    tp2A = float(targetParam1) * percMissing
    tp2B = (1.0-float(targetParam1)) * percMissing
    if(target == 1 or target == 2) :
        return(marMF(param, pedlist, target, tp2A, tp2B))
    elif(target == "P1" or target == "P0"):
        return(marPop(param, pedlist, target, tp2A, tp2B))
    else:
        print "Error .. use 1/2 for gender or P0/P1 for population"


def mnarCC(param, pedlist, target, tp2A, tp2B):
    for ped in pedlist:                                   #for each ped
        for i in range(0,2+ped.getSibNum()):              # for each person
            for j in range(0, param["numMarkers"]):       #   for each marker
                if ped.sibCaseBool():                     #      if a CASE is in the ped
                    if ranf() < tp2A:                     #         use the target specific parameter
                        ped.makeMissing(i,j)              #         one gender will have more missing data.
                else:
                    if ranf() < tp2B:
                        ped.makeMissing(i,j)
    return(pedlist)


def mnarHH(param, pedlist, target, tp2A, tp2B):
    for ped in pedlist:                                   #for each ped
        for i in range(0,2+ped.getSibNum()):              # for each person
            for j in range(0, param["numMarkers"]):       #   for each marker
                if sum(ped.getGeno(i,j)) == 1.0:          #      if het is in the ped
                    if ranf() < tp2A:                     #         use the target specific parameter
                        ped.makeMissing(i,j)              #         one gender will have more missing data.
                else:
                    if ranf() < tp2B:
                        ped.makeMissing(i,j)
    return(pedlist)


def mnar(param, pedlist):
    percMissing = float(param["missingParam1"])
    target = param["missingParam2"]
    targetParam1 = param["missingParam3"]
    tp2A = float(targetParam1) * percMissing
    tp2B = (1.0-float(targetParam1)) * percMissing
    if(target == "CC"):
        return(mnarCC(param, pedlist, target, tp2A, tp2B))
    elif(target == "HH"):
        return(mnarHH(param, pedlist, target, tp2A, tp2B))
    else:
        print "Error .. use 0/1 for gender or P0/P1 for population"
    return(pedlist)


def missingEngine(param, peds):
    if param["missingType"] == "mcar":
        return(mcar(param, peds))
    if param["missingType"] == "mar":
        return(mar(param, peds))
    if param["missingType"] == "mnar":
        return(mnar(param, peds))
