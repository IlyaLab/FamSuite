__author__ = 'davidgibbs'

import gzip
import os
import sys


def calc_benjamini_hochberg_cutoff(p_values, fdr, numTests):
    """Computes the p-value cutoff for a particular false discovery rate (fdr).

    Args:
        p_values: list of p-values (0 < x < 1)
        fdr: 0.05 for 5% FDR cutoff
        numTests: number of tests performed.

    Returns:
        A float that can be used as a threshold.

    """
    lastP = 2.220446049250313e-15
    for k in range(1,numTests+1):
            # when the cut off becomes larger than P... then we return it
            # and take only things less than that.
        if (p_values[k-1] <= (float(k) / float(numTests)) * fdr):
            lastP = ((float(k) / float(numTests)) * fdr)  # if less, then potentially use it.
        else:
            return(lastP) # if greater, take the last P
    return(lastP)
#    for k, p_value in enumerate(p_values):
#        if (p_value > (float(k + 1) / float(numTests)) * fdr):
#            return (lastP)
#        else:
#            lastP = p_value
#    return (lastP)


def intersect(a, b):
     return list(set(a) & set(b))

def setDiff(a, b):
    return set(a).difference(set(b))

def union(a, b):
    return set(a).union(set(b))


def safeFloat(x):
    if x == 'NA':
        return(1.0)
    else:
        return(float(x))


def ciFun(z, mode):
    w = z.replace('(','').replace(')','').split(",")
    try:
        y = [safeFloat(w[0].strip()), safeFloat(w[1].strip())]
        if mode == "cifbat-median":
            return((y[0]+y[1])/2.0)
        elif mode == "cifbat-min":
            return(min(y))
        else:
            return(max(y))
    except:
        return(0.99)


def isSignif(ps, psi, i, pvalCutOff, mode):
    #epsilon = sys.float_info.epsilon  # 2.220446049250313e-16
    epsilon = 2.220446049250313e-15
    a = ps[i]  - pvalCutOff
    b = psi[i] - pvalCutOff  # should be both > 0.0 and smaller than epsilon

    if mode == "fbat":
        if a < epsilon:
            return(i)
    elif mode == "cifbat":
        if a < epsilon and b < epsilon:
            return(i)
    elif mode == "cifbat-interval" or mode == "cifbat-min":
        if b < epsilon:
            return(i)
    elif mode == "cifbat-median":
        if b < epsilon:
            return(i)

    return(-1)


def parseResults(param, stri, mode):
    log = open(param["resultsDir"]+"/log_"+stri+".txt", 'r').read().strip().split("\n")
    jdx = log.index('$')
    kdx = [i for i,txt in enumerate(log) if 'Cases' in txt][0]
    caseSize = (log[kdx].split(":")[1]).strip()
    causalIdx = [int(a) for a in log[1:jdx]]
    fileInName = param["resultsDir"]+"/cifbat_"+stri+".out.gz"
    txt = gzip.open(fileInName, 'r').read().strip().split("\n")

    ps  = [(z.split("\t"))[23] for z in txt][1:]  # column 23 is the p-values
    ps  = [safeFloat(z) for z in ps]
    psi = [(z.split("\t"))[25] for z in txt][1:]  # column 25 is the p-value interval.
    psi = [ciFun(z,mode) for z in psi]  # column 25 is the p-value interval.

    pvalCutOff = calc_benjamini_hochberg_cutoff(sorted(ps), param["cutoff"], len(ps))  # FDR based on FBAT #

    predictedIndex = [isSignif(ps, psi, i, pvalCutOff, mode) for i in range(0,len(psi))]
    predictedIndex = [x for x in predictedIndex if x > -1]

    truthTable = ['NA','NA','NA','NA', caseSize]  # TT  predTcausalF  predFcausalT  FF #
    truthTable[0] = str(len(intersect(causalIdx, predictedIndex)))
    truthTable[1] = str(len(setDiff(causalIdx, predictedIndex)))
    truthTable[2] = str(len(setDiff(predictedIndex,causalIdx)))
    truthTable[3] = str(len(setDiff(set(range(0,len(ps))), union(predictedIndex,causalIdx))))
    try:
        sens = str(float(truthTable[0])/(float(truthTable[0])+float(truthTable[1])))
    except:
        sens = "0.0"
    try:
        ppv  = str(float(truthTable[0])/(float(truthTable[0])+float(truthTable[2])))
    except:
        ppv = "0.0"
    truthTable.append(sens)
    truthTable.append(ppv)
    return(truthTable)


def main(argv):
    param = dict()
    if len(argv) < 3:
        print("ARG1: fbat/cifbat  ARG2: FDR cut off  ARG3: root dir  ARG4: write dir")
        sys.exit(0)
    mode    = argv[1]
    cutoff  = argv[2]
    rootDir = argv[3]
    writDir = argv[4]
    fout = open(writDir,'w')
    for dirName, subdirList, fileList in os.walk(rootDir, topdown=False):
        print('Found directory: %s' % dirName)
        if 'pedNum' in dirName:
            param["resultsDir"] = dirName
            param["cutoff"] = float(cutoff)
            label = dirName.split("/")
            missing = label[len(label)-1]
            missing = missing.split("_")[1]
            scenario = label[len(label)-2]
            for i in range(0,10):
                try:
                    tt = parseResults(param, str(i), mode)
                    outputLine = [scenario, missing] + tt
                    fout.write("\t".join(outputLine)+"\n")
                except:
                    pass
        elif 'Trials' in dirName:
            param["resultsDir"] = dirName
            label = dirName.split("/")
            label = label[len(label)-1]
            label = "trial_" + label.split("_")[1]+ "_" + label.split("_")[2]
            for i in range(0,10):
                try:
                    tt = parseResults(param, str(i), mode)
                    outputLine = [label] + tt
                    fout.write("\t".join(outputLine)+"\n")
                except:
                    pass

    fout.close()


if __name__ == "__main__":
   main(sys.argv)
