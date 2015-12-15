__author__ = 'davidgibbs'

import gzip
#from statsmodels.sandbox.stats.multicomp import fdrcorrection0


def calc_benjamini_hochberg_cutoff(p_values, fdr, numTests):
    # find the largest K such that the test in the if statement is true
    lastP = 2.220446049250313e-15
    for k, p_value in enumerate(p_values):
        if (p_value > (float(k+1)/float(numTests))*fdr):
            return(lastP)
        else:
            lastP = p_value
    return(lastP)


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


def minner(z):
    w = z.replace('(','').replace(')','').split(",")
    try:
        y = [float(w[0].strip()), float(w[1].strip())]
    except:
        y = [0.99,0.99]
    return(max(y))


def ciFun(z, mode):
    w = z.replace('(','').replace(')','').split(",")
    try:
        y = [float(w[0].strip()), float(w[1].strip())]
    except:
        y = [0.99,0.99]
    if mode == "cifbat-median":
        return((y[0]+y[1])/2.0)
    else:
        return(max(y))


def isSignif(ps, psi, i, pvalCutOff, mode):
    #epsilon = sys.float_info.epsilon  # 2.220446049250313e-16
    epsilon = 2.220446049250313e-15
    a = ps[i]  - pvalCutOff
    b = psi[i] - pvalCutOff  # should be both >= 0.0 and smaller than epsilon

    if mode == "fbat":
        if a <= epsilon:
            return(i)
    elif mode == "cifbat":
        if a <= epsilon and b <= epsilon:
            return(i)
    elif mode == "cifbat-interval":
        if b <= epsilon:
            return(i)
    elif mode == "cifbat-median":
        if b <= epsilon:
            return(i)

    return(-1)


def parseResults(peds, param, stri, mode):
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

    pvalCutOff = calc_benjamini_hochberg_cutoff(sorted(ps), 0.05, len(ps))  # FDR based on FBAT #

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



def parseResultsORIG(peds, param, stri):
        log = open(param["resultsDir"]+"/log_"+stri+".txt", 'r').read().strip().split("\n")
        jdx = log.index('$')
        causalIdx = [int(a) for a in log[1:jdx]]
        txt = gzip.open(param["resultsDir"]+"/cifbat_"+stri+".out.gz", 'r').read().strip().split("\n")
        ps  = [(z.split("\t"))[23] for z in txt][1:]  # column 23 is the p-values
        ps  = [safeFloat(z) for z in ps]
        psi = [(z.split("\t"))[25] for z in txt][1:]  # column 25 is the p-value interval.
        psi = [minner(z) for z in psi]  # column 25 is the p-value interval.
        psi = [safeFloat(z) for z in psi]
        #qs = fdrcorrection0(ps, alpha=0.05, method='indep', is_sorted=False)
        pvalCutOff = calc_benjamini_hochberg_cutoff(sorted(ps), 0.05, param["numMarkers"])  # FDR based on FBAT #
        predictedIndex = [isSignif(ps, psi, i, pvalCutOff) for i in range(0,len(psi))]
        predictedIndex = [x for x in predictedIndex if x > -1]
        truthTable = ['NA','NA','NA','NA']  # TT  predTcausalF  predFcausalT  FF #
        truthTable[0] = str(len(intersect(causalIdx, predictedIndex)))
        truthTable[1] = str(len(setDiff(causalIdx, predictedIndex)))
        truthTable[2] = str(len(setDiff(predictedIndex,causalIdx)))
        truthTable[3] = str(len(setDiff(set(range(0,len(ps))), union(predictedIndex,causalIdx))))
        return(truthTable)
