__author__ = 'davidgibbs'


from numpy.random import choice


def getPopPrev(param, peds):
    K = []
    for pi in range(0,param['numPopulations']):
        K.append([])
    for i, pi in enumerate(peds):
        thisK = pi.getPopInt()
        K[thisK] += pi.getAllSibPhenos()
    PopPrev = []
    for j in range(0,len(K)):
        PopPrev.append(float(K[j].count(2.0))/float(len(K[j])))
        print("Pop"+str(j)+ "  Size: " +str(len(K[j])) + " K: " + str(float(K[j].count(2))/float(len(K[j]))))
    return(PopPrev)


def balancer(peds, param):
    popPrev = getPopPrev(param, peds)
    param['popPrev'] = popPrev
    newPeds = []
    unaffectedList = []
    N = param['numPeds']
    numCases = sum([int(2.0 in pi.getAllSibPhenos()) for pi in peds])
    numControls = int(numCases * float(param['ProportionControls']))
    param['sampleSize'] = [numCases, numControls]
    for i, pi in enumerate(peds):
        if 2.0 in pi.getAllSibPhenos(): # get all the affected peds
            newPeds.append(pi)
        else:                             # there's a chance to add the unaffected ped
            unaffectedList.append(i)
    if numControls > len(unaffectedList) :
        print ("!!!  WARNING: size of the unaffected list < numControls requested. !!!")
        numControls = min(len(unaffectedList), numControls)

    theseUnaffecteds = choice(unaffectedList, size=numControls, replace=False)
    for ti in theseUnaffecteds:
        newPeds.append(peds[ti])

    return((newPeds, param))
