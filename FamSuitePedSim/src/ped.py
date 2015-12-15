__author__ = 'davidgibbs'


import numpy as np
from numpy.random import ranf,randint

# pedigree object
# represents one family.


class Ped:

    def genOneMarker(self, param, i):
        probs = param["markerFreqs"]
        popProbs = probs[self.popStr]
        y = ranf(2)
        if y[0] < popProbs[i]:
            a = 1
        else:
            a = 0
        if y[1] < popProbs[i]:
            b = 1
        else:
            b = 0
        return(a,b)

    def generateMarkers(self, param):
        # for each marker, get the population frequency, and draw the number of minor alleles
        return([self.genOneMarker(param, i) for i in range(0, param["numMarkers"])])

    def generateOffspring(self, idstr, param):
        genotype = []
        for i in range(0, param["numMarkers"]):
            a = [self.markers["F"][i][0], self.markers["F"][i][1]][randint(0,2)]
            b = [self.markers["M"][i][0], self.markers["M"][i][1]][randint(0,2)]
            genotype.append( (a,b) )
        return(genotype)

    def famMarkers(self, param):
        self.markers["F"] = self.generateMarkers(param)
        self.markers["M"] = self.generateMarkers(param)
        for j in range(0, self.numSibs):
            idstr = "NB"+str(j)
            self.genders[idstr] = randint(low=1,high=3,size=1)[0]
            self.markers[idstr] = self.generateOffspring(idstr, param)

    def __init__(self, pop, numSibs, idnum, param):
        # num peds
        #np.random.seed()
        self.popID = pop        # index into params for the marker freqs for pop K
        self.popStr = "P"+str(pop)
        self.numSibs = numSibs
        self.idnum = idnum
        self.markers = dict()
        self.genders = dict()
        self.famMarkers(param)
        self.nbPheno = dict()
        self.nbPhenoCall = dict()
        self.causalIdx = []
        self.nbCausalCount = dict()

    def __str__(self):
        return ("fam ID: "+ str(self.idnum) +
                "\nnumber of children: " + str(self.numSibs) +
                "\npopulation: " + self.popStr + "\n")

    def addNBPheno(self, ys):
        for j in range(0, len(ys)):
            idstr = "NB"+str(j)
            self.nbPheno[idstr] = ys[j]

    def getCausalSet(self):
        return(self.causalIdx)

    def getSibs(self):
        return(self.numSibs)

    def getPop(self):
        return(self.popStr)

    def getPopInt(self):
        return(self.popID)

    def getID(self):
        return(str(self.idnum))

    def getSibGeno(self, j, k):
        idstr = "NB"+str(j)
        return(self.markers[idstr][k])

    def getSibNum(self):
        return(self.numSibs)

    def getSibPheno(self, j):
        nbIdx = "NB"+str(j)
        return(self.nbPheno[nbIdx])

    def getAllSibPhenos(self):
        return(self.nbPhenoCall.values())

    def getSibGender(self, j):
        nbIdx = "NB"+str(j)
        return(self.genders[nbIdx])

    def getGender(self, i):
        if i == 0:
            return(1)  # Mother
        if i == 1:
            return(2)  # Father
        if i > 1:
            return(self.getSibGender(i-2))


    def setSibPhenoCall(self, j, y):
        nbIdx = "NB"+str(j)
        self.nbPhenoCall[nbIdx] = y
        return()

    def getSibPhenoCall(self, j):
        nbIdx = "NB"+str(j)
        return(self.nbPhenoCall[nbIdx])

    def sibCaseBool(self):
        if 2.0 in self.nbPhenoCall.values():
            return(True)
        else:
            return(False)

    def countCausal(self, idstr):
        idx = 0
        self.nbCausalCount[idstr] = 0.0
        for g in self.markers[idstr]:
            if (sum(g) > 0.0) and (idx in self.causalIdx):
                self.nbCausalCount[idstr] += 1
            idx += 1

    def saveCausalIdx(self, causalIdx):
        self.causalIdx = causalIdx
        for j in range(0, self.numSibs):
            idstr = "NB"+str(j)
            self.countCausal(idstr)

    def getSibCausalCount(self, j):
        nbIdx = "NB"+str(j)
        return(self.nbCausalCount[nbIdx])

    def getPedID(self):
        return(self.idnum)

    def getIndvID(self, i):
        if i == 0:
            return('-M')
        if i == 1:
            return('-F')
        if i > 1:
            return('-NB' + str(i-2))

    def getVarSum(self, jstr, vi):
        var = self.markers[jstr][vi]
        if (sum(var) == -2):
            return("NA")
        return(sum(var))

    def getGeno(self, i, j):
        if (i == 0):
            return(self.markers["M"][j])
        if (i == 1):
            return(self.markers["F"][j])
        if (i > 1):
            nbID = "NB"+str(i-2)
            (self.markers[nbID][j])
        return((-1,-1))


    def makeMissing(self, i, j):
        if (i == 0):
            self.markers["M"][j] = (-1,-1)
        if (i == 1):
            self.markers["F"][j] = (-1,-1)
        if (i > 1):
            nbID = "NB"+str(i-2)
            self.markers[nbID][j] = (-1,-1)
        return()
