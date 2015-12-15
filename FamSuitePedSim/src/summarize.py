__author__ = 'davidgibbs'


def summaryByVarCount (peds, param):
    Aff = [0,0,0]
    Non = [0,0,0]
    for pi in peds:
        for j in range(0, pi.getSibNum()):
            k = pi.getSibCausalCount(j)
            if k == 0:
                l = 0
            elif k == 1:
                l = 1
            else:
                l = 2
            if pi.getSibPhenoCall(j) == 1.0:
                Aff[l] += 1
            else:
                Non[l] += 1
    print("Affected 0: " + str(Aff[0]))
    print("Non      0: " + str(Non[0]) +"\n")
    print("Affected 1: " + str(Aff[1]))
    print("Non      1: " + str(Non[1]) +"\n")
    print("Affected N: " + str(Aff[2]))
    print("Non      N: " + str(Non[2]) +"\n")

