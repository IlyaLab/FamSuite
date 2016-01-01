__author__ = 'davidgibbs'

import sys

def stringToList(x):
    x1 = [float(z) for z in x.split(",") if z not in [' ', ',']]
    x2 = [  [x1[0],x1[1]], [x1[2],x1[3]]   ]
    return(x2)


def readTheParams(tofile):
    param = dict()
    ptxt = open(tofile, 'r').read().strip().split("\n")
    ptxt = [z.strip() for z in ptxt if z[0] != '#']
    for pi in ptxt:
        if ',' in pi:
            pip = pi.split("\t")
            param[pip[0]] = stringToList(pip[1])
        else:
            pip = pi.split("\t")
            try:
                if '.' in pip[1]:
                    param[pip[0]] = float(pip[1])
                else:
                    param[pip[0]] = int(pip[1])
            except:
                param[pip[0]] = pip[1]
    return(param)
