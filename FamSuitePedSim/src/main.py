__author__ = 'davidgibbs'

import os
import sys
import subprocess
from pedsim import simped
from paramReader import readTheParams
from simDataOut import writeData
from parseResults import parseResults
from missingEngine import missingEngine
from balancer import balancer

# This bit was for running on our cluster #
# warnings were causing the run to exit #
try:
    from Crypto.pct_warnings import PowmInsecureWarning
    import warnings
    warnings.simplefilter("ignore", PowmInsecureWarning)
except:
    pass


def runsim(param, stri):
    (peds, varParams) = simped(param)
    (peds, param)  = balancer(peds, param)
    peds = missingEngine(param, peds)
    writeData(peds, param, varParams, stri)
    print("data written")
    return((peds, varParams))


def runCifbat(param, stri):
    cmd = (param["pathToPython"] + " " + param["pathToCIFBAT"] +
           " -fm=" + param["resultsDir"]+"/genotype_"+stri+".tsv" +
           " -phenotype=" + param["resultsDir"]+"/pheno_"+stri+".tsv" +
           " -gender=" +  param["resultsDir"]+"/gender_"+stri+".tsv" +
           " -runs=" + str(param["cifbatTrials"]) +
           " -out=" + param["resultsDir"]+"/cifbat_"+stri+".out")
    print(cmd)
    subprocess.call(cmd, shell=True)


def main(argv):
    # build param from argv .. or read a parameter file
    if len(argv) <= 1:
        print("Please use a config file, and run the program like:")
        print("python2.7 src/main.py configFile.txt")
        sys.exit(1)
    param  = readTheParams(argv[1])
    outdir = param["resultsDir"]
    label  = param["outputFileLabel"]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.chmod(outdir, 0755)
    resOut = open(outdir+"out.txt",'w')
    for ti in range(0,param["totalTrials"]):
        (peds, varParams) = runsim(param, str(ti))
        runCifbat(param, str(ti))
        truthTable = parseResults(peds, param, str(ti), "cifbat")
        truthTable += [str(param['popPrev'][0]), str(param['popPrev'][1]), str(param['sampleSize'][0]), str(param['sampleSize'][1])]
        outputLine = [param['outputFileLabel']] + truthTable
        print(outputLine)
        resOut.write("\t".join(outputLine)+"\n")
    resOut.close()


if __name__ == "__main__":
   main(sys.argv)
