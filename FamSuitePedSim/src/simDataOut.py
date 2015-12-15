__author__ = 'davidgibbs'



def writeData(peds, param, varParams, stri):

    fout = open(param["resultsDir"]+"/log_"+stri+".txt", 'w')
    logstr = "#Causal Markers\n"
    causalSet = peds[0].getCausalSet()
    for i in range(0, param["numCausal"]):
        logstr += str(causalSet[i]) + "\n"
    fout.write(logstr.strip()+"\n")
    fout.write("$\n")

    fout.write("\nCases   : "+str(param['sampleSize'][0])+"\n")
    fout.write("\nControls: "+str(param['sampleSize'][1])+"\n")
    for pop in range(0,param["numPopulations"]):
        fout.write("\nPOP"+str(pop)+"\n")
        fout.write("   Pop Prev: " + str(param['popPrev'][pop])+"\n")
        for k in range(0, param["numCausal"]):
            fout.write("marker_" + str(causalSet[varParams[("P"+str(pop))][k]["ID"]]) +"\n")
            fout.write("f0: "+str(varParams[("P"+str(pop))][k]["f0"]) +"\n")
            fout.write("f1: "+str(varParams[("P"+str(pop))][k]["f1"]) +"\n")
            fout.write("f2: "+str(varParams[("P"+str(pop))][k]["f2"]) +"\n")
            fout.write("beta1: "+str(varParams[("P"+str(pop))][k]["beta1"]) +"\n")
            fout.write("beta2: "+str(varParams[("P"+str(pop))][k]["beta2"]) +"\n")

    fout.close()

    fout = open(param["resultsDir"]+"/genotype_"+stri+".tsv", 'w')
    # write PED ID line
    fout.write("SAMPLE_ID\t")
    idstr = ""
    for i in range(0,len(peds)):
        for j in range(0,2+peds[i].getSibNum()):
            idstr += ("PED" + str(peds[i].getPedID()) + peds[i].getIndvID(j) +"\t")
    fout.write(idstr.strip() + "\n")

    fout.write("F_ID\t")  # father of the NBs
    idstr = ""
    for i in range(0,len(peds)):
        for j in range(0,2+peds[i].getSibNum()):
            if j > 1:
                idstr += ("PED" + str(peds[i].getPedID()) + "-F" +"\t")
            else:
                idstr += ("NA"+"\t")
    fout.write(idstr.strip() + "\n")

    fout.write("M_ID\t")  # father of the NBs
    idstr = ""
    for i in range(0,len(peds)):
        for j in range(0,2+peds[i].getSibNum()):
            if j > 1:
                idstr += ("PED" + str(peds[i].getPedID()) + "-M" +"\t")
            else:
                idstr += ("NA"+"\t")
    fout.write(idstr.strip() + "\n")

    # write out the var dat
    for vi in range(0,param["numMarkers"]):  # for each marker
        vard = "marker_" + str(vi) + "\t"
        for i in range(0,len(peds)):            # for each pedigree
            vard += (str(peds[i].getVarSum("F", vi)) + "\t")
            vard += (str(peds[i].getVarSum("M", vi)) + "\t")
            for j in range(0, peds[i].getSibNum()):
                x=peds[i].getVarSum("NB"+str(j), vi)
                vard += (str(peds[i].getVarSum("NB"+str(j), vi)) + "\t")
        fout.write(vard.strip()+"\n")
    fout.close()

    fout = open(param["resultsDir"]+"/pheno_"+stri+".tsv",'w')
    for i in range(0,len(peds)):
        for j in range(0, peds[i].getSibNum()):
            phen = "PED" + str(peds[i].getPedID()) + peds[i].getIndvID(2+j) +"\t"+ str(int(peds[i].getSibPhenoCall(j))) + "\n"
            fout.write(phen)
    fout.close()

    fout = open(param["resultsDir"]+"/gender_"+stri+".tsv",'w')
    for i in range(0,len(peds)):
        for j in range(0, peds[i].getSibNum()):
            gend = "PED" + str(peds[i].getPedID()) + peds[i].getIndvID(2+j) +"\t"+ str(peds[i].getSibGender(j)) + "\n"
            fout.write(gend)
    fout.close()
