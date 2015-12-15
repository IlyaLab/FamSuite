__author__ = 'davidgibbs'


import sys
import os

# writing sets of config files to sweep parameter choices #

def main(argv):

    if len(argv) < 3:
        print("python makeConfig.py path_to_cifbat path_to_write_results path_to_python")
        sys.exit(0)

    thispath = os.getcwd()

    path1 = argv[1]  # path to CIFBAT
    path2 = argv[2]  # path for writing results
    path3 = argv[3]  # path to python
    path4 = argv[4]  # path to missinfoSim/src/main.py


    configs = 24

    totalTrials = 5
    parts = 5


    pedNum = [100,200,300,400,600,800,1000,1200,1400,1600,1800,2000,2200,
              2400,2600,2800,3000,3400,3800,4200,4600,5200,5800,6000]

    label="pedNum"

    runout = open("runlist.txt",'w')
    for i in range(0,configs):
        for j in range(0,parts):
            runout.write("1  " + path3+" "+path4+" "+thispath+"/config"+"_"+label+"_"+str(i)+"_"+str(j)+".txt\n")  # each line needs unique config.
            fout = open("config"+"_"+label+"_"+str(i)+"_"+str(j)+".txt",'w')                                       # config file
            fout.write("outputFileLabel\t"+label+"_"+str(i)+"\n")                                                  # same label for collecting results
            fout.write("pathToCIFBAT\t"+path1+"\n")
            fout.write("resultsDir\t"+path2+label+"_"+str(i)+"_"+str(j)+"\n")
            fout.write("pathToPython\t"+path3+"\n")
            fout.write("totalTrials\t"+str(totalTrials)+"\n")
            fout.write("cifbatTrials\t"+ str(200)+"\n")
            fout.write("numPeds\t"+	str(pedNum[i])+"\n")
            fout.write("ProportionControls\t"+ str(1.0)+"\n")
            fout.write("numMarkers\t"+str(300)+"\n")
            fout.write("numCausal\t"+str(3)+"\n")
            fout.write("numPopulations\t"+str(2)+"\n")
            fout.write("popProportions\t"+str(0.33)+"\n")
            fout.write("minChildren\t"+str(1)+"\n")
            fout.write("maxChildren\t"+str(3)+"\n")
            fout.write("backgroundPenetrance\t"+ "0.001,0.001,0.001,0.001"+"\n")
            fout.write("penetranceParams\t"+ "0.2,0.01,0.2,0.01"+"\n")
            fout.write("riskratios\t"+ "3.0,2.0,3.0,2.0"+"\n")
            fout.write("missingType\t"+ "mcar"+"\n")
            fout.write("missingParam1\t"+str(0.00)+"\n")
            fout.write("missingParam2\t"+"CC"+"\n")
            fout.write("missingParam3\t"+str(0.8)+"\n")
            fout.close()
    runout.close()


if __name__ == "__main__":
   main(sys.argv)


# python ../../../../Code/missInfoSim/src/makeConfig.py ~/Dropbox/Research/Projects/Sim_Gen_CIFBAT/FamSuite1/scanFBAT.py /Volumes/StorageDisk/CIFBAT_DEPOT/tmp2/ python2.7
# python /users/dgibbs/missInfoSim/src/makeConfig.py users/dgibbs/FamSuite1/scanFBAT.py /users/dgibbs/CIFBAT_DEPOT/ /tools/bin/python2.7 /users/dgibbs/missInfoSim/src/main.py
