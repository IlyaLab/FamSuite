#! /usr/bin/env python
import sys
import string
import gzip
import numpy as np
from classMarker import AutosomalMarker,ChrXMarker

#USAGE:  python scanFBAT.py 
#-fm=<featurematrix.txt> <required>
#-phenotype=<phenotype.txt> <required>
#-offset=<a number between 0 and 1> DEFAULT is 0.5 to give equal and opposite weightage to cases and controls (used only for FBAT)
#-test=tdt OR fbat OR omit to run both
#-version=ext OR std OR omit to run both (NOTE: standard score is computed and reported even with extended version)
#-models=a(dditive) OR d(ominant) OR r(ecessive) OR DEFAULT to additive 
#-out=<output file path> or DEFAULT to tdt.out.gz
#-gender=<NB gender file path> <required> (column 1 is pedID of the NB, column 2 is '1' for male, '2' for female)
#-runs=<integer; no. of runs of random sampling for scanFBAT> OPTIONAL. DEFAULT is 100
#-alpha=<confidence level;floating point number>. OPTIONAL. DEFAULT is 0.95
#INPUT FORMATS
#feature matrix - first row contains trio ids, second row indicates member type: 1 (father), 2(mother), 3(offspring). 
#Cell values set to 0 (ref homozygous), 1 (heterozygous), 2 (non ref homozygous), or NA (missing). First column is row ids. For marker rows, it is <chr:position>
#phenotype file - 2 column file with pedigree ids in first column, affection status in second column (1 = control, 2 = case, NA = unknown)


#ASSUMPTIONS:
#For trios only, multiple offsprings and larger pedigrees are not handled
#NOTE: TDT statistic is always positive, so ref/alt annotation doesn't matter. FBAT computes over or under transmission of the alternate allele(coded as 1)
#NOTE: phenoFile may contain NA too. those pedigrees will not be used in the analysis
#####################################
##GLOBAL VARIABLES and LOOK-UP TABLES
headerColumns = []
sampleID = []
fIDs = []
mIDs = []
NBPhenoDict = {}
NBGenderDict = {}
excludeNB = []
excludeNBForChrX = []

FM_FILENAME=""                    #required
PHENO_FILENAME=""                 #required
OFFSET = 0.5			  #optional	
TEST=""                           #optional
VERSION=""                        #optional
MODELS=[]                         #optional
OUTPUT_FILENAME = "tdt.out.gz"    #optional
GENDER_FILENAME = ""
RUNS = 100
ALPHA = 0.05

fmFile = None
phenoFile = None
outputFile = None
genderFile = None

############################################
#####FUNCTION DEFINITIONS
def readInputArguments():
	global FM_FILENAME
	global PHENO_FILENAME
	global OFFSET
	global TEST
	global VERSION
	global MODELS
	global OUTPUT_FILENAME
	global GENDER_FILENAME
	global RUNS
	global ALPHA
        #PARSE COMMAND LINE ARGS. (3 of the above arguments are required)
        assert (len(sys.argv) >=4 ), 'Insufficient number of arguments.'

        while len(sys.argv) > 1:
                # Command line arguments containing '=' are implicitly options.
                thisArg = sys.argv.pop(1)
                if thisArg.find("=")==-1:
                        print 'Unrecognised argument: '+thisArg
                        sys.exit(1)
                else:
                        name,value = thisArg.split("=")  #split name value pairs
			if __debug__: # ...just to see what's going on.
                                print( "{},{}".format( name, value ) )
                        
			name = name.lower().strip("- ")  #strip hyphens and white spaces
                        if name == "fm":
                                FM_FILENAME = value.strip(" ")
                        elif name == "phenotype":
                                PHENO_FILENAME = value.strip(" ")
			elif name == "offset":
				OFFSET = float(value.strip(" "))
                        elif name == "test":
                                TEST = value.lower().strip(" ")
                        elif name == "version":
                                VERSION = value.lower().strip(" ")
                        elif name == "model":
                                MODELS = value.lower().strip(" ").split(",")
                        elif name == "out":
                                OUTPUT_FILENAME = value.strip(" ")
			elif name == "gender":
				GENDER_FILENAME = value.strip(" ")
			elif name == "runs":
				RUNS = int(value.strip(" "))	
			elif name == "alpha":
				ALPHA = float(value.strip(" "))
                        else:
                                print "unrecognized option:", name
                                sys.exit(1)

        assert (FM_FILENAME <> ""), 'Feature matrix path was not provided'
        assert (PHENO_FILENAME <> ""), 'Phenotype file path was not provided'
	assert (GENDER_FILENAME <> ""), 'Gender file path was not provided'

def createFileObjects():
	global FM_FILENAME
        global PHENO_FILENAME
	global VERSION
        global OUTPUT_FILENAME
	global GENDER_FILENAME

	global fmFile
	global phenoFile
	global outputFile
	global genderFile
	
	#INPUT: feature matrix - first row contains sample ids, second row contains father ids, third row contains mother ids
	#fourth row onwards - Cell values set to 0 (ref homozygous), 1 (heterozygous), 2 (non ref homozygous), or NA (missing)
	if FM_FILENAME.endswith("gz"):
	        fmFile = gzip.open(FM_FILENAME,"r")
	else:
        	fmFile = open(FM_FILENAME,"r")

	#INPUT: phenotype file - 2 column file with pedigree ids in first column, affection status in second column (1 = control, 2 = case)
	phenoFile = open(PHENO_FILENAME,"r")

	#INPUT: GENDER FILE
	if GENDER_FILENAME <> "":
                genderFile = open(GENDER_FILENAME,"r")

	#OUTPUT: output file
	if not OUTPUT_FILENAME.endswith(".gz"):
	        OUTPUT_FILENAME = OUTPUT_FILENAME+".gz"

	outputFile = gzip.open(OUTPUT_FILENAME,"w")

	

def createNBPhenoDict():
	global phenoFile
	global NBPhenoDict
	global excludeNB
        #get these pedigrees from phenotype file into a dictionary
        for line in phenoFile:
                columns = line.strip().split()
                if columns[1].isdigit():
			NBPhenoDict[columns[0]] = int(columns[1])
		else:
			print "Missing phenotype for ",columns[0],". This NB will not be analysed."
			excludeNB.append(columns[0])

def createNBGenderDict():
	global genderFile
	global NBGenderDict
	global excludeNBForChrX	
        #get these pedigrees from phenotype file into a dictionary
        for line in genderFile:
                columns = line.strip().split()
                if columns[1].isdigit():
                	NBGenderDict[columns[0]] = int(columns[1])
                else:
  			print "Missing gender for ",columns[0],". This NB will not be analysed for chrX."
                        excludeNBForChrX.append(columns[0])        

###################################################
######### PROCESSING STARTS HERE ##################

#parse input arguments
readInputArguments()

###open input files for reading and output files for writing
createFileObjects()

#PRINT HEADER line according to options selected by user
outputColumns = ['MarkerID','MAF','n[0/0,0/1,1/1,./.]','nCompleteInformativeCases_MaleNB','nCompleteInformativeCases_FemaleNB','nCompleteInformativeControls_MaleNB','nCompleteInformativeControls_FemaleNB','nCompleteNonInformativeCases','nCompleteNonInformativeControls','nIncompleteInformativeCases_MaleNB','nIncompleteInformativeCases_FemaleNB','nIncompleteInformativeControls_MaleNB','nIncompleteInformativeControls_FemaleNB','nIncompleteNonInformativeCases','nIncompleteNonInformativeControls','nMIE']
if TEST == "" or TEST == "tdt": #both or TDT
	if not MODELS or "a" in MODELS:
		outputColumns.extend(['ChiSq_TDT_Additive','P-value_TDT_Additive'])
		if VERSION == "ext" or VERSION == "":
			outputColumns.extend(['min_ChiSq_rTDT_Additive','min_P-value_rTDT_Additive','max_ChiSq_rTDT_Additive','max_P-value_rTDT_Additive'])
	if "d" in MODELS:    
		outputColumns.extend(['ChiSq_TDT_Dominant','P-value_TDT_Dominant'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['min_ChiSq_rTDT_Dominant','min_P-value_rTDT_Dominant','max_ChiSq_rTDT_Dominant','max_P-value_rTDT_Dominant'])	
	if "r" in MODELS:
		outputColumns.extend(['ChiSq_TDT_Recessive','P-value_TDT_Recessive'])
		if VERSION == "ext" or VERSION == "":	
			outputColumns.extend(['min_ChiSq_rTDT_Recessive','min_P-value_rTDT_Recessive','max_ChiSq_rTDT_Recessive','max_P-value_rTDT_Recessive'])	
	if MODELS and "a" not in MODELS and "d" not in MODELS and "r" not in MODELS: 
		print "Unrecognised models: "+",".join(MODELS)
       		sys.exit(1)


if TEST == "" or TEST == "fbat": #both or FBAT
	if not MODELS or "a" in MODELS:
                outputColumns.extend(['Z_FBAT_Additive','P-value_FBAT_Additive'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['ConfidenceIntervalZ'+str(ALPHA),'ConfidenceIntervalPValues'+str(ALPHA)])
        if "d" in MODELS:
                outputColumns.extend(['Z_FBAT_Dominant','P-value_FBAT_Dominant'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['Confidence_Interval_'+str(ALPHA),'ConfidenceIntervalPValues'+str(ALPHA)])
        if "r" in MODELS:
                outputColumns.extend(['Z_FBAT_Recessive','P-value_FBAT_Recessive'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['Confidence_Interval_'+str(ALPHA),'ConfidenceIntervalPValues'+str(ALPHA)])
        if MODELS and "a" not in MODELS and "d" not in MODELS and "r" not in MODELS:
                print "Unrecognised models: "+",".join(MODELS)
                sys.exit(1)	
outputColumns.append('TotalIncompleteInformativeTrios')
outputFile.write('\t'.join(outputColumns)+'\n')

	
#READ header line containing sampleIDs 
sampleIDs = fmFile.readline().strip().split('\t')[1:]  #1st row of the feature matrix - sample ids (skip 1st column)
fIDs = fmFile.readline().strip().split('\t')[1:]  #second row of the feature matrix is father IDs(skip 1st column)
mIDs = fmFile.readline().strip().split('\t')[1:] #third row of the feature matrix is mother IDs(skip 1st column)


#generate pedigree-phenotype dictionary
createNBPhenoDict()
#pedigree-NBgender dictionary
createNBGenderDict()

## ALL GLOBAL CHECKS and ASSERTS HERE
## TDT must be provided with only case trios, FBAT must have both case and controls
if TEST == "tdt" and 1 in NBPhenoDict.values():
        print 'Pedigree file must contain only case pedigrees for TDT.'
        sys.exit(1)
if TEST == "fbat" and 1 not in NBPhenoDict.values():
        print 'Pedigree file must contain control pedigrees for FBAT.'
        sys.exit(1)



#read feature matrix one line at a time. First two lines have been read above for pedigree ids and member type
for line in fmFile:
	#create new marker object
	#chrM and chrY are not tested
	if line.startswith('chrM') or line.startswith('chr25') or line.startswith('chrY') or line.startswith('chr24'):
		continue	
	elif line.startswith('chrX') or line.startswith('chr23'): 
		thisMarker = ChrXMarker()
	else:  #TODO: more thorough check for validity of data format, like chr numbers??
		thisMarker = AutosomalMarker()
	
	#get vcf columns
 	vcfValues = line.strip().split('\t')

	#assert that all genotype values are numeric #TODO: verify that this assert works
	assert(all(v.isdigit() or v=="NA" for v in vcfValues[1:]))   #1st column is variant id, 2nd onwards are sample genotypes

	
	#set this marker object's sample values 
	thisMarker.markerID = vcfValues[0]
	thisMarker.getPedGenotypes(vcfValues[1:],sampleIDs,fIDs,mIDs)
	
	#Note: thought of checking if chrX marker has heterozygous males, but it's not possible since the feature matrix comes in with encoded genotypes. 
	#So all you can check is whether autosomal chrs have any genotypes other than 0/1/2/NA and chrX has any genotypes other that 0/1/NA for males, 0/1/2/NA for females
	#for chrX
	#TODO: verify that hasValidGenotypes() works
	if isinstance(thisMarker,ChrXMarker) and not thisMarker.hasValidGenotypes(NBGenderDict):
		print 'Invalid genotype found at ',thisMarker.markerID,'. This marker will not be tested.'
		continue
	#for autosomal chromosomes
	elif not isinstance(thisMarker,ChrXMarker) and not thisMarker.hasValidGenotypes():
		print 'Invalid genotype found at ',thisMarker.markerID,'. This marker will not be tested.'
                continue
		
	#COMPUTE ALLELE FREQUENCY
	#compute regardless of whether 'mi' has been selected or not, because allele frequency will be reported in the output
	if isinstance(thisMarker,ChrXMarker):
		thisMarker.computeMAF(NBGenderDict)
	else:
		thisMarker.computeMAF()
	
	try:
		#assert that MAF is always positive
		assert(thisMarker.maf >= 0)
	except(AssertionError):
		print 'Marker ',thisMarker.markerID,', MAF=',thisMarker.maf
		exit(1)	

        #DEBUG
#	print 'Ref and Alt Frequencies'
#       print vcfValues[0],thisMarker.maf
#       raw_input('continue')
	
	#get variant distribution
	if isinstance(thisMarker,ChrXMarker):	
		thisMarker.getVariantDistribution(NBGenderDict)
	else:
		thisMarker.getVariantDistribution()

	#count complete and incomplete case and control trio types and populate corresponding vectors
	thisMarker.populateTrioTypeCountVectors(NBPhenoDict,NBGenderDict)
	
	#concatenate output string and print to output file
        outputColumns = [thisMarker.markerID,str(thisMarker.maf),str(thisMarker.nVariantType),str(thisMarker.nCompleteInformativeCaseTrio_MaleNB),str(thisMarker.nCompleteInformativeCaseTrio_FemaleNB),str(thisMarker.nCompleteInformativeControlTrio_MaleNB),str(thisMarker.nCompleteInformativeControlTrio_FemaleNB),str(thisMarker.nCompleteNonInformativeCaseTrio),str(thisMarker.nCompleteNonInformativeControlTrio),str(thisMarker.nIncompleteInformativeCaseTrio_MaleNB),str(thisMarker.nIncompleteInformativeCaseTrio_FemaleNB),str(thisMarker.nIncompleteInformativeControlTrio_MaleNB),str(thisMarker.nIncompleteInformativeControlTrio_FemaleNB),str(thisMarker.nIncompleteNonInformativeCaseTrio),str(thisMarker.nIncompleteNonInformativeControlTrio),str(thisMarker.nMIE)]			
	# run appropriate tests based on options selected by the user, add appropriate output columns	
	#**************TDT***************************************************************************
	if TEST == "" or TEST == "tdt":
		#ADDITIVE std. TDT-------------------------------------------------------------------
		if not MODELS or "a" in MODELS:
			thisMarker.stdTDT("a")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
			#ADDITIVE mi-TDT--------------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedTDT("a")
				outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])

		#DOMINANT std. TDT----------------------------------------------------------------------------------
		if "d" in MODELS:   
			thisMarker.stdTDT("d")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
			#DOMINANT mi-TDT----------------------------------------------------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedTDT("d")
				outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])
		
		#RECESSIVE std. TDT----------------------------------------------------------------------------------
                if "r" in MODELS:   
                        thisMarker.stdTDT("r")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
                        #RECESSIVE mi-TDT----------------------------------------------------------------------------
			if VERSION == "ext" or VERSION == "":	
                                thisMarker.extendedTDT("r")
                                outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])
	
	if TEST == "" or TEST == "fbat":
		#ADDITIVE std. FBAT-------------------------------------------------------------------
                if not MODELS or "a" in MODELS:
			thisMarker.stdFBAT("a",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #ADDITIVE ext-FBAT------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedFBAT("a",OFFSET,RUNS,ALPHA)
				outputColumns.extend([str(thisMarker.ConfIntervalZ),str(thisMarker.ConfIntervalPValues)])
		
		#DOMINANT std. FBAT-------------------------------------------------------------------
                if "d" in MODELS:
                        thisMarker.stdFBAT("d",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #DOMINANT ext-FBAT------------------------------
			if VERSION == "ext" or VERSION == "":
                                thisMarker.extendedFBAT("d",OFFSET,RUNS,ALPHA)
                                outputColumns.extend([str(thisMarker.ConfIntervalZ),str(thisMarker.ConfIntervalPValues)])	

		#RECESSIVE std. FBAT-------------------------------------------------------------------
                if "r" in MODELS:
                        thisMarker.stdFBAT("r",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #RECESSIVE ext-FBAT-----------------------------
			if VERSION == "ext" or VERSION == "":
                                thisMarker.extendedFBAT("r",OFFSET,RUNS,ALPHA)
                                outputColumns.extend([str(thisMarker.ConfIntervalZ),str(thisMarker.ConfIntervalPValues)])

	#print to outputfile
	outputColumns.append(str(thisMarker.totalIncompleteInformativeTrios))	
	outputFile.write('\t'.join(outputColumns)+'\n')

#TODO: close all files
fmFile.close()
phenoFile.close()
outputFile.close()
genderFile.close()
