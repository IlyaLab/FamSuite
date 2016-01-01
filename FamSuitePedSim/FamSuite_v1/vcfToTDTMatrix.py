import sys
import string
import gzip

if sys.argv[1].endswith("gz"):
	vcfFile = gzip.open(sys.argv[1],"r")
else:
	vcfFile = open(sys.argv[1],"r")	
#vcfIDFile = open(sys.argv[2],"r")

#vcfIDList = vcfIDFile.read().strip().split()

#memberType = []

def getNumericCode(genotype):
	if genotype == "0/0" or genotype == "0":
		return "0"
	elif genotype == "0/1" or genotype == "1/0" or genotype == "1":
		return "1"
	elif genotype == "1/1":
		return "2"
	elif genotype == "./." or genotype == ".":
		return "NA"
sampleIndex = []
sampleIDs = []
fIds = []
mIds = []
#hdrColumns = []	
for line in vcfFile:
	vcfColumns = line.strip().split()
	if line.startswith("chr"):
		variantID = vcfColumns[0]+":"+vcfColumns[1]
		numericVector = [getNumericCode(vcfColumns[9+idx]) for idx in sampleIndex]
		print variantID,'\t','\t'.join(numericVector)
	elif line.startswith("#CHROM"):   #header row
		sampleIDs = vcfColumns[9:]
		sampleIndex = [idx for idx,val in enumerate(vcfColumns[9:])]
		
		for idx in sampleIndex:
			pedId = vcfColumns[9+idx][vcfColumns[9+idx].find("-")+1:]
			if vcfColumns[9+idx].find("F") <> -1 or vcfColumns[9+idx].find("M") <> -1:
				fIds.append("NA")
				mIds.append("NA")
			elif vcfColumns[9+idx].find("NB") <> -1:
				thisFather = "F-"+pedId
				thisMother = "M-"+pedId
				if thisFather in sampleIDs:	
					fIds.append(thisFather)
				else:
					fIds.append("NA")
				if thisMother in sampleIDs:
					mIds.append(thisMother)
				else:
					mIDs.append("NA")
			else:
				print "unrecognized id format : ",x
		print "SAMPLE_ID","\t","\t".join(sampleIDs)
		print "F_ID","\t","\t".join(fIds)
		print "M_ID","\t","\t".join(mIds)

		#testing
#		for id in sampleIDs:
#			if id not in vcfIDList:
#				print 'Mismatch: ',id	
#				sys.exit(0)

#		pedIDs = [x[0:7] if x.endswith(tuple(['F','M','NB'])) else x[x.find('-')+1:len(x)] for x in sampleIDs]
#		print "PED_ID",'\t','\t'.join(pedIDs)		 
#		print "MEMBER_TYPE",'\t','\t'.join(memberType)


vcfFile.close()	
