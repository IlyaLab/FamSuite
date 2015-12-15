import sys
import string
import gzip

vcfFile = gzip.open(sys.argv[1],"r")
memberType = []

def getNumericCode(genotype):
	if genotype == "0/0" or genotype == "0":
		return "0"
	elif genotype == "0/1" or genotype == "1/0" or genotype == "1":
		return "1"
	elif genotype == "1/1":
		return "2"
	elif genotype == "./." or genotype == ".":
		return "NA"
	
for line in vcfFile:
	vcfColumns = line.strip().split()
	if line.startswith("chr"):
		variantID = vcfColumns[0]+":"+vcfColumns[1]
		numericVector = [getNumericCode(x) for x in vcfColumns[9:]]
		print variantID,'\t','\t'.join(numericVector)
	elif line.startswith("#CHROM"):   #header row
		sampleIDs = vcfColumns[9:]
		for x in sampleIDs:
			if x.find("F") <> -1:
				memberType.append("1")
			elif x.find("M") <> -1:
				memberType.append("2")
			elif x.find("NB") <> -1:
				memberType.append("3")
			else:
				print "unrecognized id format : ",x
	
		print "PED_ID",'\t','\t'.join([x[0:7] for x in sampleIDs])		 
		print "MEMBER_TYPE",'\t','\t'.join(memberType)


vcfFile.close()	
