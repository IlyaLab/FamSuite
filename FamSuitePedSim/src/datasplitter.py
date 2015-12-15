# data splitter #


fin = open("/Volumes/StorageDisk/Epistasis/combined_data/gnmc_data_combined_regions.txt",'r')
thisChr = ''

outfilePrefix = "/Volumes/StorageDisk/Epistasis/combined_data/"
outfilePostfix = ".txt"

for line in fin:
    bits = split(line[1:128], "\t")
    bit1 = split(bits[0], ':')
    bit2 = split(bit[3], '_')
    if thisChr != bit2:  # if we move onto a new chromosome
        thisChr = bit2
        print(thisChr)
        outfile = outfilePrefix+thisChr+outfilePostfix
        fout = open(outfile,'w')
    fout.write(line)
