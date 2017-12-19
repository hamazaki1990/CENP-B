import pybedtools

for i in range(1, 6):
    bed = pybedtools.BedTool("D11Z1_monomer"+str(i)+"_modified.bed")
    inputfa = "D11Z1.fa"
    outputfa = "D11Z1_monomer"+str(i)+".fa"
    bed.sequence(fi=inputfa).save_seqs(outputfa)
