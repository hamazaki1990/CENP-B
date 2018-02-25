import pybedtools

for i in range(1, 13):
    bed = pybedtools.BedTool("DXZ1_monomer"+str(i)+"_2000.bed")
    inputfa = "DXZ1.fa"
    outputfa = "DXZ1_monomer"+str(i)+"_2000.fa"
    bed.sequence(fi=inputfa).save_seqs(outputfa)
