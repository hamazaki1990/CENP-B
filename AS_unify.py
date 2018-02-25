from Bio import SeqIO


def unify_seq(file1, format1):
        rec = SeqIO.parse(file1, format1)
        uniqueseq = [next(rec)]
        while True:
            try:
                seq = next(rec)
            except StopIteration:
                break
            else:
                for i in uniqueseq:
                    if str(i.seq) == str(seq.seq):
                        break
                else:
                    uniqueseq.append(seq)
        return uniqueseq


for i in range(1, 6):
    inputfile1 = "D11Z1_monomer"+str(i)+".fa"
    outputf1 = "D11Z1_monomer"+str(i)+"unique.fa"
    uniqueseq = unify_seq(inputfile1, "fasta")
    with open(outputf1, "w") as outfile:
        for j in uniqueseq:
            SeqIO.write(j, outfile, "fasta")
