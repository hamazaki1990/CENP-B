from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Seq

CENPB = Seq.Seq("NTTCGNNNNANNCGGGN")
seqA = Seq.Seq("ATTCGGGCCATTCGGGA")

print(CENPB[1:5])

def test(seq):
    if seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]:
        return "True"
    else:
        return "False"

print(test(seqA))

def check_CENPB(inputf, inputfmt):
    CENPB = Seq.Seq("NTTCGNNNNANNCGGGN")
    for seq_record in SeqIO.parse(inputf, inputfmt):
        i = 0
        while True:
            seq = seq_record.seq[i:i+17]
            if seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]:
                yield [i, seq]
            i += 1


#search_CENPB = next(check_CENPB("DXZ1_HOR_2000.fa", "fasta"))
#print(search_CENPB)
