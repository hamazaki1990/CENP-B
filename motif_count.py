from Bio import SeqIO
from Bio import motifs


def length_sieve(inputf, form):
    rec_iter = SeqIO.parse(inputf, form)
    instances = []
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            if not instances:
                instances.append([seq_record.seq])
            else:
                for i in instances:
                    if len(seq_record.seq) == len(i[0]):
                        i.append(seq_record.seq)
                        break
                else:
                    instances.append([seq_record.seq])
    return instances


for i in range(1, 13):
    inputf = "DXZ1_monomer"+str(i)+".fa"
    instances = length_sieve(inputf, "fasta")
    m = [motifs.create(i) for i in instances]
    print(inputf)
    for x in m:
        print(len(x))
        print(x.counts)
