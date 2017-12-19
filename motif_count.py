from Bio import SeqIO
from Bio import motifs


def motifs_count(inputf):
    rec_iter = SeqIO.parse(inputf, "fasta")
    instances = []
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            if not instances:
                instances.append([seq_record])
            else:
                for i in instances:
                    if len(seq_record) == len(i[0]):
                        i.append(seq_record)
                        break
                    else:
                        instances.append([seq_record])
    print(instances)
    m = motifs.create(instances[0])
    print(m)


#for i in range(1, 13):
#    inputf = "DXZ1_monomer"+str(i)+".fa"
motifs_count("/Users/hamazaki/github/mkreads/test_mfa2.fa")
