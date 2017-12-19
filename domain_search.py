from Bio import SeqIO
from Bio import Seq

CENPB = Seq.Seq("NTTCGNNNNANNCGGGN")
CENPB_3CtoT = Seq.Seq("NTTTGNNNNANNCGGGN")
CENPB_8GtoA = Seq.Seq("NTTCGNNNNANNCGAGN")


def check_CENPB(inputf, inputfmt):
    rec_iter = SeqIO.parse(inputf, inputfmt)
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            i = 0
            while True:
                seq = seq_record.seq[i:i+17]
                if len(seq) < 17:
                    break
                elif seq[1] == CENPB[1]:
                    i += 1
                elif (seq[2:5] == CENPB[2:5] and
                      seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_CENPBrev(inputf, inputfmt):
    rec_iter = SeqIO.parse(inputf, inputfmt)
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            i = 0
            while True:
                seq = seq_record.seq[i:i+17]
                seqr = seq[::-1]
                if len(seq) < 17:
                    break
                elif (seqr[1:5] == CENPB[1:5] and
                      seqr[9] == CENPB[9] and seqr[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_CENPBrevcom(inputf, inputfmt):
    rec_iter = SeqIO.parse(inputf, inputfmt)
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            i = 0
            while True:
                seq = seq_record.seq[i:i+17]
                seqrc = seq.reverse_complement()
                if len(seq) < 17:
                    break
                elif seqrc[1] == CENPB[1]:
                    i += 1
                elif (seqrc[2:5] == CENPB[2:5] and
                      seqrc[9] == CENPB[9] and seqrc[12:15] == CENPB[12:15]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_CENPBcom(inputf, inputfmt):
    rec_iter = SeqIO.parse(inputf, inputfmt)
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            i = 0
            while True:
                seq = seq_record.seq[i:i+17]
                seqc = seq.complement()
                if len(seq) < 17:
                    break
                elif (seqc[1:5] == CENPB[1:5] and
                      seqc[9] == CENPB[9] and seqc[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def main():
    inputf = "DXZ1_HOR_2000.fa"
    print(inputf, "forward")
    search_CENPB = check_CENPB(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB)
        except StopIteration:
            break
        else:
            print(found[0], found[1], found[2])
    print("reverse_complement")
    search_CENPBrevcon = check_CENPBrevcom(inputf, "fasta")
    while True:
        try:
            foundrc = next(search_CENPBrevcon)
        except StopIteration:
            break
        else:
            print(foundrc)
    print("complement")
    search_CENPBcom = check_CENPBcom(inputf, "fasta")
    while True:
        try:
            foundc = next(search_CENPBcom)
        except StopIteration:
            break
        else:
            print(foundc)


if __name__ == '__main__':
    main()
