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
                elif (seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and
                      seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
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
                elif (seqrc[1:5] == CENPB[1:5] and seqrc[9] == CENPB[9] and
                      seqrc[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
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


def check_CENPB_1T(inputf, inputfmt):
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
                elif (seq[2:5] == CENPB[2:5] and seq[9] == CENPB[9] and
                      seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_2T(inputf, inputfmt):
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
                elif seq[2] == CENPB[2]:
                    i += 1
                elif (seq[1] == CENPB[1] and seq[3:5] == CENPB[3:5] and
                      seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_3C(inputf, inputfmt):
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
                elif seq[3] == CENPB[3]:
                    i += 1
                elif (seq[1:3] == CENPB[1:3] and seq[4] == CENPB[4] and
                      seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_4G(inputf, inputfmt):
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
                elif seq[4] == CENPB[4]:
                    i += 1
                elif (seq[1:4] == CENPB[1:4] and
                      seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_5A(inputf, inputfmt):
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
                elif seq[9] == CENPB[9]:
                    i += 1
                elif seq[1:5] == CENPB[1:5] and seq[12:16] == CENPB[12:16]:
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_6C(inputf, inputfmt):
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
                elif seq[12] == CENPB[12]:
                    i += 1
                elif (seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and
                      seq[13:16] == CENPB[13:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_7G(inputf, inputfmt):
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
                elif seq[13] == CENPB[13]:
                    i += 1
                elif (seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and
                      seq[12] == CENPB[12] and seq[14:16] == CENPB[14:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_8G(inputf, inputfmt):
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
                elif seq[14] == CENPB[14]:
                    i += 1
                elif (seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and
                      seq[12:14] == CENPB[12:14] and seq[15] == CENPB[15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_CENPB_9G(inputf, inputfmt):
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
                elif seq[15] == CENPB[15]:
                    i += 1
                elif (seq[1:5] == CENPB[1:5] and seq[9] == CENPB[9] and
                      seq[12:15] == CENPB[12:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def main():
    inputf = "/users/hamazaki/db/blast/human_ASs_2017.fa"
    print(inputf, "intact")
    search_CENPB = check_CENPB(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_1T")
    search_CENPB_1T = check_CENPB_1T(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_1T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_2T")
    search_CENPB_2T = check_CENPB_2T(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_2T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_3C")
    search_CENPB_3C = check_CENPB_3C(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_3C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_4G")
    search_CENPB_4G = check_CENPB_4G(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_4G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_5A")
    search_CENPB_5A = check_CENPB_5A(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_5A)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_6C")
    search_CENPB_6C = check_CENPB_6C(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_6C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_7G")
    search_CENPB_7G = check_CENPB_7G(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_7G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_8G")
    search_CENPB_8G = check_CENPB_8G(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_8G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("CENPB_9G")
    search_CENPB_9G = check_CENPB_9G(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB_9G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])


if __name__ == '__main__':
    main()
