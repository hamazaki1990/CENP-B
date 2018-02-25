from Bio import SeqIO
from Bio import Seq

pJalpha = Seq.Seq("NNNNNNNGPuAAAAGGNN")


def check_pJalpha(inputf, inputfmt):
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
                elif (seq[7] == pJalpha[7] and (seq[8] == "A" or seq[8] == "G")
                      and seq[9:15] == pJalpha[9:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpharev(inputf, inputfmt):
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
                elif (seqr[1:5] == pJalpha[1:5] and
                      seqr[9] == pJalpha[9] and seqr[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_pJalpharevcom(inputf, inputfmt):
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
                elif (seqrc[1:5] == pJalpha[1:5] and seqrc[9] == pJalpha[9] and
                      seqrc[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalphacom(inputf, inputfmt):
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
                elif (seqc[1:5] == pJalpha[1:5] and
                      seqc[9] == pJalpha[9] and seqc[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_pJalpha_1T(inputf, inputfmt):
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
                elif seq[1] == pJalpha[1]:
                    i += 1
                elif (seq[2:5] == pJalpha[2:5] and seq[9] == pJalpha[9] and
                      seq[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_2T(inputf, inputfmt):
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
                elif seq[2] == pJalpha[2]:
                    i += 1
                elif (seq[1] == pJalpha[1] and seq[3:5] == pJalpha[3:5] and
                      seq[9] == pJalpha[9] and seq[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_3C(inputf, inputfmt):
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
                elif seq[3] == pJalpha[3]:
                    i += 1
                elif (seq[1:3] == pJalpha[1:3] and seq[4] == pJalpha[4] and
                      seq[9] == pJalpha[9] and seq[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_4G(inputf, inputfmt):
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
                elif seq[4] == pJalpha[4]:
                    i += 1
                elif (seq[1:4] == pJalpha[1:4] and
                      seq[9] == pJalpha[9] and seq[12:16] == pJalpha[12:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_5A(inputf, inputfmt):
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
                elif seq[9] == pJalpha[9]:
                    i += 1
                elif seq[1:5] == pJalpha[1:5] and seq[12:16] == pJalpha[12:16]:
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_6C(inputf, inputfmt):
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
                elif seq[12] == pJalpha[12]:
                    i += 1
                elif (seq[1:5] == pJalpha[1:5] and seq[9] == pJalpha[9] and
                      seq[13:16] == pJalpha[13:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_7G(inputf, inputfmt):
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
                elif seq[13] == pJalpha[13]:
                    i += 1
                elif (seq[1:5] == pJalpha[1:5] and seq[9] == pJalpha[9] and
                      seq[12] == pJalpha[12] and seq[14:16] == pJalpha[14:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_8G(inputf, inputfmt):
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
                elif seq[14] == pJalpha[14]:
                    i += 1
                elif (seq[1:5] == pJalpha[1:5] and seq[9] == pJalpha[9] and
                      seq[12:14] == pJalpha[12:14] and seq[15] == pJalpha[15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJalpha_9G(inputf, inputfmt):
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
                elif seq[15] == pJalpha[15]:
                    i += 1
                elif (seq[1:5] == pJalpha[1:5] and seq[9] == pJalpha[9] and
                      seq[12:15] == pJalpha[12:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def main():
    inputf = "/users/hamazaki/db/blast/human_ASs_2017.fa"
    print(inputf, "pJalpha intact")
    search_pJalpha = check_pJalpha(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_1")
    search_pJalpha_1T = check_pJalpha_1T(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_1T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_2T")
    search_pJalpha_2T = check_pJalpha_2T(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_2T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_3C")
    search_pJalpha_3C = check_pJalpha_3C(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_3C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_4G")
    search_pJalpha_4G = check_pJalpha_4G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_4G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_5A")
    search_pJalpha_5A = check_pJalpha_5A(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_5A)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_6C")
    search_pJalpha_6C = check_pJalpha_6C(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_6C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_7G")
    search_pJalpha_7G = check_pJalpha_7G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_7G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_8G")
    search_pJalpha_8G = check_pJalpha_8G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_8G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJalpha_9G")
    search_pJalpha_9G = check_pJalpha_9G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJalpha_9G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])


if __name__ == '__main__':
    main()
