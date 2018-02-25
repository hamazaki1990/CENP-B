from Bio import SeqIO
from Bio import Seq

pJa = Seq.Seq("NNNNNNNGPAAAAGGNN")


def check_pJa(inputf, inputfmt):
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
                elif (seq[7] == pJa[7] and (seq[8] == "A" or seq[8] == "G")
                      and seq[9:15] == pJa[9:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_1G(inputf, inputfmt):
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
                elif seq[7] == pJa[7]:
                    i += 1
                elif ((seq[8] == "A" or seq[8] == "G")
                      and seq[9:15] == pJa[9:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_2Pu(inputf, inputfmt):
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
                elif seq[8] == "A" or seq[8] == "G":
                    i += 1
                elif (seq[7] == pJa[7] and seq[9:15] == pJa[9:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_3A(inputf, inputfmt):
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
                elif seq[9] == pJa[9]:
                    i += 1
                elif (seq[7] == pJa[7] and (seq[8] == "A" or seq[8] == "G")
                      and seq[10:15] == pJa[10:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_4A(inputf, inputfmt):
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
                elif seq[10] == pJa[10]:
                    i += 1
                elif (seq[7] == pJa[7] and (seq[8] == "A" or seq[8] == "G")
                      and seq[9] == pJa[9] and seq[11:15] == pJa[11:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_5A(inputf, inputfmt):
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
                elif seq[11] == pJa[11]:
                    i += 1
                elif (seq[7] == pJa[7] and (seq[8] == "A" or seq[8] == "G")
                      and seq[9:11] == pJa[9:11] and seq[12:15] == pJa[12:15]):
                    i += 1
                else:
                    i += 1


def check_pJa_6C(inputf, inputfmt):
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
                elif seq[12] == pJa[12]:
                    i += 1
                elif (seq[1:5] == pJa[1:5] and seq[9] == pJa[9] and
                      seq[13:16] == pJa[13:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_7G(inputf, inputfmt):
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
                elif seq[13] == pJa[13]:
                    i += 1
                elif (seq[1:5] == pJa[1:5] and seq[9] == pJa[9] and
                      seq[12] == pJa[12] and seq[14:16] == pJa[14:16]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_8G(inputf, inputfmt):
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
                elif seq[14] == pJa[14]:
                    i += 1
                elif (seq[1:5] == pJa[1:5] and seq[9] == pJa[9] and
                      seq[12:14] == pJa[12:14] and seq[15] == pJa[15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def check_pJa_9G(inputf, inputfmt):
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
                elif seq[15] == pJa[15]:
                    i += 1
                elif (seq[1:5] == pJa[1:5] and seq[9] == pJa[9] and
                      seq[12:15] == pJa[12:15]):
                    yield [seq_record.id, i, seq, seq_record.description]
                    i += 1
                else:
                    i += 1


def main():
    inputf = "/users/hamazaki/db/blast/human_ASs_2017.fa"
    print(inputf, "pJa intact")
    search_pJa = check_pJa(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_1")
    search_pJa_1T = check_pJa_1T(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_1T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_2T")
    search_pJa_2T = check_pJa_2T(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_2T)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_3C")
    search_pJa_3C = check_pJa_3C(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_3C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_4G")
    search_pJa_4G = check_pJa_4G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_4G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_5A")
    search_pJa_5A = check_pJa_5A(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_5A)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_6C")
    search_pJa_6C = check_pJa_6C(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_6C)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_7G")
    search_pJa_7G = check_pJa_7G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_7G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_8G")
    search_pJa_8G = check_pJa_8G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_8G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])
    print("pJa_9G")
    search_pJa_9G = check_pJa_9G(inputf, "fasta")
    while True:
        try:
            found = next(search_pJa_9G)
        except StopIteration:
            break
        else:
            print(found[3], found[1], found[2])


if __name__ == '__main__':
    main()
