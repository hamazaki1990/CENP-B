from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Seq


def check_CENPB(inputf, inputfmt):
    CENPB = Seq.Seq("NTTCGNNNNANNCGGGN")
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
                elif (seq[1:5] == CENPB[1:5] and
                      seq[9] == CENPB[9] and seq[12:16] == CENPB[12:16]):
                    yield [seq_record.id, i, seq]
                    i += 1
                else:
                    i += 1


def check_CENPBrev(inputf, inputfmt):
    CENPB = Seq.Seq("NTTCGNNNNANNCGGGN")
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


def main():
    inputf = "/Users/hamazaki/db/blast/GenData2015_cons.fa"
    print("forward")
    search_CENPB = check_CENPB(inputf, "fasta")
    while True:
        try:
            found = next(search_CENPB)
        except StopIteration:
            break
        else:
            print(found)
    print("reverse")
    search_CENPBrev = check_CENPBrev(inputf, "fasta")
    while True:
        try:
            foundr = next(search_CENPBrev)
        except StopIteration:
            break
        else:
            print(foundr)


if __name__ == '__main__':
    main()
