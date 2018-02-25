import csv
from domain_CENPB import check_CENPB, check_CENPBrevcom


for i in range(1, 5):
    inputf = "/users/hamazaki/db/blast/D1Z5_peri_5end_"+str(i)+".fa"
    out_f = "CENPBf_9GtoX_D1Z5_peri"+str(i)+".tsv"
    out_r = "CENPBr_9GtoX_D1Z5_peri"+str(i)+".tsv"

    print("forward"+str(i))

    with open(out_f, "w") as outfile:
        search_CENPB = check_CENPB(inputf, "fasta")
        while True:
            try:
                found = next(search_CENPB)
            except StopIteration:
                break
            else:
                row = [found[0], found[1], found[1]+17, found[2]]
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(row)

    print("reverse"+str(i))

    with open(out_r, "w") as outfile:
        search_CENPBrc = check_CENPBrevcom(inputf, "fasta")
        while True:
            try:
                foundrc = next(search_CENPBrc)
            except StopIteration:
                break
            else:
                rowrc = [foundrc[0], foundrc[1], foundrc[1]+17, foundrc[2]]
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(rowrc)
