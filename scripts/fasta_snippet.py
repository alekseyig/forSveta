#!/usr/bin/env python
import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(params):
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--length", type=int, help="snippet length")
    parser.add_argument("fasta_file", action="store")
    args = parser.parse_args(params)

    fasta_file = args.fasta_file
    n = args.length
    if not os.path.exists(fasta_file):
        print(f"Could not find {fasta_file}")

    for record in SeqIO.parse(fasta_file, "fasta"):
        id = "|".join(record.description.split("|")[-3:])
        seq = str(record.seq)
        length = len(seq)

        if length <= n:
            id = "|".join((id, "ORF1ab", f"1-{length}"))
            sr = SeqRecord(Seq(data=seq), id=id, description='')
            SeqIO.write(sr, sys.stdout, "fasta")
        else:
            for i in range(0, length - n + 1):
                _id = "|".join((id, "ORF1ab", f"{i+1}-{i+n+1}"))
                sr = SeqRecord(Seq(data=seq[i:i+n]), id=_id, description='')
                SeqIO.write(sr, sys.stdout, "fasta")
        #break


if __name__ == '__main__':
    main(sys.argv[1:])
