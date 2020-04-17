#!/usr/bin/env python
import os
import sys
import csv
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

import warnings
warnings.filterwarnings("ignore")


class ORF(object):
    def __init__(self):
        self._key = None
        self.parts = None
        self.strand = None

        # self.genome_id = None
        # self.genome_name = None
        # self.location = None
        # self.annotation = None
        # self.contig = None
        # self.translation = None


class Organism:
    def __init__(self):
        self._key = None
        self.name = None
        self.taxonomy = None

def main(args):
    path = args[0]
    annotations = []
    genes = []

    seen_org_name = False

    orf_id = 1
    with open('org.jsonl', 'w') as orgf, open("orfs.jsonl", 'w') as orff:
        for seq_record in SeqIO.parse(path, "genbank"):  # type: SeqRecord
            genome_id = None
            contig = seq_record.name
            # seen_org_name = seq_record.annotations['organism']

            if not seen_org_name:
                o = Organism()
                o._key = seq_record.id
                o.name = seq_record.annotations['organism']
                o.taxonomy = seq_record.annotations['taxonomy']
                id = seq_record.id
                orgf.write(json.dumps(o.__dict__))
                seen_org_name = True


            for feature in seq_record.features:  # type: SeqFeature
                if feature.type == 'source' or feature.type != 'CDS':
                    continue

                orf = ORF()
                orf._key = f"{id}:{orf_id}"
                orf.parts = [(fl.start, fl.end) for fl in feature.location.parts]
                orf.strand = feature.location.strand
                orff.write(json.dumps(orf.__dict__))
                orff.write('\n')

                orf_id += 1
                continue
    return



    #             gene_annot = feature.qualifiers['product'][0]  # type: str
    #             for (annot, _) in annotations:  # type: str
    #                 if annot in gene_annot:
    #                     gene = ORF()
    #                     gene.genome_id = genome_id
    #                     gene.contig = contig
    #                     gene.genome_name = seen_org_name
    #                     gene.id = [_id.lstrip("SEED:")
    #                                for _id in feature.qualifiers['db_xref'] if _id.startswith("SEED")][0]
    #                     gene.location = feature.location
    #                     gene.translation = feature.qualifiers['translation'][0]
    #                     gene.annotation = annot
    #                     genes.append(gene)
    #
    # print(">>> Genomes", len(genes))


if __name__ == "__main__":
    args = sys.argv[1:]
    sys.exit(main(args))
