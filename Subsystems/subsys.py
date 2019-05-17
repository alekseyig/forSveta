#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from xlwt import Workbook, Worksheet, Formula, easyxf

import warnings
warnings.filterwarnings("ignore")

RUST_JOB_ID = 'RUST job ID'
GENOME_ID = 'Genome ID'
ORG_NAME = 'Organism name'


def get_annotations(fname):
    result = []
    with open(fname) as f:
        for ln in f:
            short_name, annotation = ln.strip().split('\t')
            result.append((annotation, short_name))
    return result


def write_header(sheet, annotations):
    sheet.write(0, 0, RUST_JOB_ID)
    sheet.write(0, 1, GENOME_ID)
    sheet.write(0, 2, ORG_NAME)

    i = 2
    for _, abbr in annotations:
        sheet.write(0, i+1, '{a}, gen ID'.format(a=abbr))
        sheet.write(0, i+2, '{a}, size (aa)'.format(a=abbr))
        i += 2


def get_record(path, annotations):
    record = {}
    record[RUST_JOB_ID] = os.path.split(path)[-1].split('_')[0]

    for seq_record in SeqIO.parse(path, "genbank"):  # type: SeqRecord
        record[ORG_NAME] = seq_record.annotations['organism']

        for feature in seq_record.features:  # type: SeqFeature
            if feature.type == 'source':
                record[GENOME_ID] = feature.qualifiers['genome_id'][0]
                continue

            if feature.type != 'CDS':
                continue

            gene_annot = feature.qualifiers['product'][0]  # type: str
            for i, (anot,_) in enumerate(annotations):
                if anot in gene_annot:
                    gene_id = [_id for _id in feature.qualifiers['db_xref'] if _id.startswith("SEED")][0] \
                        .lstrip("SEED:").split('.')[-1]
                    # todo: need a hyperlink
                    #print gene_id
                    gene_ids = record.get(i, ([], []))[0]
                    gene_ids.append(gene_id)

                    gene_lengths = record.get(i, ([], []))[1]
                    gene_lengths.append(str(len(feature.qualifiers['translation'][0])))
                    record[i] = (gene_ids, gene_lengths)
    return record


def write_record(sheet,  # type: Worksheet
                 record, annot_count, row_num):
    sheet.write(row_num, 0, record[RUST_JOB_ID])
    sheet.write(row_num, 1, record[GENOME_ID])
    sheet.write(row_num, 2, record[ORG_NAME])

    pos = 3
    for i in xrange(annot_count):
        gene_ids, gene_lengths = record.get(i, ('', ''))
        f = Formula('HYPERLINK("http://stackoverflow.com/"; "SO3"), HYPERLINK("http://stackoverflow.com/"; "SO2")')
        #sheet.write(row_num, pos, ', '.join(gene_ids))
        sheet.write(row_num, pos, f)
        #print f.get_references()
        #print f.patch_references
        #print f.rpn()
        #print f.text()
        sheet.write(row_num, pos + 1, ', '.join(gene_lengths))
        pos += 2


def main(args):
    input_file, input_dir = args
    annotations = get_annotations(input_file)

    wb = Workbook()
    sheet = wb.add_sheet('Sheet 1')  # type: Worksheet

    write_header(sheet, annotations)

    for dirpath, _, filenames in os.walk(input_dir):
        for i, fname in enumerate(filenames):
            path = os.path.join(dirpath, fname)
            write_record(sheet, get_record(path, annotations), len(annotations), i+1)

    wb.save("/Users/alekseyzhukov/a.xls")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
