#!/usr/bin/env python

import os
import sys
from colorutils.convert import *

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from typing import List

import xlsxwriter

import warnings
warnings.filterwarnings("ignore")

RAST_JOB_ID = 'RAST job ID'
GENOME_ID = 'Genome ID'
GENOME_NAME = 'Organism name'

BLUE_BG = '#ABD8F2'
WHITE_BG = '#FFFFFF'
BLACK_FG = '#000000'
RED_FG = '#FF0000'

NK_DIFF = 5000


class Gene(object):
    def __init__(self):
        self.id = None
        self.genome_id = None
        self.genome_name = None
        self.location = None
        self.annotation = None
        self.rast_job_id = None
        self.contig = None
        self.translation = None
        self.bg_color = WHITE_BG
        self.fg_color = BLACK_FG


class Subsystems(object):
    def __init__(self, input_file, input_dir):
        self.input_file = input_file
        self.input_dir = input_dir

        self.genes = []  # type: List[Gene]
        self.annotations = []  # type: List[(str, str)]

        self.wb = xlsxwriter.Workbook(os.path.join(os.path.expanduser("~"), 'raw_sybsystems.xlsx'))
        self.wsheet = self.wb.add_worksheet()

    def close(self):
        self.wb.close()

    def collect_genes(self):
        for dirpath, _, filenames in os.walk(self.input_dir):
            filenames.sort()
            for fname in filenames:
                path = os.path.join(dirpath, fname)
                print "Processing genome:", fname

                rast_job_id = os.path.split(path)[-1].split('_')[0] if "_" in path else ""

                for seq_record in SeqIO.parse(path, "genbank"):  # type: SeqRecord
                    genome_id = None
                    contig = seq_record.name
                    org_name = seq_record.annotations['organism']

                    for feature in seq_record.features:  # type: SeqFeature
                        if feature.type == 'source':
                            genome_id = feature.qualifiers['genome_id'][0]
                            continue

                        if feature.type != 'CDS':
                            continue

                        gene_annot = feature.qualifiers['product'][0]  # type: str
                        for (annot, _) in self.annotations:  # type: str
                            if annot in gene_annot:
                                gene = Gene()
                                gene.genome_id = genome_id
                                gene.contig = contig
                                gene.genome_name = org_name
                                gene.rast_job_id = rast_job_id
                                gene.id = [_id.lstrip("SEED:")
                                           for _id in feature.qualifiers['db_xref'] if _id.startswith("SEED")][0]
                                gene.location = feature.location
                                gene.translation = feature.qualifiers['translation'][0]
                                gene.annotation = annot
                                self.genes.append(gene)
        self._colorize_genes()

    def __next_color(self, hex_color):
        h, s, v = hex_to_hsv(hex_color)
        if h == 0:
            h = 20
        else:
            h += (h/100. * 20)
        h %= 360
        return hsv_to_hex((h, s, v)).upper()

    def _colorize_genes(self):
        genomes = set([g.genome_name for g in self.genes])

        # setting bg color
        for genome in genomes:
            current_bg_color = None
            last_seen_gene = None

            genes = [g for g in self.genes if g.genome_name == genome]
            genes.sort(key=lambda x: int(x.id.split('.')[-1]))

            for gene in genes:
                if last_seen_gene is None:
                    last_seen_gene = gene
                    continue

                if gene.id == last_seen_gene.id:
                    bg_color = last_seen_gene.bg_color
                    if bg_color == WHITE_BG and current_bg_color is None:
                        last_seen_gene.bg_color = BLUE_BG
                        gene.bg_color = BLUE_BG
                        current_bg_color = BLUE_BG
                    elif bg_color == WHITE_BG and current_bg_color is not None:
                        current_bg_color = self.__next_color(current_bg_color)
                        gene.bg_color = current_bg_color
                        last_seen_gene.bg_color = current_bg_color
                    else:
                        gene.bg_color = bg_color
                last_seen_gene = gene

        # setting fg color
        for genome in genomes:
            current_fg_color = RED_FG
            groups = []

            genes = [g for g in self.genes if g.genome_name == genome]
            genes.sort(key=lambda x: int(x.location.start))

            last_end = None
            last_start = None
            for g in genes:
                if last_end is None:
                    last_end = g.location.end
                    last_start = g.location.start
                    groups.append({g})
                    continue

                if abs(g.location.start - last_end) <= NK_DIFF or\
                        (last_start <= g.location.end and g.location.start <= last_end):
                    groups[-1].add(g)
                else:
                    groups.append({g})
                last_start, last_end = g.location.start, g.location.end

            for group in groups:
                if len(set([g.id for g in group])) == 1:
                    continue
                else:
                    for g in group:
                        g.fg_color = current_fg_color
                    current_fg_color = self.__next_color(current_fg_color)

    def process(self):
        self.get_annotations()
        self.collect_genes()
        self.write_genes()

    def write_genes(self):
        black_url = self.wb.add_format({'font_color': 'black', 'underline':  1})
        bold = self.wb.add_format({'bold': True})

        genomes = set([g.genome_name for g in self.genes])

        # write rast_job_id column
        row, col = 0, 0
        self.wsheet.write(row, col, RAST_JOB_ID, bold)
        row += 1
        for genome in genomes:
            genes = [g for g in self.genes if g.genome_name == genome]
            self.wsheet.write(row, col, genes[0].rast_job_id)
            row += 1

        # write genome_id column
        row, col = 0, 1
        self.wsheet.write(row, col, GENOME_ID, bold)
        row += 1
        for genome in genomes:
            genes = [g for g in self.genes if g.genome_name == genome]
            self.wsheet.write_url(row, col,
                                  'http://rast.theseed.org/FIG/seedviewer.cgi?page=Organism&organism={genome_id}'
                                  .format(genome_id=genes[0].genome_id), string=genes[0].genome_id, cell_format=black_url)
            row += 1

        # write genome_name column
        row, col = 0, 2
        self.wsheet.write(row, col, GENOME_NAME, bold)
        row += 1
        for genome in genomes:
            genes = [g for g in self.genes if g.genome_name == genome]
            self.wsheet.write(row, col, genes[0].genome_name)
            row += 1

        # write annotation columns
        col = 3
        for annot, short_name in self.annotations:
            genes = [g for g in self.genes if annot in g.annotation]
            col_num = 0
            genes_by_genome = {}
            for genome in genomes:
                _genes = [g for g in genes if g.genome_name == genome]
                col_num = max(col_num, len(_genes))
                genes_by_genome[genome] = _genes

            for i in xrange(col_num):
                row = 0
                self.wsheet.write(row, col, short_name, bold)
                row += 1
                for genome in genomes:
                    genes = genes_by_genome[genome]
                    gene = genes.pop(0) if len(genes) else None  # type: Gene
                    if gene is None:
                        self.wsheet.write(row, col, "")
                    else:
                        cell_format = self.wb.add_format({'font_color': gene.fg_color, 'underline':  1, 'bg_color': gene.bg_color})
                        text = gene.id.split('.')[-1] + '(' + str(len(gene.translation)) + ')'
                        self.wsheet.write_url(row, col,
                                              "http://rast.theseed.org/FIG/seedviewer.cgi?page=Annotation&feature={fig_id}"
                                              .format(fig_id=gene.id), string=text, cell_format=cell_format)
                    row += 1
                col += 1

    def get_annotations(self):
        with open(self.input_file) as f:
            for ln in f:
                if not len(ln):
                    continue
                try:
                    short_name, annotation = ln.strip().split('\t')  # type: str
                except Exception:
                    print "Please check the format of your Annotations input file, some lines are mis-formatted"
                    sys.exit(1)
                self.annotations.append((annotation, short_name))


if __name__ == '__main__':
    input_file, input_dir = sys.argv[1:]
    subsystems = Subsystems(input_file, input_dir)
    subsystems.process()
    subsystems.close()
