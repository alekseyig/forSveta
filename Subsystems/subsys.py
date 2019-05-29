#!/usr/bin/env python

import os
import sys
from colorutils.convert import *

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

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

NK_DIFF = 50000


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

        self.wb = xlsxwriter.Workbook('/Users/alekseyzhukov/b.xlsx')
        self.wsheet = self.wb.add_worksheet()

    def close(self):
        self.wb.close()

    def collect_genes(self):
        self.genes = []  # type: list[Gene]
        for dirpath, _, filenames in os.walk(self.input_dir):
            for i, fname in enumerate(filenames):
                path = os.path.join(dirpath, fname)

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
                        for i, (annot,_) in enumerate(self.annotations):
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
        h %= 360
        if h == 0:
            h = 20
        else:
            h += (h/100. * 20)
        return hsv_to_hex((h, s, v)).upper()

    def _colorize_genes(self):
        genomes = set([g.genome_name for g in self.genes])

        # setting bg color
        current_bg_color = None
        seen_genes = []
        for genome in genomes:
            genes = [g for g in self.genes if g.genome_name == genome]
            genes.sort(key=lambda x: int(x.id.split('.')[-1]))

            for gene in genes:
                seen_gene = None
                for sg in seen_genes:
                    if sg.id == gene.id:
                        seen_gene = sg

                if seen_gene is None:
                    seen_genes.append(gene)
                else:
                    bg_color = seen_gene.bg_color
                    if bg_color == WHITE_BG and current_bg_color is None:
                        seen_gene.bg_color = BLUE_BG
                        gene.bg_color = BLUE_BG
                        current_bg_color = BLUE_BG
                    elif bg_color == WHITE_BG and current_bg_color is not None:
                        current_bg_color = self.__next_color(current_bg_color)
                        gene.bg_color = current_bg_color
                        seen_gene.bg_color = current_bg_color
                    else:
                        gene.bg_color = bg_color

        # setting fg color
        current_fg_color = RED_FG
        for genome in genomes:
            genes = [g for g in self.genes if g.genome_name == genome]
            genes.sort(key=lambda x: int(x.location.start))
            groups = []
            last_end = None
            for g in genes:
                if last_end is None:
                    last_end = g.location.end
                    groups.append([g])
                    continue

                if abs(g.location.start - last_end) <= NK_DIFF:
                    groups[-1].append(g)
                else:
                    groups.append([g])

            for group in groups:
                if len(group) == 1:
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
        self.annotations = []  # type: list[(str,str)]
        with open(self.input_file) as f:
            for ln in f:
                if not len(ln):
                    continue
                try:
                    short_name, annotation = ln.strip().split('\t')
                except Exception as e:
                    print "Please check the format of your Annotations input file, some lines are mis-formatted"
                    sys.exit(1)
                self.annotations.append((annotation, short_name))


if __name__ == '__main__':
    input_file, input_dir = sys.argv[1:]
    subsystems = Subsystems(input_file, input_dir)
    subsystems.process()
    subsystems.close()
