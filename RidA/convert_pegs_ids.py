#!/usr/bin/env python
from __future__ import with_statement
from Bio import SeqIO
import os

pegs_fasta_file = os.path.join("supportive_data", "pegs_not_found.fasta")
#pegs_fasta_file = os.path.join("supportive_data", "pegs_not_found.fasta.short")
pegs_fasta_converted_file = os.path.join("supportive_data", "pegs_not_found.fasta.converted")
pegs_to_new_id_file = os.path.join("supportive_data", "pegs_to_new_id")

pegs_fasta_converted_file_handler = open(pegs_fasta_converted_file, 'w')
pegs_to_new_id_file_handler = open(pegs_to_new_id_file, 'w')


count_id = 1
for seq_record in SeqIO.parse(pegs_fasta_file, "fasta"):
  pegs_fasta_converted_file_handler.write(">%s\n" % count_id)
  pegs_fasta_converted_file_handler.write("%s\n" % seq_record.seq)

  pegs_to_new_id_file_handler.write("%s\t%s\n" % (seq_record.id, count_id))
  count_id += 1


pegs_fasta_converted_file_handler.close()
pegs_to_new_id_file_handler.close()
