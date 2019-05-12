#!/usr/bin/env python
from __future__ import with_statement
import os

hits_file = os.path.join("supportive_data", "hitdata.txt")
pegs_to_new_id_file = os.path.join("supportive_data", "pegs_to_new_id")

#pegs_to_new_id_file_handler = open(pegs_to_new_id_file, 'w')

id_to_peg = {}
with open(pegs_to_new_id_file) as f:
  for ln in f:
    peg, id = ln.strip().split('\t')
    id_to_peg[id] = peg


#Q#1 - >1        cd00448 YjgF_YER057c_UK114_family
with open(hits_file) as f:
  for ln in f:
    id, cd, name = ln.strip().split('\t')
    id = id.split('>')[1]

    #if not (cd.startswith('cd') or cd.startswith('cl')):
    if cd.startswith('cl'):
      print "%s\t\t%s:%s" % (id_to_peg[id], cd, name)

