#!/usr/bin/env python
from __future__ import with_statement
from Bio import Entrez
import os

pegs_to_ncbi_file = os.path.join("supportive_data", "pegs_to_ncbi.short")
#pegs_to_ncbi_file = os.path.join("supportive_data", "pegs_to_ncbi")

families_file = os.path.join("supportive_data", "cdd_families.txt")
super_families_file = os.path.join("supportive_data", "cdd_super_families.txt")
other_families_file = os.path.join("supportive_data", "cdd_other_families.txt")


Entrez.tool = __file__
Entrez.email = "alekseyig@hotmail.com"


aliases = {}
with open(pegs_to_ncbi_file) as f:
  for ln in f:
    ln = ln.strip()
    try:
      peg, ncbi_ids = ln.split('\t')
      ncbi_ids = ncbi_ids.split(',')
      ncbi_ids.sort()
      alias_id = ncbi_ids[0]
      alias_id = alias_id.split('|')[1]
      aliases.setdefault(alias_id, []).append(peg)
    except ValueError:
      print "%s\tNO ALIAS" %(ln)
    except IndexError:
      print ln
      raise IndexError


families_file_handler = open(families_file, 'w')
super_families_file_handler = open(super_families_file, 'w')
other_families_file_handler = open(other_families_file, 'w')

fam_info = {}
for id in aliases:
  results = Entrez.read(Entrez.elink(dbfrom="protein", id=id, db='cdd'))[0]['LinkSetDb']
  if len(results):
    fam_ids = {}
    #families = {}
    #super_families = {}
    for links_dic in results:
      if links_dic['LinkName'] == 'protein_cdd_concise_2':
        #families = dict(((d['Id'],None) for d in links_dic['Link']))
        fam_ids = dict(((d['Id'],None) for d in links_dic['Link']))
      if links_dic['LinkName'] == 'protein_cdd_superfamily_2':
        #super_families = dict(((d['Id'],None) for d in links_dic['Link']))
        fam_ids.update(dict(((d['Id'],None) for d in links_dic['Link'])))

    if not len(fam_ids):
      for peg in aliases[id]:
        print "%s\tgi|%s" %(peg,id)
      continue
    #print id
    #print results
    #print "Fam:", families
    #print "SupFam:", super_families
    #print "All fams:", fam_ids
    #print

    cls, cds, others = [], [], []
    for fam_id in fam_ids:
      if fam_id not in fam_info:
        result = Entrez.read( Entrez.esummary(db="cdd", id=fam_id))[0]
        fam_info[fam_id] = (result['Accession'], result['Title'])

      accession, title = fam_info[fam_id]
      if accession.startswith('cl'):
        cls.append((accession,title))
      elif accession.startswith('cd'):
        cds.append((accession,title))
      else:
        others.append((accession,title))

    for peg in aliases[id]:
      if len(cls):
        super_families_file_handler.write("%s\tgi|%s" % (peg,id))
        for cl, title in cls:
          super_families_file_handler.write("\t%s:%s" % (cl,title))
        super_families_file_handler.write("\n")

      if len(cds):
        families_file_handler.write("%s\tgi|%s" % (peg,id))
        for cd, title in cds:
          families_file_handler.write("\t%s:%s" % (cd,title))
        families_file_handler.write("\n")

      if len(others):
        other_families_file_handler.write("%s\tgi|%s" % (peg,id))
        for other, title in others:
          other_families_file_handler.write("\t%s:%s" % (other,title))
        other_families_file_handler.write("\n")

  else:
    #print id, "pusto"
    #print

    for peg in aliases[id]:
      print "%s\tgi|%s" %(peg,id)

families_file_handler.close()
super_families_file_handler.close()
other_families_file_handler.close()


