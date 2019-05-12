#!/usr/bin/env python
from __future__ import with_statement

import os, sys
# urllib2, logging, time, urllib, socket
#from argparse import ArgumentParser
#from xml.dom.minidom import parse

sys.path.append(os.path.join(os.path.dirname(__file__), '..',
                '../lib/biopython-1.62/build/lib.macosx-10.6-intel-2.7'))
#import html5lib
#from bs4 import BeautifulSoup
from Bio.Blast import NCBIStandalone


functions = {}
with open('../Creinhardtii_236_defline.txt') as f:
  for ln in f:
    id, func = ln.strip().split('\t')
    functions[id] = func

seed_functions = {}
with open('../Chlamy_SEED_ids.functions') as f:
  for ln in f:
    id, func = ln.strip().split('\t')
    seed_functions[id] = func

seed_subsys = {}
with open('../Chlamy_SEED_ids.subsys') as f:
  for ln in f:
    id, func = ln.strip().split('\t')
    seed_subsys.setdefault(id,[]).append(func)




print 'Query ID\tQuery Length\tHit ID\tHit Length\tE-value\tQuery Identity\tHit Identity\tPercent\tHit Function\tQuery Function\tQuery Subsys'

blast_parser = NCBIStandalone.BlastParser()

for dirname,_,fnames in os.walk('../blast_results'):
  fnames = [os.path.join(dirname,fname) for fname in fnames if fname.startswith('a')]
  for fname in fnames:
    #print "Process file: " + fname
    result_handle = open(fname)
    blast_iterator = NCBIStandalone.Iterator(result_handle, blast_parser)

    for blast_record in blast_iterator:
      query = blast_record.query
      query_length = blast_record.query_letters

      hit, hit_length, e_val = '', '', ''
      q_identity, h_identity = '', ''
      percent = ''
      if len(blast_record.alignments):
        alignment = blast_record.alignments[0]
        title = alignment.title
        desc = [desc for desc in title.split(' ') if desc.startswith('transcriptName_')][0]
        hit = desc.split('transcriptName_')[1]

        hit_length = alignment.length

        hsp = alignment.hsps[0]
        e_val = "%.2E" % float(hsp.expect)

        q_identity, h_identity = hsp.identities
        percent = "%d" % round(int(q_identity) * 1.0 / int(h_identity) * 100)

      hit_function = functions.get(hit,'')
      seed_function = seed_functions.get(query,'')

      sybsys = seed_subsys.get(query)
      if sybsys is None:
        sybsys = ''
      else:
        sybsys = '\t'.join(sybsys)

      result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (query, query_length, hit, hit_length, e_val, q_identity, h_identity, percent, hit_function, seed_function, sybsys)
      print result


      #print ">>> BlastRecord"
      #for d in dir(blast_record):
      #  if d.startswith('_'):
      #    continue
      #  print '!!!B', d, getattr(blast_record,d)

      #for descriptions in blast_record.descriptions:
      #  print ">>> Desc "
      #  for d in dir(descriptions):
      #    if d.startswith('_'):
      #      continue
      #    print '!!!D', d, getattr(descriptions,d)

      #for alignment in blast_record.alignments:
      #  print ">>> Alignment "
      #  for d in dir(alignment):
      #    if d.startswith('_'):
      #      continue
      #    print '!!!A', d, getattr(alignment,d)

      #  for hsp in alignment.hsps:
      #    print ">>> HSP "
      #    for d in dir(hsp):
      #      if d.startswith('_'):
      #        continue
#      print '!!!H', d, getattr(hsp,d)



      #for alignment in blast_record.alignments:
      #  for hsp in alignment.hsps:
      #    print '****Alignment****'
      #    print 'sequence:', alignment.title
      #    print 'length:', alignment.length
      #    print 'e value:', hsp.expect
      #    if len(hsp.query) > 75:
      #      dots = '...'
      #    else:
      #      dots = ''
      #      print hsp.query[0:75] + dots
      #      print hsp.match[0:75] + dots
      #      print hsp.sbjct[0:75] + dots



    result_handle.close()
    #break

