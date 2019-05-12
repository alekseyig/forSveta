#!/usr/bin/env python

from __future__ import with_statement
import os

blast_file = os.path.join("supportive_data", "PubSEED-ridA-Hits.txt")

orgs = {} # 7097
with open(blast_file) as f:
    for ln in f:
        fields = ln.strip().split('\t')
        try:
            if '.' in fields[1]:
                orgs[fields[1]] = None
        except IndexError:
            pass

for org in orgs:
  print org
#print len(orgs)
