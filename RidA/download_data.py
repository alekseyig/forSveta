#!/usr/bin/env python

from __future__ import with_statement
import sys, subprocess, time, os

cmd_template = "/Applications/myRAST.app/bin/%s"
data_dir = "data"
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

#print len(orgs)

have_data_for = {}
for _,_,files in os.walk(data_dir):
    for fname in files:
        have_data_for[fname] = None
    
i = 0
for org in orgs:
    i += 1
    if org in have_data_for:
        continue
    cmd = cmd_template % 'svr_all_features %s peg | svr_location_of > data/%s'
    cmd = cmd % (org,org)
    t = time.localtime()
    sys.stdout.write("%d %d:%d:%d %s\n" % (i, t.tm_hour, t.tm_min, t.tm_sec, cmd))
    sys.stdout.flush()
    subprocess.check_output(cmd, shell=True)
    


