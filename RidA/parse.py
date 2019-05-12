#!/usr/bin/env python2.7
from __future__ import with_statement
from Bio import SwissProt

#f = open('O82662') 
f = open('uniprot_sprot.dat')
for r in SwissProt.parse(f):
    if r.organism == 'Arabidopsis thaliana (Mouse-ear cress).':
        #if 'Arabidopsis thaliana (Mouse-ear cress)' in r.organism and r.entry_name == 'PUB33_ARATH':
        #print dir(r) #.description 
        #print r.__dict__

        #link = '=HYPERLINK("http://www.uniprot.org/uniprot/%s","%s")' %(r.accessions[0],r.accessions[0])
        link = '=HYPERLINK("http://www.uniprot.org/uniprot/%s","%s")' %(r.entry_name, r.entry_name)

        if 'OrderedLocusNames=' in r.gene_name:
            id = r.gene_name.split('OrderedLocusNames=')[1].split(';')[0].upper()
        else:
            id=''
            # CWP02_ARATH
            #PUB33_ARATH

        description_fileds = r.description.split(';')
        name = description_fileds[0]
        name = name.strip().split('RecName: Full=')[1]

        ec = []
        for field in description_fileds[1:]:
            if 'EC=' in field:
                ec.append(field.strip().split('EC=')[1].strip())
        if ec:
            ec = ' '.join(ec)
        else:
            ec= ''
                          
        compound, c_activity = '',''
        for comment in r.comments:
            if comment.startswith('SUBCELLULAR LOCATION:'):
                compound = comment.split('SUBCELLULAR LOCATION:')[1].strip().rstrip('.')
            if comment.startswith('CATALYTIC ACTIVITY:'):
                c_activity = comment.split('CATALYTIC ACTIVITY:')[1].strip()
        #result = [id, link, r.entry_name, name, ec, compound, c_activity]
        result = [id, link, name, ec, compound, c_activity]
        print '\t'.join(result)


f.close()

#with open('uniprot-Arabidopsis+thaliana.tab') as f:
#    for ln in f:
#        entry, entry_name, status, protein_names, gene_names, organism, length = ln.strip().split('\t')
#        #O82662	SUCB_ARATH	reviewed	Succinyl-CoA ligase [ADP-forming] subunit beta, mitochondrial (EC 6.2.1.5) (Succinyl-CoA synthetase beta chain) (SCS-beta)	At2g20420 F11A3.3	Arabidopsis thaliana (Mouse-ear cress)	421
#        if entry == 'O82662':
#            recomended_name, ec = protein_names.split('(EC ')
#            ec = ec.split(')')[0]
#            result = [gene_names.upper().split(' ')[0], recomended_name, ec]
#            print '\t'.join(result)
