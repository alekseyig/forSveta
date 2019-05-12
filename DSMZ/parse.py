#!/usr/bin/env python
import sys
import requests
from bs4 import BeautifulSoup
import csv

# Name:  Atopobium minutum,  before (
# DSM No.:   DSM-20349
# Strain designation:   everything
# Aliases   ATCC 13162, IFO 3516, NBRC 3516, NCIB 10784
# Isolated from:     flowers
# Country: country of origin unknown
# Date of sampling: before 22.08.1990
# Genbank accession numbers: whole genome shotgun sequence: JQBH00000000 (link is preferable)
#
# 'Name' 'DSM number' 'Strain designation' 'Aliases' 'Isolated from' 'Country Date o sampling' '16S rRNA' 'Genome sequence'

dt_names = {'Name: ', 'DSM No.: ', 'Strain designation: ', 'Other collection no. or WDCM no.: ',
            'Isolated from: ', 'Country: ', 'Date of sampling: ', 'Genbank accession numbers: '}

header = ['Name', 'DSM #', 'Strain designation', 'Aliases', 'Isolated from', 'Country', 'Date of sampling',
          'Genome sequence', '16S rRNA']


def main(argv):
    dsm_ids_file = argv[0]

    dsmz_ids = []
    with open(dsm_ids_file) as f:
        for ln in f:
            dsmz_ids.append(ln.strip())
    #dsmz_ids.append('DSM-1296')
    #dsmz_ids.append('DSM-29216')
    #dsmz_ids.append('DSM-20190')
    #dsmz_ids.append('DSM-7084')

    w = csv.DictWriter(open('output.csv', 'w'), fieldnames=list(header), delimiter='\t')
    w.writeheader()

    for dsmz_id in dsmz_ids:
        page = requests.get("https://www.dsmz.de/catalogues/details/culture/{id}.html".format(id=dsmz_id))
        soup = BeautifulSoup(page.text, 'html.parser')
        aa = soup.find(class_='tx-dsmzresources')
        print dsmz_id

        row = {}
        for dt_tag in aa('dt'):
            dd_tag = dt_tag.find_next('dd')

            dt = dt_tag.text
            dd = dd_tag.text

            if dt in dt_names:
                if dt == 'DSM No.: ':
                    dd = dd.split(',')[0]
                    row["DSM #"] = "DSM-" + dd.strip().encode('utf-8')

                if dt == 'Name: ':
                    dd = dd.split('(')[0]
                    row['Name'] = dd.strip().encode('utf-8')

                if dt == 'Strain designation: ':
                    row['Strain designation'] = dd.strip().encode('utf-8')

                if dt == 'Other collection no. or WDCM no.: ':
                    row['Aliases'] = dd.strip().encode('utf-8')

                if dt == 'Isolated from: ':
                    row['Isolated from'] = dd.strip().encode('utf-8')

                if dt == 'Country: ':
                    row['Country'] = dd.strip().encode('utf-8')

                if dt == 'Date of sampling: ':
                    row['Date of sampling'] = dd.strip().encode('utf-8')

                if dt == 'Genbank accession numbers: ':
                    if '16S rRNA gene: ' in dd:
                        rrnas = [rrna.strip().encode('utf-8') for rrna in dd.split('16S rRNA gene: ')[1:]]
                        row['16S rRNA'] = ','.join(rrnas)
                    elif 'whole genome shotgun sequence: ' in dd:
                        row['Genome sequence'] = dd.split('whole genome shotgun sequence: ')[1].strip().encode('utf-8')
                    else:
                        row['Genome sequence'] = dd.strip().encode('utf-8')

        if len(row):
            #print row
            w.writerow(row)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
