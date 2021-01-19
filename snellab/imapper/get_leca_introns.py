#!/usr/bin/env python

import sys

introns = {}
with open(sys.argv[1]) as intron_table:
    for line in intron_table:
        if line[:4] != '#MAP':
            continue
        fields = line.rstrip().split('\t')
        og = fields[2].split('/')[1]
        position = int(fields[3])
        phase = position % 3
        if phase == 0:
            location = f"{position}'3"
        else:
            location = f"{position + 1}'{phase}"
        introns[fields[1]] = {'Family': og, 'Location': location}

with open(sys.argv[2]) as site_histories:
    for i in range(3):
        site_histories.readline()
    check = site_histories.readline().split('\t')[317]
    if check != 'LECA-present':
        sys.exit('Check if file is correct.')
    for line in site_histories:
        fields = line.rstrip().split('\t')
        site = fields[0]
        leca_pr = float(fields[317])
        if leca_pr >= 0.5:
            info = introns[site]
            print(f"{site}\t{info['Family']}\t{info['Location']}\t{round(leca_pr, 1)}")
