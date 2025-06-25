#!/usr/bin/env python3

import sys
fasta = sys.argv[1]
recoding_scheme = sys.argv[2]
if recoding_scheme == 'SR4':
    recoding_table = str.maketrans('AGNPSTCHWYDEKQRFILMV', 'AAAAAACCCCGGGGGTTTTT')
elif recoding_scheme == 'SR6':
    recoding_table = str.maketrans('APSTCWDEGNFHYILMVKQR', 'AAAACCDDDDFFFIIIIKKK')
else:
    sys.exit(f'{recoding_scheme} not recognised. Only options are SR4 and SR6.')

with open(fasta) as fastafile:
    for line in fastafile:
        if line.startswith('>'):
            print(line, end = '')
        else:
            print(line.translate(recoding_table), end = '')
