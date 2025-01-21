### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Filter dietary proteins based on peptide analysis

#!/usr/bin/env python

import re
import sys

# dietary peptides
diet = {}
with open(sys.argv[1], 'r') as file:
    for line in file:
        line = line.strip()
        cols = line.split('\t')
        if cols[0] == "Query":
            continue
        if int(cols[1]) > 0 and float(cols[3]) >= 100:
            diet[cols[4]] = 0

# preparing output file
diet_file = open(sys.argv[3], 'w')

# dietary proteomics output
with open(sys.argv[2], 'r') as file:
    for line in file:
        line = line.strip()
        cols = line.split('\t')
        if cols[0] == "Protein Accession":
            diet_file.write(line + '\n')
        else:
            if cols[0] in diet and int(cols[9]) > 0:
                diet_file.write(line + '\n')

diet_file.close()
