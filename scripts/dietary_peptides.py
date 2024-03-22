### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Extract dietary peptides from detected dietary proteins

#!/usr/bin/env python

import re
import sys


with open(sys.argv[1], 'r') as file:
    for line in file:
        line = line.strip()
        cols = line.split('\t')
        if cols[0] == "Protein Accession":
            continue
        campos = cols[0].split()
        
        unique = cols[6].split('|')
        position = 0
        for i in range(len(unique)):
            if ( len(unique[i]) > 0):
                print(f">{campos[0]}_{position}\n{unique[i]}")
                position += 1
        
        shared = cols[7].split('|')
        for j in range(len(shared)):
            if ( len(shared[j]) > 0):
                print(f">{campos[0]}_{position}\n{shared[j]}")
                position += 1
