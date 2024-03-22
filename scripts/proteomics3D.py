### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Divide proteomics output into host, bacterial and dietary proteins

#!/usr/bin/env python

import re
import sys

# host database
host = {}
with open(sys.argv[2], 'r') as file:
    for line in file:
        line = line.strip()
        if '>' in line:
            match = re.search(r'\|([0-9A-Z]*)\|', line)
            if match:
                protein = match.group(1)
                host[protein] = 0

# bacterial database
bacteria = {}
with open(sys.argv[3], 'r') as file:
    for line in file:
        line = line.strip()
        if '>' in line:
            line = line.replace('>', '')
            bacteria[line] = 0

# dietary database
diet = {}
with open(sys.argv[4], 'r') as file:
    for line in file:
        line = line.strip()
        if '>' in line:
            line = line.replace('>', '')
            diet[line] = 0

# preparing output files
host_file = open(sys.argv[5], 'w')
bacteria_file = open(sys.argv[6], 'w')
diet_file = open(sys.argv[7], 'w')

# proteomics output
with open(sys.argv[1], 'r') as file:
    for line in file:
        line = line.strip()
        cols = line.split('\t')
        if cols[0] == "Protein Accession":
            host_file.write(line + '\n')
            bacteria_file.write(line + '\n')
            diet_file.write(line + '\n')
        else:
            if float(cols[5]) > 1:
                continue
            if cols[0] in host:
                host_file.write(line + '\n')
            if cols[0] in bacteria:
                bacteria_file.write(line + '\n')
            if cols[0] in diet:
                diet_file.write(line + '\n')

host_file.close()
bacteria_file.close()
diet_file.close()
