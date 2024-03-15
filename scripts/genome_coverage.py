### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Summarizes the genome coverage (callabilities)

#!/usr/bin/env python
import sys

convert = {}
with open(sys.argv[1], "r") as file:
    for line in file:
        col = line.strip().split("\t")
        convert[col[1]] = col[0]

count = {}
with open(sys.argv[2], "r") as file:
    for line in file:
        col = line.strip().split("\t")
        taxid = convert.get(col[0])
        count.setdefault(taxid, [0, 0])[0] += 1
        if int(col[2]) > 0:
            count.setdefault(taxid, [0, 0])[1] += 1

for tax, values in sorted(count.items()):
    print(f"{tax}\t{values[0]}\t{values[1]}")