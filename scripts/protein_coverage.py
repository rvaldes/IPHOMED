### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Summarizes protein callability

#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np

data = defaultdict(lambda: [0])

with open(sys.argv[1], "r") as file:
    for line in file:
        col = line.strip().split()
        for i in range(int(col[1]) + 1):
            data[col[0]].append(0)

while True:
    try:
        line = input()
        col = line.strip().split("\t")
        if col[1] in data:
            data[col[1]][0] += 1
            for i in range(int(col[8]), int(col[9]) + 1):
                data[col[1]][i] += 1
    except EOFError:
        break

for protein, values in data.items():
    total = values.pop(0)
    length = len(values)
    q0 = np.quantile(values, 0)
    q1 = np.quantile(values, 0.25)
    q2 = np.quantile(values, 0.5)
    q3 = np.quantile(values, 0.75)
    q4 = np.quantile(values, 1)
    print(f"{protein}\t{total}\t{length}\t{q0}\t{q1}\t{q2}\t{q3}\t{q4}")