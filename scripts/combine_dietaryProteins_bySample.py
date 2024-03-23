### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Combine dietary proteins from individual sample detection

#!/usr/bin/env python

import os
import glob
import pandas as pd
import sys

# Pattern to match files
pattern = 'iphomed/Search/Task1SearchTask/*/*-calib_ProteinGroups.tsv'

# Find files matching the pattern
files = glob.glob(pattern)
files.sort()


# Combine files into a single DataFrame
dfs = []
for file in files[0:]:
    file_rep = file.replace(" ", "\ ")
    cmd = "sed -i 's/\t$//' %s"%(file_rep)
    os.system(cmd)
    df = pd.read_csv(file, sep='\t')
    df['Sample'] = os.path.basename(file).replace('-calib_ProteinGroups.tsv', '')
    df.rename(columns={df.columns[15]: 'Intensity'}, inplace=True)
    dfs.append(df)

df_combined = pd.concat(dfs)

columns_to_remove = ['Fragment Sequence Coverage', 'Modification Info List',
                     'Sequence Coverage with Mods', 'Sequence Coverage']
df_combined.drop(columns=columns_to_remove, inplace=True)

# dietary peptides
diet = {}
with open(sys.argv[1], 'r') as file:
    for line in file:
        line = line.strip()
        cols = line.split('\t')
        diet[cols[0]] = 0

diet_set = set(diet.keys())

# Filter the df_combined DataFrame
df_filtered = df_combined[df_combined.iloc[:, 0].isin(diet_set)]

# Save the combined data to a single file
df_filtered.to_csv(sys.argv[2], index=False, sep='\t')
