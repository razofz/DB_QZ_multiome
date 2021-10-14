import scarf
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
from collections import Counter


args = str(sys.argv)
valid_sample_names = ['GFP_EPCR_LSKs_RNA', 'GFP_Viable_RNA']
if sys.argv[1] not in valid_sample_names:
    sys.exit(f'"{sys.argv[1]}" is not a valid sample name. Pick one of {valid_sample_names}.')
else:
    sample = sys.argv[1]

ds = scarf.DataStore(zarr_loc=f'data/processed/{sample}/counts_{sample}.zarr',
                     default_assay='RNA', nthreads=16)

ct = pd.read_csv('data/processed/GFP_EPCR_LSKs_RNA/celltag/v3.celltag.parsed.tsv',
            delimiter='\t')

cellsds = list(ds.cells.to_pandas_dataframe(columns=['ids'])['ids'])
len(cellsds)
cellsds = [x.split('-')[0] for x in cellsds]

cellsct = list(ct[['Cell.BC']]['Cell.BC'])
len(cellsct)
len(set(cellsct))
# nbr cells w celltag(s)
len(set(cellsds).intersection(set(cellsct)))
# nbr cells w/o celltag(s)
len(set(cellsds).difference(set(cellsct)))
# seems to be same as CT.celltag.stats.txt file

