# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Prepare our data for inputting into SAILERX
#
# ---

# %%
import datetime
import pytz

"Notebook last run at " + datetime.datetime.now(pytz.timezone("Europe/Stockholm")).strftime("%Y-%m-%d %H:%M")

# %%
import os
import re
from collections import Counter
from dataclasses import dataclass
from io import StringIO

import h5py
import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix
import collections

import scipy.sparse as sp
import tables

# %%
output_dir = os.getenv("PROJECT_PATH") + "/data/processed/notebooks/Immature/"
snakemake_analysis_path = f"{os.getenv('PROJECT_PATH')}/data/processed/snakemake/"
inputs = {
    "their_h5": f"{os.getenv('PROJECT_PATH')}/data/external/pbmc10k.hdf5",
    "cellranger_h5": f"{os.getenv('PROJECT_PATH')}/data/raw/aggregate/2021_083_Qinyu_agg/outs/filtered_feature_bc_matrix.h5",
    "cell_whitelist": f"{output_dir}cell_whitelist.tsv",
    "metadata": f"{output_dir}metadata.tsv",
    "rna_pca": f"{output_dir}rna_pca.tsv",
    "atac_features": f"{output_dir}atac_features.tsv",
    "atac_matrix": f"{output_dir}atac_matrix.mm",
}

outputs = {
    # "h5_file": f"{output_dir}immature.h5"
    "h5_file": f"{os.getenv('PROJECT_PATH')}/notebooks/Immature/SAILERX/data/immature.h5"
}

# %%
print(f"{inputs=}")
print(f"{outputs=}")

# %%
# ls {output_dir}

# %% [markdown]
# ## Check their h5 file, to see what target we have

# %%
filename = inputs["their_h5"]

# %%
# !h5dump -A {filename}

# %%
their_h5 = h5py.File(filename, "r")

# %%
head = {
    y: {x: their_h5[f"{y}/{x}"][:10] for x in list(their_h5[y].keys())}
    for y in list(their_h5.keys())
    if y not in ["pca"]
}

# %%
head

# %% [markdown]
# ## Import the 10X h5 file
#
# Using 10X's [method](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices), slightly altered (CSR instead of CSC).

# %%
CountMatrix = collections.namedtuple(
    "CountMatrix", ["feature_ref", "barcodes", "matrix"]
)


def get_matrix_from_h5(filename):
    with tables.open_file(filename, "r") as f:
        mat_group = f.get_node(f.root, "matrix")
        barcodes = f.get_node(mat_group, "barcodes").read()
        data = getattr(mat_group, "data").read()
        indices = getattr(mat_group, "indices").read()
        indptr = getattr(mat_group, "indptr").read()
        shape = getattr(mat_group, "shape").read()[::-1]
        matrix = sp.csr_matrix((data, indices, indptr), shape=shape)

        feature_ref = {}
        feature_group = f.get_node(mat_group, "features")
        feature_ids = getattr(feature_group, "id").read()
        feature_names = getattr(feature_group, "name").read()
        feature_types = getattr(feature_group, "feature_type").read()
        feature_ref["id"] = feature_ids
        feature_ref["name"] = feature_names
        feature_ref["feature_type"] = feature_types
        tag_keys = getattr(feature_group, "_all_tag_keys").read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()

        return CountMatrix(feature_ref, barcodes, matrix)


# %%
filtered_matrix_h5 = inputs["cellranger_h5"]
mat = get_matrix_from_h5(filtered_matrix_h5)

# %%
[x for x in dir(mat) if not x.startswith("_")]

# %% [markdown]
# ## Import barcodes for \<Sample\> cells:

# %%
# !wc -l {inputs["cell_whitelist"]}

# %%
# !head {inputs["cell_whitelist"]}

# %%
with open(inputs["cell_whitelist"], "r") as f:
    immature_cells = [x.strip() for x in f.readlines()]

# %%
barcodes_h5 = [x.decode("UTF-8") for x in mat.barcodes.tolist()]
barcodes_h5[-5:], len(barcodes_h5)

# %%
immature_cells[:5], len(immature_cells)

# %%
sum(el in barcodes_h5 for el in immature_cells)

# %%
len(set(barcodes_h5).difference(set(immature_cells)))

# %%
sorted(list(set(barcodes_h5).intersection(set(immature_cells)))) == immature_cells

# %% [markdown]
# ## Do subset of cells to only current sample

# %%
save_cells = np.where(np.in1d(barcodes_h5, immature_cells))[0]

# %%
mat.matrix[
    save_cells,
]

# %%
mat.barcodes[save_cells], len(mat.barcodes[save_cells])

# %% [markdown]
# Neat.

# %% [markdown]
# ## Import metadata for \<Sample\> cells

# %%
# !head -n 2 {inputs["metadata"]}

# %%
immature_metadata = pd.read_table(inputs["metadata"], index_col="barcode")

# %%
immature_metadata.head()

# %% [markdown]
# ## Import PCA embeddings for RNA data

# %%
pca = pd.read_table(inputs["rna_pca"], index_col="barcode")
pca

# %%
assert all(pca.index == immature_cells)

# %% [markdown]
# Same order of cells, so it should be fine to just create this as a new hdf5 Group.

# %% [markdown]
# ## Read ATAC matrix

# %%
atac = mmread(inputs["atac_matrix"])

# %%
atac

# %%
with open(inputs["atac_features"], "r") as f:
    feature_intervals = [x.strip() for x in f.readlines()]
feature_intervals[-5:]

# %% tags=[]
feature_intervals = [x.replace("-", ":", 1) for x in feature_intervals]
feature_intervals[-5:]

# %%
mat.matrix[
    save_cells,
]

# %%
mat.feature_ref["interval"]

# %%
save_indices = [True for i in range(len(mat.feature_ref["interval"]))]

# %% [markdown]
# We need to subset the featurelist into only atac features/peaks:

# %%
Counter(mat.feature_ref["feature_type"])

# %% tags=[]
bool_peaks = np.array([x.decode("UTF-8") for x in mat.feature_ref["feature_type"]]) == "Peaks"
peak_indices = np.where(
    bool_peaks
)[0]
bool_peaks, peak_indices, len(peak_indices)

# %%
Counter(mat.feature_ref["feature_type"][peak_indices])

# %%
features_h5 = [x.decode("UTF-8") for x in mat.feature_ref["interval"].tolist()]
features_h5[:5], features_h5[-5:], len(features_h5)

# %% [markdown]
# We also need to remove the non-standard chromosomes (e.g. GL...):

# %%
bool_std_chrs =  np.char.startswith(np.array(features_h5), "chr")
standard_chromosome_intervals = np.where(
    bool_std_chrs
)[0]

# %%
standard_chromosome_intervals = np.char.startswith(np.array(features_h5)[np.where(peak_indices)[0]], "chr")

# %% [markdown]
# Make sure the non-standard chromosomes are subset correctly:

# %%
Counter(np.char.startswith(np.array(features_h5)[np.where(~bool_peaks)[0]], "chr"))

# %%
Counter(np.char.startswith(np.array(features_h5)[np.where(bool_peaks)[0]], "chr"))

# %%
assert (
    Counter(
        np.char.startswith(
            np.array(features_h5)[np.where(~bool_peaks)[0]], "chr"
        )
    )[False]
    + Counter(
        np.char.startswith(
            np.array(features_h5)[np.where(bool_peaks)[0]], "chr"
        )
    )[False]
    == Counter(bool_std_chrs)[False]
)

# %% [markdown]
# If no errors, great.
#
# ---

# %% [markdown]
# Add both sets of indices together, to know which indices to subset the feature_* arrays with in the feature_ref dict:

# %%
np.equal(bool_peaks, bool_std_chrs)

# %%
Counter(bool_peaks & bool_std_chrs)

# %%
save_indices = np.where(bool_peaks & bool_std_chrs)[0]

# %%
save_indices

# %% [markdown]
# Finally......
#
# ---

# %% [markdown]
# ## Let's subset!

# %% [markdown]
# Get the target shape to aim for:

# %%
target_shape = len(save_cells), len(save_indices)
target_shape

# %%
{x: len(mat.feature_ref[x]) for x in list(mat.feature_ref.keys())}

# %% [markdown]
# So we should subset all structures in `mat.feature_ref`.

# %%
subset_matrix = mat.matrix[save_cells, ][:, save_indices]
subset_matrix

# %%
assert subset_matrix.shape == target_shape

# %%
subset_features = {k: v[save_indices] for k, v in mat.feature_ref.items()}
subset_features

# %%
for v in subset_features.values():
    assert v.shape[0] == target_shape[1] 

# %%
assert pca.shape[0] == target_shape[0]

# %% [markdown]
# If no errors, fantastic.

# %% [markdown]
# Then we are only missing `chrom_size`s?

# %% [markdown]
# ## Construct `chrom_sizes`
#
# Which is not, as I first thought, the total sizes of the chromosomes in the genome. It is the _number of peaks_ that exist in each chromosome _in the peak/feature set in this specific data set_.
#
# E.g. if feature set is 10 features long, the distribution might be:
#
# ```python
# {
#     "chr1": 5,
#     "chr2": 3,
#     "chr15": 2
# }
# ```

# %%
primary_chroms = ["chr" + x for x in [str(i) for i in range(1, 20)] + ["X", "Y"]]


# %%
def sort_chrom_wise(df):
    return df.sort_values(
        by="chrom",
        key=lambda column: column.map(
            {v: i for i, v in enumerate(primary_chroms)}
        ),
    ).reset_index(drop=True)


# %%
subset_features["interval"][:10]

# %%
chrom_feature_distribution = [
    re.sub("chr", "", z)
    for z in [
        re.sub(":.*", "", y)
        for y in [
            x.decode("UTF-8") for x in subset_features["interval"].tolist()
        ]
    ]
]

# %%
chrom_feature_distribution = {"chr" + str(k): v for k, v in dict(Counter(chrom_feature_distribution)).items()}

# %%
chrom_feature_distribution = {k: chrom_feature_distribution[k] for k in primary_chroms}

# %%
chrom_feature_distribution

# %%
assert sum(chrom_feature_distribution.values()) == target_shape[1]

# %% [markdown]
# ---

# %% [markdown]
# ## Create own, new h5 file (pyhdf5 way)

# %%
# !rm {outputs["h5_file"]}

# %%
own_h5 = h5py.File(outputs["h5_file"], "w")

# %%
own_h5

# %%
own_h5.create_dataset_like(
    "cells/labels", their_h5["cells/labels"], shape=len(immature_cells)
)

# %%
name = "cells/rna_depth"
own_h5.create_dataset_like(name, their_h5[name], shape=len(immature_cells))

# %%
name = "pca/data"
own_h5.create_dataset_like(name, their_h5[name], shape=pca.shape)

# %%
name = "atac/chrom_size"
own_h5.create_dataset_like(name, their_h5[name], shape=len(chrom_feature_distribution))

# %%
name = "atac/data"
own_h5.create_dataset_like(name, their_h5[name], shape=subset_matrix.data.shape)

# %%
name = "atac/indptr"
own_h5.create_dataset_like(name, their_h5[name], shape=subset_matrix.indptr.shape)

# %%
name = "atac/indices"
own_h5.create_dataset_like(name, their_h5[name], shape=subset_matrix.indices.shape)

# %%
own_h5["atac"].attrs.create(name="shape", data=target_shape)
# wait with this, create when either atac/data or atac/indices have been created and use their shapes

# %%
own_h5["cells/labels"][:] = list(immature_metadata["seurat_clusters"])

# %%
own_h5["cells/labels"].shape, their_h5["cells/labels"].shape

# %% [markdown]
# Seems to have worked.

# %%
own_h5["cells/rna_depth"][:] = list(immature_metadata["nCount_RNA"])

# %%
own_h5["pca/data"][:] = np.array(pca)

# %%
own_h5["atac/data"][:] = subset_matrix.data

# %%
own_h5["atac/indptr"][:] = subset_matrix.indptr

# %%
own_h5["atac/indices"][:] = subset_matrix.indices

# %%
own_h5["atac/chrom_size"][:] = list(chrom_feature_distribution.values())

# %%
own_h5["atac/chrom_size"][:]

# %%
own_h5.close()

# %%
# own_h5 = h5py.File(outputs["h5_file"], "a")

# %%
# !h5dump -A {outputs["h5_file"]}

# %% [markdown]
# Looks great! Nice!

# %% [markdown]
# # Done!

# %%
"Notebook finished at " + datetime.datetime.now(pytz.timezone("Europe/Stockholm")).strftime("%Y-%m-%d %H:%M")

# %% [markdown]
# ---
