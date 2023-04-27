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

# %%
from datetime import timedelta as td
from datetime import datetime as dt

# %%
start_time = dt.now()
print("Notebook started at " + start_time.strftime("%Y-%m-%d %H:%M:%S"))

# %% [markdown]
# # Evaluating SAILERX results

# %%
import re
import os

import numpy as np
import pandas as pd

# %%
output_dir = os.getenv("PROJECT_PATH") + "/data/adhoc/lower_features/"
sailerx_dir = f"./SAILERX/DB_QZ/DB_QZ/"

inputs = {
    "embedding": f"{sailerx_dir}embedding.npy",
    "embedding_default": f"{sailerx_dir}default_embedding.npy",
    "cells": f"{output_dir}cell_whitelist.tsv"
}

outputs = {
    "embedding_csv": f"{output_dir}sailerx_embeddings.csv"
}

# %%
# ls {inputs["cells"]}

# %%
# !wc -l {inputs["cells"]}

# %%
embedding = np.load(inputs["embedding"], allow_pickle=False)
embedding

# %%
embedding.shape

# %% [markdown] tags=[]
# ## Export embeddings to tsv

# %%
df = pd.DataFrame(embedding)
df

# %%
# !ls {os.path.dirname(inputs["cells"])}

# %%
cells = pd.read_table(inputs["cells"], header=None)

# %%
df["cells"] = cells[0]

# %%
df.set_index("cells", inplace=True)

# %%
df

# %%
df.to_csv(outputs["embedding_csv"])

# %%
# ls {os.path.dirname(outputs["embedding_csv"])}

# %%
# !head {outputs["embedding_csv"]} -n 3

# %% [markdown]
# ---

# %%
end_time = dt.now()
print("Notebook finished at " + end_time.strftime("%Y-%m-%d %H:%M:%S"))
time_diff = (end_time - start_time)
print(f"Notebook run time was {round(time_diff.seconds / 60)} minutes, {time_diff.seconds % 60} seconds.")
