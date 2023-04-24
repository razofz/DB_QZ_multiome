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
# # Evaluating SAILERX results

# %%
import datetime
import pytz

"Notebook last run at " + datetime.datetime.now(pytz.timezone("Europe/Stockholm")).strftime("%Y-%m-%d %H:%M")

# %%
import re
import os

import numpy as np
import pandas as pd

# %%
output_dir = os.getenv("PROJECT_PATH") + "/data/processed/notebooks/Immature/"
sailerx_dir = f"./SAILERX/DB_QZ/DB_QZ/"

inputs = {
    "embedding": f"{sailerx_dir}embedding.npy",
    "embedding_default": f"{sailerx_dir}default_embedding.npy",
    "immature_cells": f"{output_dir}cell_whitelist.tsv"
}

outputs = {
    "embedding_csv": f"{output_dir}sailerx_embeddings.csv"
}

# %%
# ls {inputs["immature_cells"]}

# %%
# !wc -l {inputs["immature_cells"]}

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
# !ls {os.path.dirname(inputs["immature_cells"])}

# %%
cells = pd.read_table(inputs["immature_cells"], header=None)

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

# %%
"Notebook finished at " + datetime.datetime.now(pytz.timezone("Europe/Stockholm")).strftime("%Y-%m-%d %H:%M")

# %% [markdown]
# ---
