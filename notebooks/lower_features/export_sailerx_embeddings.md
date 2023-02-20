---
jupyter:
  jupytext:
    formats: ipynb,py:percent,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Evaluating SAILERX results

```python
import re
import os
import statistics

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import umap
from IPython.display import SVG
from sklearn.neighbors import kneighbors_graph
# from sknetwork.clustering import Louvain
# from sknetwork.utils import CNNDense, KNNDense
# from sknetwork.visualization import svg_graph
```

```python
output_dir = os.getenv("PROJECT_PATH") + "/data/adhoc/lower_features/"
sailerx_dir = f"./SAILERX/DB_QZ/DB_QZ/"

inputs = {
    "embedding": f"{sailerx_dir}embedding.npy",
    "embedding_default": f"{sailerx_dir}default_embedding.npy",
    "diverse_cells": f"{output_dir}cell_whitelist.tsv"
}

outputs = {
    "embedding_csv": f"{output_dir}sailerx_embeddings.csv"
}
```

```python
embedding = np.load(inputs["embedding"], allow_pickle=False)
embedding
```

```python
embedding_default = np.load(inputs["embedding_default"], allow_pickle=False)
embedding_default
```

```python
embedding.shape
```

```python
embedding_default.shape
```

<!-- #region tags=[] -->
## Export embeddings to tsv
<!-- #endregion -->

```python
df = pd.DataFrame(embedding)
df
```

```python
!ls {os.path.dirname(inputs["diverse_cells"])}
```

```python
cells = pd.read_table(inputs["diverse_cells"], header=None)
```

```python
df["cells"] = cells[0]
```

```python
df.set_index("cells", inplace=True)
```

```python
df
```

```python
df.to_csv(outputs["embedding_csv"])
```

```python
ls {os.path.dirname(outputs["embedding_csv"])}
```

```python
!head {outputs["embedding_csv"]} -n 3
```
