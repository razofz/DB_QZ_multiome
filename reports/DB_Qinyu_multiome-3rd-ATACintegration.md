---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python slideshow={"slide_type": "skip"} tags=[]
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
from IPython.display import Image
from IPython.core.display import SVG
from IPython.display import display, Markdown, HTML
import IPython
```

```python slideshow={"slide_type": "skip"} tags=[]
import pandas as pd
import os
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
# Multiome data (RNA+ATAC), Qinyu Zhang @ David Bryder lab

Meeting 2022-02-10
<!-- #endregion -->

```python slideshow={"slide_type": "skip"} tags=[]
ls data/interim/adhoc/ATAC-integration/
```

```python slideshow={"slide_type": "skip"} tags=[]
ls data/interim/adhoc/ATAC-integration/python
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
## Focus on ATAC integration and cluster annotation

Agreed-upon focus on integrating the two samples' ATAC modalities, as well as annotating RNA clusters based on differentially expressed genes, with the help of BM genes and given list of 10 genes.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
What was accomplished was integration of ATAC modalities, as well as formalising what has been done this far, which means higher reproducibility and the data being ready for annotation efforts.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Naming

The two samples:

1. Immature (formerly EPCR)
2. Diverse (Formerly Total, formerly Viable_RNA)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Integration methods for ATAC tried

- Harmony
- Partial PCA
- Projection
    - k-Nearest Neighbours (Scarf)
    - Reciprocal LSI (Signac/Seurat), aka anchor-based approach
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Not yet integrated (baseline/starting point)
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''![UMAP](data/interim/adhoc/ATAC-integration/umap-aggregated.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "skip"} tags=[]
IPython.display.HTML("<img src='data/interim/adhoc/ATAC-integration/umap-aggregated.svg', height='80px'>")
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### rLSI/anchor-based approach
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''![UMAP](data/interim/adhoc/ATAC-integration/umap-rlsi-anchors-only.svg)'''
display(Markdown(md))
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Harmony
<!-- #endregion -->

```python slideshow={"slide_type": "slide"} tags=[]
md = '''![UMAP](data/interim/adhoc/ATAC-integration/umap-harmony.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Unintegrated
![UMAP](data/interim/adhoc/ATAC-integration/python/umap-aggregated.svg)
'''
display(Markdown(md))
```

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Partial PCA
![UMAP](data/interim/adhoc/ATAC-integration/python/umap-pPCA.svg)
'''
display(Markdown(md))
```

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### kNN-based Projection
![UMAP](data/interim/adhoc/ATAC-integration/python/umap-projection-scarf.svg)

This looks not good at all, however not entirely unprobable there is something I overlooked. Nevertheless, the Harmony results looks good enough that I think we should go forward with that.
'''
display(Markdown(md))
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Discussion

Similarly to for the RNA data, I think Harmony gave the best result, and is what we should move forward with.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
![UMAP](data/interim/adhoc/ATAC-integration/umap-harmony.svg)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
by Rasmus Olofzon
---
<!-- #endregion -->

```python slideshow={"slide_type": "skip"} tags=[]
#!jupyter nbconvert --to slides DB_Qinyu_multiome-2nd-integration.ipynb --config slides-config.py
```

```python slideshow={"slide_type": "skip"} tags=[]
#!jupyter nbconvert --to slides DB_Qinyu_multiome-2nd-integration.ipynb --config slides-config.py --post serve
# and add `?print-pdf` to url, then print to pdf file
```
