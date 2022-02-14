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
from IPython.display import display, Markdown
```

```python slideshow={"slide_type": "skip"} tags=[]
import pandas as pd
import os
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
# Multiome data (RNA+ATAC), Qinyu Zhang @ David Bryder lab

Meeting 2022-02-08
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
hello
<!-- #endregion -->

```python slideshow={"slide_type": "skip"} tags=[]
ls data/interim/adhoc
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
## Focus on integration

Focus on integrating/merging the two samples, agreed upon in our meeting last week.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
The two samples:

1. EPCR
2. Total (aka Viable_RNA)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
There are two conceptual approaches to integrating these multi-modal data.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
- Modality-based
    - Integrate the samples modality-wise
    - (EPCR RNA with Total RNA, and EPCR ATAC with Total ATAC)
    - then integrate the result's modalities (RNA & ATAC) with each other.
    - (joint-RNA with joint-ATAC)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
- Sample-based
    - Integrate the modalities sample-wise,
    - (EPCR RNA with EPCR ATAC, and Total RNA with Total ATAC)
    - and then integrate the two samples with each other.
    - (EPCR-joint with Total-joint)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
For now I went with the first approach, modality-based. Specifically, integrating EPCR RNA with Total RNA.
<!-- #endregion -->

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Not yet integrated (baseline/starting point)
![UMAP](data/interim/adhoc/umap_origin.svg)'''
display(Markdown(md))
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
I tried a couple of approaches:
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
- Integration through projection (EPCR on total)
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
- Integration through Partial PCA
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
- Integration with Harmony
<!-- #endregion -->

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Partial PCA
![UMAP](data/interim/adhoc/umap-partial-pca.svg)

Not very good. Barely merged.'''
display(Markdown(md))
```

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''#### Partial PCA, clustering
![UMAP](data/interim/adhoc/umap-partial-pca-clusters.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''#### Cell cycle stage
![UMAP](data/interim/adhoc/pumap-cc.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Projection (scarf), "unified UMAP"
![UMAP](data/interim/adhoc/unified-umap.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''#### Unified UMAP with Total clusters
![UMAP](data/interim/adhoc/unified-umap-total-w-total-clusters.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "subslide"} tags=[]
md = '''#### Unified UMAP with Total clusters, only EPCR cells shown
![UMAP](data/interim/adhoc/unified-umap-EPCR-transferred-labels-from-total.svg)'''
display(Markdown(md))
```

```python slideshow={"slide_type": "slide"} tags=[]
md = '''### Harmony
![UMAP](data/interim/adhoc/R/umap-harmony-integrated.svg)'''
display(Markdown(md))
```

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
### Discussion

Projection gave a more evenly distribution of EPCR cells over Total cells, while Harmony placed them mostly to one end of the data. Biologically, this seems to make more sense than from the projection approach, considering that EPCR cells should be mostly immature cells while Total should be more diverse. Based on this, I would say we go with Harmony.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "fragment"} tags=[] -->
Of course, there are a lot of integration methods (and also conceptual approaches as I mentioned before). We can try others, but from our discussions so far I think the Harmony result align with the biological interpretation, which should be good to base further analysis on.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} tags=[] -->
---
<!-- #endregion -->

```python slideshow={"slide_type": "skip"} tags=[]
#!jupyter nbconvert --to slides DB_Qinyu_multiome-2nd-integration.ipynb --config slides-config.py
```

```python
#!jupyter nbconvert --to slides DB_Qinyu_multiome-2nd-integration.ipynb --config slides-config.py --post serve
# and add `?print-pdf` to url, then print to pdf file
```
