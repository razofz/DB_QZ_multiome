# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + slideshow={"slide_type": "skip"} tags=[]
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
from IPython.display import Image
from IPython.core.display import SVG
from IPython.display import display, Markdown

# + slideshow={"slide_type": "skip"} tags=[]
import pandas as pd
import os

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# # Multiome data (RNA+ATAC), Qinyu @ David Bryder lab

# + slideshow={"slide_type": "skip"} tags=[]
# ls data/processed/GFP_Viable_RNA/

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Experiment layout

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# Two samples:
#
# 1. _EPCR_
# 2. *Viable_RNA*

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Overview of analysis

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# _(After Cellranger pipeline, i. e. aligning, peak calling etc have already been done. In other words, this analysis starts from_ count matrices_.)_

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### For both modalities (RNA & ATAC)

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Initial filtering**
#     - (e. g. for a feature to be included it has to have been expressed/present in at least 10 cells)

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Filtering** based on metadata
#     - (e. g. number of RNA counts per cell, % mitochondrial RNA per cell, number of ATAC features per cell)

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Normalization**
#     - library-size normalization for RNA and tf-idf for ATAC 

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - finding **Highly Variable Genes** for RNA, and **Prevalent Peaks** for ATAC
#     - (using the normalized values)

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Linear dimensionality reduction**
#     - i. e. PCA for RNA and LSI for ATAC

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Non-linear dimensionality reduction**
#     - i. e. UMAP

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Clustering**

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### RNA-specific

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - **Differentially expressed genes**, aka **DE** genes
#     - for each cluster, identify the genes with a markedly higher expression compared to all the other clusters

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### ATAC-specific

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# - Adding **Gene Scores**
#     - use GRCm38 (mouse reference) annotations for gene bodies and their corresponding promoter region, and identify which peaks fall into these regions. This gives a Gene Score, which can be used to better understand the dataset, since more is known about marker genes than non-coding regions.

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Filtering results 
# <small><i>(metadata)</i></small>

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### Viable_RNA

# + slideshow={"slide_type": "fragment"} tags=[]
Image('data/processed/GFP_Viable_RNA/cell_dists_prefiltering.png')
Image('data/processed/GFP_Viable_RNA/cell_dists_postfiltering.png')

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### EPCR_LSKs_RNA

# + slideshow={"slide_type": "fragment"} tags=[]
Image('data/processed/GFP_EPCR_LSKs_RNA/cell_dists_prefiltering.png')
Image('data/processed/GFP_EPCR_LSKs_RNA/cell_dists_postfiltering.png')

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Highly variable genes (hvgs) selection

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### EPCR
![UMAP](data/processed/GFP_EPCR_LSKs_RNA/hvgs_volcano.png)'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### Viable_RNA
![UMAP](data/processed/GFP_Viable_RNA/hvgs_volcano.png)'''
display(Markdown(md))

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## UMAPs for Viable_RNA
# *Non-linear dimensionality reduction*

# + slideshow={"slide_type": "subslide"} tags=[]
SVG('data/processed/GFP_Viable_RNA/umap_RNA.svg')

# + slideshow={"slide_type": "subslide"} tags=[]
SVG('data/processed/GFP_Viable_RNA/umap_ATAC.svg')

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## UMAPs for EPCR
# *Non-linear dimensionality reduction*

# + slideshow={"slide_type": "subslide"} tags=[]
SVG('data/processed/GFP_EPCR_LSKs_RNA/umap_RNA.svg')

# + slideshow={"slide_type": "subslide"} tags=[]
SVG('data/processed/GFP_EPCR_LSKs_RNA/umap_ATAC.svg')

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## DE genes (RNA)

# + slideshow={"slide_type": "skip"} tags=[]
from IPython.display import display, Markdown
#display(Markdown('**_some_ markdown** and an !'))

# + slideshow={"slide_type": "skip"} tags=[]
md = '''<table>
    <tr>
        <th> [GFP_Viable_RNA] </th>
        <th> [GFP_EPCR] </th>
    </tr>
    <tr style="height:40%">
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_Viable_RNA/DE_genes_heatmap.png" /> </div> </td>
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_EPCR_LSKs_RNA/DE_genes_heatmap.png"  /> </div> </td>
    </tr>
</table>'''

# + slideshow={"slide_type": "subslide"} tags=[]
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=[]
print("Viable_RNA")
pd.read_csv('data/processed/GFP_Viable_RNA/top_deg_genes.csv', index_col=0)

# + slideshow={"slide_type": "subslide"} tags=[]
print("EPCR")
pd.read_csv('data/processed/GFP_EPCR_LSKs_RNA/top_deg_genes.csv', index_col=0)

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## DE genes visualized on UMAPs

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### Viable_RNA

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_Viable_RNA/deg_plots/cluster-1_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 1, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_Viable_RNA/deg_plots/cluster-2_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 2, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_Viable_RNA/deg_plots/cluster-3_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 3, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_Viable_RNA/deg_plots/cluster-4_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 4, 1st de gene

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### EPCR

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_EPCR_LSKs_RNA/deg_plots/cluster-1_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 1, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_EPCR_LSKs_RNA/deg_plots/cluster-2_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 2, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_EPCR_LSKs_RNA/deg_plots/cluster-3_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 3, 1st de gene

# + slideshow={"slide_type": "subslide"} tags=[]
fn = "data/processed/GFP_EPCR_LSKs_RNA/deg_plots/cluster-4_de-0.svg"
if os.path.exists(fn): SVG(fn)
else: "DE gene for cluster not found"
# Cluster 4, 1st de gene

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Gene scores (ATAC)

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### Viable_RNA
![UMAP](data/processed/GFP_Viable_RNA/umaps_ATAC_GeneScores.svg)'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### EPCR
![UMAP](data/processed/GFP_EPCR_LSKs_RNA/umaps_ATAC_GeneScores.svg)'''
display(Markdown(md))

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Comparison of RNA and ATAC modalities

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# #### Clusters

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''<table class="r-stretch" style="height:80%">
    <tr>
        <th> ATAC clusters <i>[Viable_RNA]</i> </th>
        <th>  </th>
    </tr>
    <tr style="height:40%">
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_Viable_RNA/umap_ATAC.svg" /> </div> </td>
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_Viable_RNA/umap_RNA-w-ATAC-clusters.svg" /> </div> </td>
    </tr>
</table>
'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=[]
md = '''<table class="r-stretch" style="height:80%">
    <tr>
        <th> RNA clusters <i>[Viable_RNA]</i> </th>
        <th>  </th>
    </tr>
    <tr style="height:40%">
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_Viable_RNA/umap_ATAC-w-RNA-clusters.svg" /> </div> </td>
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_Viable_RNA/umap_RNA.svg"  /> </div> </td>
    </tr>
</table>
'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''<table>
    <tr>
        <th> ATAC clusters <i>[EPCR]</i> </th>
        <th> RNA clusters </th>
    </tr>
    <tr style="height:40%">
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_EPCR_LSKs_RNA/umap_ATAC.svg" /> </div> </td>
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_EPCR_LSKs_RNA/umap_RNA-w-ATAC-clusters.svg" /> </div> </td>
    </tr>
</table>
'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=[]
md = '''<table>
    <tr>
        <th> RNA clusters <i>[EPCR]</i> </th>
        <th>  </th>
    </tr>
    <tr style="height:40%">
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_EPCR_LSKs_RNA/umap_ATAC-w-RNA-clusters.svg" /> </div> </td>
        <td style="height:60%"> <div style="height:80%"> <img src="data/processed/GFP_EPCR_LSKs_RNA/umap_RNA.svg"  /> </div> </td>
    </tr>
</table>
'''
display(Markdown(md))

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### Viable_RNA

# + slideshow={"slide_type": "fragment"} tags=[]
pd.read_csv('data/processed/GFP_Viable_RNA/cells-in-clusters_RNA.csv', index_col=0, names=['Number of cells'])
pd.read_csv('data/processed/GFP_Viable_RNA/cells-in-clusters_ATAC.csv', index_col=0, names=['Number of cells'])

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# #### Overlap between RNA and ATAC clusters
# Columns: ATAC clusters   
# Rows: RNA clusters

# + slideshow={"slide_type": "fragment"} tags=[]
pd.read_csv('data/processed/GFP_Viable_RNA/clusters_RNA-ATAC_crosstab.csv', index_col=0)

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# ### EPCR

# + slideshow={"slide_type": "fragment"} tags=[]
pd.read_csv('data/processed/GFP_EPCR_LSKs_RNA/cells-in-clusters_RNA.csv', index_col=0, names=['Number of cells'])
pd.read_csv('data/processed/GFP_EPCR_LSKs_RNA/cells-in-clusters_ATAC.csv', index_col=0, names=['Number of cells'])

# + [markdown] slideshow={"slide_type": "subslide"} tags=[]
# #### Overlap between RNA and ATAC clusters
# Columns: ATAC clusters   
# Rows: RNA clusters

# + slideshow={"slide_type": "fragment"} tags=[]
pd.read_csv('data/processed/GFP_EPCR_LSKs_RNA/clusters_RNA-ATAC_crosstab.csv', index_col=0)

# + [markdown] slideshow={"slide_type": "slide"} tags=[]
# ## Comparison of RNA and ATAC modalities

# + [markdown] slideshow={"slide_type": "fragment"} tags=[]
# #### Gene expression

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### Viable_RNA
![UMAP](data/processed/GFP_Viable_RNA/umaps_RNA-ATAC_geneexp.svg)'''
display(Markdown(md))

# + slideshow={"slide_type": "subslide"} tags=["to-remove"]
md = '''### Viable_RNA
![UMAP](data/processed/GFP_EPCR_LSKs_RNA/umaps_RNA-ATAC_geneexp.svg)'''
display(Markdown(md))

# + [markdown] slideshow={"slide_type": "skip"} tags=[]
# ---

# + [markdown] slideshow={"slide_type": "skip"} tags=[]
# <div style="height:75%">
#     <img src="data/processed/GFP_Viable_RNA/umaps_RNA-ATAC_clusters.svg"  />
# </div>
#
# <small><i> [GFP_Viable_RNA] </i></small>

# + [markdown] slideshow={"slide_type": "skip"} tags=[]
# <div style="height:75%">
#     <img src="data/processed/GFP_EPCR_LSKs_RNA/umaps_RNA-ATAC_clusters.svg"  />
# </div>
#
# <small><i> [GFP_EPCR_LSKs_RNA] </i></small>

# + [markdown] slideshow={"slide_type": "skip"} tags=[]
# <div style="height:40%">
#     <img src="data/processed/GFP_Viable_RNA/DE_genes_heatmap.png"  />
# </div>
#
# <small><i> [GFP_Viable_RNA] </i></small>

# + [markdown] slideshow={"slide_type": "skip"} tags=[]
# | ![UMAP](data/processed/GFP_Viable_RNA/umaps_RNA-ATAC_clusters.svg) | . |
# |---|---|
# | UMAPs for RNA and ATAC datasets |  |

# + slideshow={"slide_type": "skip"} tags=[]
# !jupyter nbconvert --to slides DB_Qinyu_multiome.ipynb --config slides-config.py
