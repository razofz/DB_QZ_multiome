import scarf
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math
import numpy as np


args = str(sys.argv)
valid_sample_names = ["GFP_EPCR_LSKs_RNA", "GFP_Viable_RNA"]
if len(sys.argv) == 1:
    sys.exit(
        f'No sample name specified. Pick one of {valid_sample_names}.'
    )
elif sys.argv[1] not in valid_sample_names:
    sys.exit(
        f'"{sys.argv[1]}" is not a valid sample name. Pick one of {valid_sample_names}.'
    )
else:
    sample = sys.argv[1]

try:
    os.makedirs(f"data/processed/{sample}/deg_plots/", exist_ok=False)
except OSError:
    os.system(f"rm -rf data/processed/{sample}/deg_plots/")
    os.makedirs(f"data/processed/{sample}/deg_plots/", exist_ok=False)
# os.makedirs(f"data/processed/{sample}/leiden/hvg_plots/", exist_ok=True)
# os.makedirs(f"data/processed/{sample}/paris/hvg_plots/", exist_ok=True)

if sample == "GFP_EPCR_LSKs_RNA":
    cmap_RNA = 'nipy_spectral'
    cmap_ATAC = 'terrain'
else:
    cmap_RNA = 'tab20'
    cmap_ATAC = 'turbo'

if len(sys.argv) >= 3 and sys.argv[2] == "cache":
    pass
else:
    reader = scarf.CrH5Reader(
        f"data/raw/count/{sample}/outs/filtered_feature_bc_matrix.h5"
    )
    print(reader.assayFeats)
    writer = scarf.CrToZarr(
        reader,
        zarr_fn=f"data/processed/{sample}/counts_{sample}.zarr",
        chunk_size=(6000, 2000),
    )
    writer.dump(batch_size=12000)

ds = scarf.DataStore(
    zarr_loc=f"data/processed/{sample}/counts_{sample}.zarr",
    default_assay="RNA",
    nthreads=16,
)
ds
ds.RNA.feats.head()
ds.ATAC.feats.head()
# ds.plot_cells_dists(cols=['ATAC_nCounts', 'ATAC_nFeatures', 'RNA_nCounts',
#                           'RNA_nFeatures', 'RNA_percentMito',
#                           'RNA_percentRibo'])
# ds.plot_cells_dists(cols=['ATAC_nCounts', 'ATAC_nFeatures', 'RNA_nCounts',
#                           'RNA_nFeatures', 'RNA_percentMito',
#                           'RNA_percentRibo'],
#                     cell_key='I')
ds.auto_filter_cells(attrs=['ATAC_nCounts', 'ATAC_nFeatures', 'RNA_nCounts',
                            'RNA_nFeatures', 'RNA_percentMito',
                            'RNA_percentRibo'], show_qc_plots=False)
# if sample == "GFP_Viable_RNA":
#     ds.filter_cells(attrs=['ATAC_nCounts', 'ATAC_nFeatures', 'RNA_nCounts',
#                            'RNA_nFeatures', 'RNA_percentMito',
#                            'RNA_percentRibo'],
#                     lows=[0, 0, 0, 0, 0, 0],
#                     highs=[125000, 35000, 30000, 7000, 30, 35])
# ds.plot_cells_dists(cols=['ATAC_nCounts', 'ATAC_nFeatures', 'RNA_nCounts',
#                           'RNA_nFeatures', 'RNA_percentMito',
#                           'RNA_percentRibo'],
#                     cell_key='I')

# set all sample-specific parameters
if sample == "GFP_Viable_RNA":
    # hvg_search
    max_var = 5
    max_mean = 3
    # graph making
    k = 11
    dims = 11
    n_centroids = 500
    # clustering
    resolution_RNA = 0.6
    resolution_ATAC = 0.7
elif sample == "GFP_EPCR_LSKs_RNA":
    # hvg_search
    max_var = 5
    max_mean = 9999
    # graph making
    k = 15
    dims = 5
    n_centroids = 800
    # clustering
    resolution_RNA = 0.7
    resolution_ATAC = 0.7

ds.mark_hvgs(
    min_cells=20,
    top_n=500,
    max_var=max_var,
    max_mean=max_mean,
    show_plot=False)

ds.make_graph(
    feat_key="hvgs",
    k=k,
    dims=dims,
    n_centroids=n_centroids
)

ds.run_umap(
    # n_epochs=200,
    # spread=5,
    # min_dist=1,
    parallel=True,
)

print("running leiden clustering..")

ds.run_leiden_clustering(resolution=resolution_RNA)

leiden_clusters = ds.cells.to_pandas_dataframe(
    columns=["RNA_leiden_cluster"], key="I")
n_clusters = len(set(leiden_clusters['RNA_leiden_cluster']))

# ds.run_marker_search(group_key="RNA_leiden_cluster", threshold=0.25)
ds.run_marker_search(group_key="RNA_leiden_cluster", threshold=0.20)

length_of_DE_genes = {
    str(i): len(ds.get_markers(
        group_key='RNA_leiden_cluster',
        group_id=str(i)))
    for i in range(1, n_clusters+1)}
print(f"{length_of_DE_genes=}")

ds.plot_marker_heatmap(
    group_key="RNA_leiden_cluster",
    topn=5,
    figsize=(5, 9),
    show_fig=False,
    savename=f"data/processed/{sample}/DE_genes_heatmap.png",
)

ds.plot_layout(
    layout_key="RNA_UMAP",
    color_by="RNA_leiden_cluster",
    show_fig=False,
    cmap=cmap_RNA,
    savename=f"data/processed/{sample}/umap_RNA.svg",
)

# # This should fail
# assert False == True

cells_in_clusters = leiden_clusters.value_counts()
clusters = len(leiden_clusters.value_counts())
print(f"{leiden_clusters.nunique()[0]=}")
print(f"{clusters=}")
print(f"{cells_in_clusters=}")

for cl in range(1, leiden_clusters.nunique()[0] + 1):
    markers = ds.get_markers(group_key="RNA_leiden_cluster", group_id=cl)["names"]

    for de in range(4 if len(markers) >= 4 else len(markers)):
        print(f"{cl=}, {de=}")
        ds.plot_layout(
            layout_key="RNA_UMAP",
            color_by=markers.iloc[de],
            show_fig=False,
            savename=f"data/processed/{sample}/deg_plots/cluster-{cl}_de-{de}.svg",
        )
    plt.close("all")

# table with top 5 de genes for each cluster (if present)
RNA_clusters = max(ds.cells.to_pandas_dataframe(
    columns=['RNA_leiden_cluster']).RNA_leiden_cluster)
n_deg_shown = 20
x = np.zeros((n_deg_shown, RNA_clusters), dtype=np.int)
df = pd.DataFrame(x, columns=[f'cl{i}' for i in range(1, RNA_clusters+1)])
df.index = ['de_gene_' + str(x+1) for x in df.index]
for x in range(1, RNA_clusters+1):
    deg_list = list(ds.get_markers(group_key='RNA_leiden_cluster',
                                  group_id=x)['names'])
    n_deg_genes = n_deg_shown
    if len(deg_list) < n_deg_shown:
        n_deg_genes = len(deg_list)
    if len(deg_list) > 0:
        df[f"cl{x}"] = deg_list[:n_deg_genes] + [0 for i in range(n_deg_shown - n_deg_genes)]
df
df.to_csv(f'data/processed/{sample}/top_deg_genes.csv')


# print("running paris clustering..")
# ds.run_clustering(n_clusters=leiden_clusters.nunique()[0])

# ds.run_marker_search(group_key="RNA_cluster", threshold=0.25)

# ds.plot_marker_heatmap(
#     group_key="RNA_cluster",
#     topn=5,
#     figsize=(5, 9),
#     show_fig=False,
#     savename=f"data/processed/{sample}/paris/marker_heatmap.svg",
# )

# ds.plot_layout(
#     layout_key="RNA_UMAP",
#     color_by="RNA_cluster",
#     show_fig=False,
#     savename=f"data/processed/{sample}/paris/umap.svg",
# )

# cells_in_clusters = dict(
#     sorted(dict(Counter(ds.RNA.cells.fetch("RNA_cluster"))).items())
# )
# print("cells_in_clusters for Paris clustering:")
# for key, value in cells_in_clusters.items():
#     print(f"{key}: {value}")

# for cl in range(1, leiden_clusters.nunique()[0] + 1):
#     markers = ds.get_markers(group_key="RNA_cluster", group_id=cl)["names"]

#     for hvg in range(4 if len(markers) >= 4 else len(markers)):
#         print(f"{cl=}, {hvg=}")
#         ds.plot_layout(
#             layout_key="RNA_UMAP",
#             color_by=markers.iloc[hvg],
#             show_fig=False,
#             savename=f"data/processed/{sample}/paris/hvg_plots/cluster-{cl}_hvg-{hvg}.svg",
#         )
#     plt.close("all")

print("ATAC stuff")

ds.mark_prevalent_peaks(from_assay='ATAC',
                        top_n=math.floor(ds.ATAC.feats.N*.25))

ds.make_graph(
    feat_key='prevalent_peaks',
    n_centroids=math.floor(ds.ATAC.cells.N*.25),
    lsi_skip_first=True,
    k=21,
    dims=50,
    from_assay='ATAC'
)

ds.run_umap(from_assay='ATAC', parallel=True)
ds.run_leiden_clustering(from_assay='ATAC', resolution=resolution_ATAC)

ds.plot_layout(
    layout_key='ATAC_UMAP',
    color_by='ATAC_leiden_cluster',
    cmap=cmap_ATAC,
    show_fig=False,
    savename=f'data/processed/{sample}/umap_ATAC.svg'
)

ds.plot_layout(
    layout_key='ATAC_UMAP',
    color_by='RNA_leiden_cluster',
    show_fig=False,
    cmap=cmap_RNA,
    savename=f'data/processed/{sample}/umap_ATAC-w-RNA-clusters.svg'
)
ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='ATAC_leiden_cluster',
    show_fig=False,
    cmap=cmap_ATAC,
    savename=f'data/processed/{sample}/umap_RNA-w-ATAC-clusters.svg'
)


ds.plot_layout(
    layout_key=['RNA_UMAP', 'ATAC_UMAP'],
    color_by=['ATAC_leiden_cluster', 'RNA_leiden_cluster'],
    cmap='tab20',
    width=4,
    height=4,
    n_columns=2,
    point_size=5,
    legend_onside=False,
    shuffle_df=True,
    show_fig=False,
    savename=f'data/processed/{sample}/umaps_RNA-ATAC_clusters.svg'
)

df = pd.crosstab(
    ds.cells.fetch('RNA_leiden_cluster'),
    ds.cells.fetch('ATAC_leiden_cluster')
)
df

df.to_csv(f'data/processed/{sample}/clusters_RNA-ATAC_crosstab.csv')
ds.cells.to_pandas_dataframe(columns=["ATAC_leiden_cluster"],
                             key="I").value_counts().to_csv(
                                 f'data/processed/{sample}/cells-in-clusters_ATAC.csv')
ds.cells.to_pandas_dataframe(columns=["RNA_leiden_cluster"],
                             key="I").value_counts().to_csv(
                                 f'data/processed/{sample}/cells-in-clusters_RNA.csv')


if not os.path.exists('data/external/scarf_annotations'):
    scarf.fetch_dataset(dataset_name='annotations',
                        save_path='data/external/scarf_annotations')

ds.add_melded_assay(from_assay='ATAC',
                    external_bed_fn=(
                        'data/external/scarf_annotations/annotations/mouse_GRCm38_gencode_vM25_gene_body.bed.gz'
                    ),
                    peaks_col='ids', renormalization=False,
                    assay_label='GeneScores', assay_type='RNA')

df = ds.GeneScores.feats.to_pandas_dataframe(columns=['I', 'names', 'dropOuts', 'nCells'])
df = df.loc[df.I == True].sort_values(by='nCells', ascending=False)
ds.plot_layout(
    layout_key='ATAC_UMAP', from_assay='GeneScores',
    color_by=list(df['names'])[:5],
    clip_fraction=0.01,
    n_columns=3,
    width=3,
    height=3,
    point_size=5,
    scatter_kwargs={'lw': 0.01},
    show_fig=False,
    savename=f'data/processed/{sample}/umaps_ATAC_GeneScores.svg'
)

df = ds.RNA.feats.to_pandas_dataframe(
    columns=ds.RNA.feats.columns
).sort_values(by='nCells', ascending=False)
ds.plot_layout(
    layout_key=['RNA_UMAP', 'ATAC_UMAP'],
    color_by=list(df[df.I__hvgs == True]['names'])[:2],
    from_assay='RNA',
    width=4,
    height=4,
    n_columns=2,
    point_size=5,
    show_fig=False,
    savename=f'data/processed/{sample}/umaps_RNA-ATAC_geneexp.svg'
)

plt.close("all")

print("all done!")

# ds.plot_layout(
#     layout_key='RNA_UMAP',
#     color_by='RNA_leiden_cluster',
#     show_fig=False,
#     savename=f'data/processed/{sample}/umap2.svg',
# )

# how to check PCA dimensions? e. g. elbow plot
# A: there's a elbow_plot param for .. that function

# for make_graph: if feat_key.startswith('I__'): do not prepend 'I__'

# UX for cells and feat MetaData object, looks like DataFrames, need to dig
# through docs to find out *to_dataframe method

# get KeyError for plot_cluster_tree, 9996

# add to scarf FAQ: can I export data from Seurat/scanpy etc?

# zarr_tree: intuitive to access variables by dot notation, but not possible.
# Possible to make that more obvious?

# ?ds.run_marker_search refers to scarf.markers.find_markers_by_rank, perhaps
# make that more clear?
# P.S.: should probably write documentation for that as well.
# P.P.S.: should probably write better documentation for how markers are found.

# how can the clustering assign a cluster -1??

# maybe rethink only public methods on scarf docs, don't get any documentation
# of attributes contained in __init__ methods

# docs, data org, mulit-omics
# also, in zarr trees, "synchronize the" what?
# "The sparse matrix above will load the If you would like to load a graph[..]"
# what?

# can one remove parts of the zarr tree?
# How about doing things only read mode from zarr file (with ds) and if the
# results are good, THEN save to zarr hierarchy?

# talked w para about idempotency, if same command is re-run w same params, the
# old output is replaced. This is good for workflow systems, if one works with
# marker files on knows that the results aren't changed under the hood.

# ds.plot_cells_dists doesn't automatically actually use the 'I' column, it
# always prints everything

# scarf.assay.normMethod "Performs library size normalization on the data.",
# Inappropriate name for what the function does.
