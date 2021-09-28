import scarf
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
from collections import Counter


args = str(sys.argv)
valid_sample_names = ['EPCR_LSKs_RNA', 'Viable_RNA']
if sys.argv[1] not in valid_sample_names:
    sys.exit(f'"{sys.argv[1]}" is not a valid sample name. Pick one of {valid_sample_names}.')
else:
    sample = sys.argv[1]

os.makedirs(f'data/processed/{sample}/leiden/hvg_plots/', exist_ok=True)
os.makedirs(f'data/processed/{sample}/paris/hvg_plots/', exist_ok=True)

if len(sys.argv) >= 3 and sys.argv[2] == 'cache':
    pass
else:
    reader = scarf.CrH5Reader(f'data/raw/count/GFP_{sample}/outs/filtered_feature_bc_matrix.h5')
    print(reader.assayFeats)
    writer = scarf.CrToZarr(reader, zarr_fn=f'data/processed/{sample}/counts_{sample}.zarr',
                            chunk_size=(6000,2000))
    writer.dump(batch_size=12000)

ds = scarf.DataStore(zarr_loc=f'data/processed/{sample}/counts_{sample}.zarr',
                     default_assay='RNA', nthreads=16)
ds
ds.RNA.feats.head()
ds.ATAC.feats.head()
ds.auto_filter_cells(show_qc_plots=False)
ds.mark_hvgs(
        min_cells=20,
        top_n=500,
        show_plot=False
)

ds.make_graph(
        feat_key='hvgs',
        k=11,
        dims=15,
        n_centroids=100
)

ds.run_umap(
    n_epochs=250,
    spread=5,
    min_dist=1,
    parallel=True,
)

print('running leiden clustering..')

if sample == 'Viable_RNA':
    resolution=0.6
else:
    resolution=0.5

ds.run_leiden_clustering(resolution=resolution)

leiden_clusters = ds.cells.to_pandas_dataframe(
    columns=['RNA_leiden_cluster'],
    key='I'
)

ds.run_marker_search(
        group_key='RNA_leiden_cluster',
        threshold=0.25
)

ds.plot_marker_heatmap(
    group_key='RNA_leiden_cluster',
    topn=5,
    figsize=(5, 9),
    show_fig=False,
    savename=f'data/processed/{sample}/leiden/marker_heatmap.svg'
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='RNA_leiden_cluster',
    show_fig=False,
    savename=f'data/processed/{sample}/leiden/umap.svg',
)

cells_in_clusters = leiden_clusters.value_counts()
clusters = len(leiden_clusters.value_counts())
print(f'{leiden_clusters.nunique()[0]=}')
print(f'{clusters=}')
print(f'{cells_in_clusters=}')

for cl in range(1,leiden_clusters.nunique()[0]+1):
    markers = ds.get_markers(
            group_key='RNA_leiden_cluster',
            group_id=cl)['names']

    for hvg in range(4 if len(markers) >= 4 else len(markers)):
        print(f'{cl=}, {hvg=}')
        ds.plot_layout(layout_key='RNA_UMAP',
                       color_by=markers.iloc[hvg],
                       show_fig=False,
                       savename=f'data/processed/{sample}/leiden/hvg_plots/cluster-{cl}_hvg-{hvg}.svg'
        )
    plt.close('all')

print('running paris clustering..')
ds.run_clustering(n_clusters=leiden_clusters.nunique()[0])

ds.run_marker_search(
        group_key='RNA_cluster',
        threshold=0.25
)

ds.plot_marker_heatmap(
    group_key='RNA_cluster',
    topn=5,
    figsize=(5, 9),
    show_fig=False,
    savename=f'data/processed/{sample}/paris/marker_heatmap.svg'
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='RNA_cluster',
    show_fig=False,
    savename=f'data/processed/{sample}/paris/umap.svg',
)

cells_in_clusters = dict(sorted(dict(Counter(ds.RNA.cells.fetch('RNA_cluster'))).items()))
print('cells_in_clusters for Paris clustering:')
for key, value in cells_in_clusters.items():
    print(f'{key}: {value}')

for cl in range(1,leiden_clusters.nunique()[0]+1):
    markers = ds.get_markers(
            group_key='RNA_cluster',
            group_id=cl)['names']

    for hvg in range(4 if len(markers) >= 4 else len(markers)):
        print(f'{cl=}, {hvg=}')
        ds.plot_layout(layout_key='RNA_UMAP',
                       color_by=markers.iloc[hvg],
                       show_fig=False,
                       savename=f'data/processed/{sample}/paris/hvg_plots/cluster-{cl}_hvg-{hvg}.svg'
        )
    plt.close('all')

print('all done!')

# ds.plot_layout(
#     layout_key='RNA_UMAP',
#     color_by='RNA_leiden_cluster',
#     show_fig=False,
#     savename=f'data/processed/{sample}/umap2.svg',
# )

# how to check PCA dimensions? e. g. elbow plot

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
