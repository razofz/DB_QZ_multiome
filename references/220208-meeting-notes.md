The total cells clusterd together with EPCR should express EPCR

DB: scarf/nabo tends to project onto edges of datasets?

Cells from fresh BM don't proliferate much

More RNA in culture cells, since they are encouraged to proliferate, so
more activity, which results in more transcripts while sequencing

look into integrating ATAC-seq

annotation of clusters from harmony data

cc on harmony-integrated

atac integration

important to ensure that the same features are measured in each datasets

quantify dataset A's peaks in dataset B to ensure common features
between the two datasets

above is if peak calling was done separately for each dataset, for our
multiome data we have aggregated cellranger run. Which means this is not
a problem for us(?). Though should probably check overlap of peaks
between EPCR and Total.

Could call peaks anew w/ macs2..
<https://satijalab.org/signac/articles/pbmc_multiomic.html>

RNA integrate (normally) on normalised data matrix (though harmony on
PCs), for atac integrate low-dimensional cell embeddings (the LSI
coordinates)

"This is much better suited to scATAC-seq data, as we typically have a
very sparse matrix with a large number of features." *(most of above is
also quotes)*

mapping approach (seurat and scarf)

seurat finds anchors, but like harmony also works on low-dimensional
embeddings (for atac: LSI)
