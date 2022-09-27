import os


configfile: "multiome-config.yaml"


PLOT_PATH = config["output_dir"] + "plots/"


rule all:
    input:
        config["output_dir"] + "seurat_object_with_gene_signatures.rds",
        PLOT_PATH + "is_HSC.svg",


rule seurat_preprocessing:
    input:
        h5_10X=config["cellranger_dir"] + "filtered_feature_bc_matrix.h5",
        fragments_file=config["cellranger_dir"] + "atac_fragments.tsv.gz",
    output:
        seurat_object=config["output_dir"] + "seurat_object_preprocessed.rds",
    conda:
        "envs/seurat_smk_doctored.yaml"
    script:
        "src/snakemake/seurat_preprocessing.R"


rule seurat_processing:
    input:
        seurat_object=rules.seurat_preprocessing.output.seurat_object,
    output:
        seurat_object=config["output_dir"] + "seurat_object_processed.rds",
    conda:
        "envs/seurat_smk_doctored.yaml"
    script:
        "src/snakemake/seurat_processing.R"


rule plot_integrated:
    input:
        seurat_object=rules.seurat_processing.output.seurat_object,
    output:
        plot_unintegrated_RNA_origin=PLOT_PATH + "unintegrated_RNA_origin.svg",
        plot_unintegrated_RNA_clusters=PLOT_PATH + "unintegrated_RNA_clusters.svg",
        plot_unintegrated_RNA_clusters_split=PLOT_PATH
        + "unintegrated_RNA_clusters_split.svg",
        plot_unintegrated_ATAC_origin=PLOT_PATH + "unintegrated_ATAC_origin.svg",
        plot_unintegrated_ATAC_clusters=PLOT_PATH + "unintegrated_ATAC_clusters.svg",
        plot_unintegrated_ATAC_clusters_split=PLOT_PATH
        + "unintegrated_ATAC_clusters_split.svg",
        plot_integrated_RNA_origin=PLOT_PATH + "integrated_RNA_origin.svg",
        plot_integrated_RNA_clusters=PLOT_PATH + "integrated_RNA_clusters.svg",
        plot_integrated_RNA_clusters_split=PLOT_PATH
        + "integrated_RNA_clusters_split.svg",
        plot_integrated_ATAC_origin=PLOT_PATH + "integrated_ATAC_origin.svg",
        plot_integrated_ATAC_clusters=PLOT_PATH + "integrated_ATAC_clusters.svg",
        plot_integrated_ATAC_clusters_split=PLOT_PATH
        + "integrated_ATAC_clusters_split.svg",
    conda:
        "envs/seurat_smk_doctored.yaml"
    script:
        "src/snakemake/visualization/plot_integrated.R"


rule add_gene_signatures:
    input:
        seurat_object=rules.seurat_processing.output.seurat_object,
        genesig_HSC=config["raw_dir"] + "genesig_HSCs.csv",
        genesig_earlydiff=config["raw_dir"] + "genesig_early_differentiation.csv",
    output:
        seurat_object=config["output_dir"] + "seurat_object_with_gene_signatures.rds",
    conda:
        "envs/seurat_smk_doctored.yaml"
    script:
        "src/snakemake/add_gene_signatures.R"


rule plot_gene_signatures:
    input:
        seurat_object=rules.add_gene_signatures.output.seurat_object,
    output:
        plot_is_HSC=PLOT_PATH + "is_HSC.svg",
        plot_is_earlydiff=PLOT_PATH + "is_earlydiff.svg",
    conda:
        "envs/seurat_smk_doctored.yaml"
    script:
        "src/snakemake/visualization/plot_gene_signatures.R"


# rule gather_plot_integrated:
#     input:
#         in_file=rules.plot_integrated.output,
#     output:
#         out_file=touch(PLOT_PATH + ".plot_integrated.done"),
