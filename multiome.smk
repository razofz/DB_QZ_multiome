import os


configfile: "multiome-config.yaml"


rule seurat_preprocessing:
    input:
        h5_10X=config["raw_dir"] + "/filtered_feature_bc_matrix.h5",
        fragments_file=config["raw_dir"] + "/atac_fragments.tsv.gz",
    conda:
        "envs/seurat_smk_doctored.yaml"
    output:
        seurat_object=config["output_dir"] + "/seurat_object_preprocessed.rds",
    script:
        "src/snakemake/1.1_seurat_preprocessing.R"


rule seurat_processing:
    input:
        seurat_object=rules.seurat_preprocessing.output.seurat_object,
    conda:
        "envs/seurat_smk_doctored.yaml"
    output:
        seurat_object=config["output_dir"] + "/seurat_object_processed.rds",
    script:
        "src/snakemake/1.2_seurat_processing.R"


plot_path = config["output_dir"] + "/plots/"


rule plot_integrated:
    input:
        seurat_object=rules.seurat_processing.output.seurat_object,
    conda:
        "envs/seurat_smk_doctored.yaml"
    output:
        plot_unintegrated_RNA_origin=plot_path + "unintegrated_RNA_origin.svg",
        plot_unintegrated_RNA_clusters=plot_path + "unintegrated_RNA_clusters.svg",
        plot_unintegrated_RNA_clusters_split=plot_path + "unintegrated_RNA_clusters_split.svg",
        plot_unintegrated_ATAC_origin=plot_path + "unintegrated_ATAC_origin.svg",
        plot_unintegrated_ATAC_clusters=plot_path + "unintegrated_ATAC_clusters.svg",
        plot_unintegrated_ATAC_clusters_split=plot_path + "unintegrated_ATAC_clusters_split.svg",
        plot_integrated_RNA_origin=plot_path + "integrated_RNA_origin.svg",
        plot_integrated_RNA_clusters=plot_path + "integrated_RNA_clusters.svg",
        plot_integrated_RNA_clusters_split=plot_path + "integrated_RNA_clusters_split.svg",
        plot_integrated_ATAC_origin=plot_path + "integrated_ATAC_origin.svg",
        plot_integrated_ATAC_clusters=plot_path + "integrated_ATAC_clusters.svg",
        plot_integrated_ATAC_clusters_split=plot_path + "integrated_ATAC_clusters_split.svg",
    script:
        "src/snakemake/1.3_plot_integrated.R"


rule gather_plot_integrated:
    input:
        in_file=rules.plot_integrated.output,
    output:
        out_file=touch(plot_path + ".plot_integrated.done"),
