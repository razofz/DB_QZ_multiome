import os


configfile: "multiome-config.yaml"


rule seurat_processing:
    input:
        h5_10X=config["raw_dir"] + "/filtered_feature_bc_matrix.h5",
        fragments_file=config["raw_dir"] + "/atac_fragments.tsv.gz",
    conda:
        "envs/seurat_smk_doctored.yaml"
    output:
        seurat_object=config["output_dir"] + "/seurat_object.rds",
    script:
        "src/snakemake/1.1_seurat_processing.R"


rule plot_integrated:
    input:
        seurat_object=rules.seurat_processing.output.seurat_object,
    output:
        # plot_unintegrated_RNA=config["output_dir"] + "/unintegrated_RNA.svg",
        # plot_unintegrated_ATAC=config["output_dir"] + "/unintegrated_ATAC.svg",
        plot_integrated_RNA=config["output_dir"] + "/integrated_RNA.svg",
        plot_integrated_ATAC=config["output_dir"] + "/integrated_ATAC.svg",
        plot_integrated_RNA_w_clusters=config["output_dir"]
        + "/integrated_RNA_w_clusters.svg",
        plot_integrated_ATAC_w_clusters=config["output_dir"]
        + "/integrated_ATAC_w_clusters.svg",
    script:
        "src/snakemake/1.2_plot_integrated.R"
