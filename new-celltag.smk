import os


configfile: "celltag-config.yaml"


rule all:
    input:
        ".smk_markers/all.done",


rule start_over:
    output:
        touch(".smk_markers/start_over.marker"),


rule clone_workflow:
    output:
        out_dir=directory("src/celltag/BiddyetalWorkflow")
    shell:
        "git clone https://github.com/morris-lab/BiddyetalWorkflow {output.out_dir}"


rule unzip_10X_barcode_tsv:
    input:
        barcodes_10X=lambda wildcards: f"{config['raw_dir']}/count/{config['samples_orig_names'][wildcards.sample]}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        tsv=expand(
            "{output_dir}/{sample}/10X_filtered_barcodes.tsv",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    shell:
        "gunzip {input.barcodes_10X} --keep --to-stdout > {output.tsv}"


# also installing the dependencies ‘GenomeInfoDbData’, ‘GenomeInfoDb’, ‘GenomicRanges’, ‘Biostrings’, ‘BiocGenerics’, ‘S4Vectors’, ‘IRanges’, ‘XVector’, ‘zlibbioc’, ‘BiocParallel’, ‘Rhtslib’

# mamba install bioconductor-genomeinfodb bioconductor-genomeinfodbdata bioconductor-genomicranges bioconductor-biostrings bioconductor-biocgenerics bioconductor-s4vectors bioconductor-iranges bioconductor-xvector bioconductor-zlibbioc bioconductor-biocparallel bioconductor-rhtslib

rule filter_celltags_and_prepare_collapsing:
    input:
        bam=lambda wildcards: f"{config['raw_dir']}/count/{config['samples_orig_names'][wildcards.sample]}/outs/gex_possorted_bam.bam",
        barcodes_tsv=lambda wildcards: f"{config['interim_dir']}/{wildcards.sample}/10X_filtered_barcodes.tsv",
    output:
        bam_obj=expand(
            "{output_dir}/{sample}/celltag/bam_obj_pre_collapsing.rds",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
        collapsing_file=expand(
            "{output_dir}/{sample}/celltag/collapsing/for-collapsing.txt",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    script:
        "src/celltag/1.1_create_celltag_object.R"


rule collapse:
    input:
        collapsing_file=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/collapsing/for-collapsing.txt"
    output:
        collapsing_result=expand(
            "{output_dir}/{sample}/celltag/collapsing/collapsing_result.txt",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    conda:
        "envs/starcode.yaml"
    shell:
        "starcode -s --print-clusters {input.collapsing_file} > {output.collapsing_result}"


rule collapsed_matrix_and_filtering:
    input:
        bam_obj=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/bam_obj_pre_collapsing.rds",
        collapsing_result=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/collapsing/collapsing_result.txt",
        # whitelist=f"{config['raw_dir']}/barcodes_reverse_complementary.csv",
        whitelist=f"src/celltag/BiddyetalWorkflow/whitelist/V3.CellTag.Whitelist.csv"
    output:
        bam_obj=expand(
            "{output_dir}/{sample}/celltag/bam_obj_post_collapsing.rds",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
        metric_plots_pre=expand(
            "{output_dir}/{sample}/celltag/metric_plots_pre_filtering.svg",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
        metric_plots_post_whitelist=expand(
            "{output_dir}/{sample}/celltag/metric_plots_post_whitelist.svg",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
        metric_plots_post_metric_filtering=expand(
            "{output_dir}/{sample}/celltag/metric_plots_post_metric_filtering.svg",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
    script:
        "src/celltag/1.2_.R"


rule clone_calling:
    input:
        bam_obj=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/bam_obj_post_collapsing.rds",
    output:
        bam_obj=expand(
            "{output_dir}/{sample}/celltag/bam_obj.rds",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
        clones_csv=expand(
            "{output_dir}/{sample}/celltag/clones.csv",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
    script:
        "src/celltag/1.3_.R"


