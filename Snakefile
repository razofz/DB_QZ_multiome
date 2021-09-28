import os


configfile: "config.yaml"


report: os.path.join(config["docs_dir"], "workflow.rst")


rule all:
    input:
        ".smk_markers/all.done",


rule start_over:
    output:
        touch(".smk_markers/start_over.marker"),


rule extract_celltag_reads:
    input:
        # start_over_marker="start_over.marker",
        possorted_bam=lambda wildcards: f"{config['raw_dir']}/count/{wildcards.sample}/outs/gex_possorted_bam.bam",
    output:
        celltag_reads=expand(
            "{output_dir}/{sample}/celltag/{CT_version}.celltag.reads.out",
            output_dir=config["output_dir"],
            CT_version=config["celltag_version"],
            allow_missing=True,
        ),
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    shell:
        "samtools view {input.possorted_bam} | grep -P -f src/celltagv3-pattern > {output.celltag_reads}"


rule parse_celltag_reads:
    input:
        celltag_reads=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/{config['celltag_version']}.celltag.reads.out",
    output:
        parsed_tsv=expand(
            "{output_dir}/{sample}/celltag/{CT_version}.celltag.parsed.tsv",
            output_dir=config["output_dir"],
            CT_version=config["celltag_version"],
            allow_missing=True,
        ),
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    shell:
        "src/BiddyetalWorkflow/scripts/celltag.parse.reads.10x.sh -v tagregex=`cat src/celltagv3-tagregex` {input.celltag_reads} > {output.parsed_tsv}"


rule unzip_10X_barcode_tsv:
    input:
        barcodes_10X=lambda wildcards: f"{config['raw_dir']}/count/{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        tsv=expand(
            "{output_dir}/{sample}/10X_filtered_barcodes.tsv",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
    shell:
        "gunzip {input.barcodes_10X} --keep --to-stdout > {output.tsv}"


rule quantify_celltags:
    input:
        celltag_reads=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/{config['celltag_version']}.celltag.reads.out",
        parsed_tsv=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/{config['celltag_version']}.celltag.parsed.tsv",
        barcodes_10X=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/10X_filtered_barcodes.tsv",
    output:
        quantified=expand(
            "{output_dir}/{sample}/celltag/{prefix}{filetypes}",
            output_dir=config["output_dir"],
            prefix="CT",
            filetypes=[".celltag.stats.txt", ".matrix.tsv",
                       ".celltag.matrix.Rds", ""],
            allow_missing=True,
        ),
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    shell:
        "Rscript src/BiddyetalWorkflow/scripts/matrix.count.celltags.R {input.barcodes_10X} {input.parsed_tsv}  {output.quantified[3]} ; touch {output.quantified[3]}"


rule gather_celltag:
    input:
        quantified=expand(
            "{output_dir}/{sample}/celltag/{prefix}{filetypes}",
            output_dir=config["output_dir"],
            prefix="CT",
            filetypes=[".celltag.stats.txt", ".matrix.tsv",
                       ".celltag.matrix.Rds", ""],
            sample=[config['samples'][x] for x in ['EPCR', 'Viable']]
        ),
    output:
        touch('.smk_markers/celltag.done')


rule gather_all:
    input:
        start_over='.smk_markers/start_over.marker',
        celltag='.smk_markers/celltag.done',
    output:
        touch('.smk_markers/all.done')

