import os


configfile: "celltag-config.yaml"


#report: os.path.join(config["docs_dir"], "workflow.rst")


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


rule filter_unmapped_reads:
    input:
        start_over_marker=".smk_markers/start_over.marker",
        possorted_bam=lambda wildcards:
            f"{config['raw_dir']}/count/{config['samples_orig_names'][wildcards.sample]}/outs/gex_possorted_bam.bam",
    output:
        filtered_bam=expand(
            "{output_dir}/{sample}/celltag/filtered_and_possorted_unmapped.bam",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    shell:
        "samtools view -b -f 4 -@ 4 {input.possorted_bam} > {output.filtered_bam}"


rule filter_GFP:
    input:
        start_over_marker=".smk_markers/start_over.marker",
        possorted_bam=lambda wildcards:
            f"{config['raw_dir']}/count/{config['samples_orig_names'][wildcards.sample]}/outs/gex_possorted_bam.bam",
    output:
        filtered_bam=expand(
            "{output_dir}/{sample}/celltag/filtered_and_possorted_GFP.bam",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    shell:
        "samtools view -b -@ 4 {input.possorted_bam} GFP > {output.filtered_bam}"


rule concatenate_filtered:
    input:
        unmapped=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/filtered_and_possorted_unmapped.bam",
        gfp=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/celltag/filtered_and_possorted_GFP.bam",
    output:
        filtered_bam=expand(
            "{output_dir}/{sample}/celltag/filtered_and_possorted.bam",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
    shell:
        "cat {input.unmapped} {input.gfp} > {output.filtered_bam}"


# rule construct_filtered_bam_index:
#     input:
#         filtered_bam=lambda wildcards:
#             f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam",
#     output:
#         filtered_bam=expand(
#             "{output_dir}/{sample}/celltag/filtered_and_possorted.bam.bai",
#             output_dir=config["output_dir"],
#             allow_missing=True,
#         ),
#     shell:
#         "samtools index -@ 0 {input.filtered_bam} {output.filtered_bam}"


rule construct_filtered_bam_index:
    input:
        filtered_bam=lambda wildcards:
            f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam",
    output:
        filtered_bai=expand(
            "{output_dir}/{sample}/celltag/filtered_and_possorted.bam.bai",
            output_dir=config["output_dir"],
            allow_missing=True,
        ),
    threads:
        4  # This value - 1 will be sent to -@
    wrapper:
        "v1.3.2/bio/samtools/index"


rule unzip_10X_barcode_tsv:
    input:
        barcodes_10X=lambda wildcards:
            f"{config['raw_dir']}/count/{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        tsv=expand(
            "{output_dir}/{sample}/10X_filtered_barcodes.tsv",
            output_dir=config["interim_dir"],
            allow_missing=True,
        ),
    shell:
        "gunzip {input.barcodes_10X} --keep --to-stdout > {output.tsv}"


rule extract_celltag_reads:
    input:
        # filtered_bam=lambda wildcards:
        #     f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam",
        # filtered_bai=lambda wildcards:
        #     f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam.bai",
        possorted_bam=lambda wildcards:
            f"{config['raw_dir']}/count/{wildcards.sample}/outs/gex_possorted_bam.bam",
        barcodes_tsv=lambda wildcards:
            f"{config['interim_dir']}/{wildcards.sample}/10X_filtered_barcodes.tsv",
    output:
        celltag_reads=expand(
            "{output_dir}/{sample}/celltag/v{CT_version}.celltag.reads.out",
            output_dir=config["output_dir"],
            CT_version=config["celltag_version"],
            allow_missing=True,
        ),
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    shell:
        "samtools view {input.possorted_bam} | grep -P -f src/celltag/celltagv3-pattern > {output.celltag_reads}"


##########################################


# rule extract_celltag_reads:
#     input:
#         # filtered_bam=lambda wildcards:
#         #     f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam",
#         # filtered_bai=lambda wildcards:
#         #     f"{config['output_dir']}/{wildcards.sample}/celltag/filtered_and_possorted.bam.bai",
#         possorted_bam=lambda wildcards:
#             f"{config['raw_dir']}/count/{wildcards.sample}/outs/gex_possorted_bam.bam",
#     output:
#         celltag_reads=expand(
#             "{output_dir}/{sample}/celltag/v{CT_version}.celltag.reads.out",
#             output_dir=config["output_dir"],
#             CT_version=config["celltag_version"],
#             allow_missing=True,
#         ),
#     conda:
#         "envs/DB_Qinyu-multiome_snakemake_R.yaml"
#     shell:
#         "samtools view {input.possorted_bam} | grep -P -f src/celltagv3-pattern > {output.celltag_reads}"

# GAATTCGATGACAGGCGCAGCTTCCGAGGGATTTGAGATCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCAAAACTAGAATGCAGTGAAAAAAATGCCTTATTTGTGAAATTTGTGATGCTATTGCCTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACA


rule parse_celltag_reads:
    input:
        celltag_reads=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.reads.out",
    output:
        parsed_tsv=expand(
            "{output_dir}/{sample}/celltag/v{CT_version}.celltag.parsed.tsv",
            output_dir=config["output_dir"],
            CT_version=config["celltag_version"],
            allow_missing=True,
        ),
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    shell:
        "src/celltag/BiddyetalWorkflow/scripts/celltag.parse.reads.10x.sh -v tagregex=`cat src/celltagv3-tagregex` {input.celltag_reads} > {output.parsed_tsv}"


# rule unzip_10X_barcode_tsv:
#     input:
#         barcodes_10X=lambda wildcards: f"{config['raw_dir']}/count/{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#     output:
#         tsv=expand(
#             "{output_dir}/{sample}/10X_filtered_barcodes.tsv",
#             output_dir=config["output_dir"],
#             allow_missing=True,
#         ),
#     shell:
#         "gunzip {input.barcodes_10X} --keep --to-stdout > {output.tsv}"


rule quantify_celltags:
    input:
        celltag_reads=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.reads.out",
        parsed_tsv=lambda wildcards: f"{config['output_dir']}/{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.parsed.tsv",
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


rule clone_calling:
    input:
        matrix_rds=lambda wildcards:
            f"{config['output_dir']}/{wildcards.sample}/celltag/CT.celltag.matrix.Rds",
    output:
        clones_csv=expand("{output_dir}/{sample}/celltag/{fn}",
            output_dir=config["output_dir"],
            fn='clones.csv',
            allow_missing=True,
        ),
        clones_size_csv=expand("{output_dir}/{sample}/celltag/{fn}",
            output_dir=config["output_dir"],
            fn='clones_size.csv',
            allow_missing=True,
        ),
        jaccard_mtx=expand("{output_dir}/{sample}/celltag/d13{fn}",
            output_dir=config["output_dir"],
            fn='_Jaccard_mtx.RDS',
            allow_missing=True,
        ),
        # jaccard_plot=expand("{output_dir}/{sample}/celltag/d13{fn}",
            # output_dir=config["output_dir"],
            # fn='_Jaccard_correlation_plot.pdf',
            # allow_missing=True,
        # ),
	#params:
		#biddypath="./BiddyetalWorkflow",
    conda:
        "envs/DB_Qinyu-multiome_snakemake_R.yaml"
    script:
        "src/celltag-clonecalling.R"


# rule gather_celltag:
#     input:
#         quantified=expand(
#             "{output_dir}/{sample}/celltag/{files}",
#             output_dir=config["output_dir"],
#             files=expand(
#                 "{prefix}{filetypes}",
#                 prefix="CT",
#                 filetypes=[".celltag.stats.txt", ".matrix.tsv",
#                            ".celltag.matrix.Rds", "",]) +
#                 ['clones.csv', 'clones_size.csv', 'd13_Jaccard_mtx.RDS'],
#             sample=[config['samples'][x] for x in ['EPCR', 'Viable']]
#         ),
#     output:
#         touch('.smk_markers/celltag.done')


# rule gather_all:
#     input:
#         start_over='.smk_markers/start_over.marker',
#         celltag=expand(
#             "{output_dir}/{sample}/celltag/{files}",
#             output_dir=config["output_dir"],
#             files=expand(
#                 "{prefix}{filetypes}",
#                 prefix="CT",
#                 filetypes=[".celltag.stats.txt", ".matrix.tsv",
#                            ".celltag.matrix.Rds", "",]) +
#                 ['clones.csv', 'clones_size.csv', 'd13_Jaccard_mtx.RDS'],
#             sample=[config['samples'][x] for x in ['EPCR', 'Viable']]
#         ),
#     output:
#         touch('.smk_markers/all.done')

