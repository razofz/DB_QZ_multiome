import os

configfile: "config.yaml"
report: os.path.join(config['docs_dir'], 'workflow.rst')


rule all:
    input:
        'all.done'

rule start_over:
    output: touch('start_over.marker')


rule extract_celltag_reads:
    input:
        possorted_bam = lambda wildcards:
            f"{config['raw_dir']}count/{wildcards.sample}/outs/gex_possorted_bam.bam",
        start_over_marker = 'start_over.marker'
    output:
        celltag_reads=expand(
            "{output_dir}/{sample}/celltag/v{config['celltag_version']}.celltag.reads.out",
            output_dir=config['output_dir'],
            allow_missing=True),
    conda:
        'envs/DB_Qinyu-multiome_snakemake_R.yaml',
    shell:
        "samtools view {input.possorted_bam} | grep -P 'TGTACG[ACTG]{8}GAATTC' > {output.celltag_reads}"


rule parse_celltag_reads:
    input:
        celltag_reads = lambda wildcards:
            f"{config['output_dir']}{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.reads.out",
    output:
        parsed_tsv=expand(
            "{output_dir}/{sample}/celltag/v{config['celltag_version']}.celltag.parsed.tsv",
            output_dir=config['output_dir'],
            allow_missing=True),
    conda:
        'envs/DB_Qinyu-multiome_snakemake_R.yaml',
    shell:
        "./src/BiddyetalWorkflow/scripts/celltag.parse.reads.10x.sh -v tagregex='TGTACG([ACTG]{8})GAATTC' {input.celltag_reads} > {output.parsed_tsv}"


rule unzip_10X_barcode_tsv:
    input:
        barcodes_10X=f"{config['raw_dir']}count/{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        expand("{output_dir}{sample}/10X_barcodes.tsv",
               output_dir=config['output_dir'],
              allow_missing=True)
    shell:
        "gunzip {input.barcodes_10X} {output}"

rule quantify_celltags:
    input:
        celltag_reads = lambda wildcards:
            f"{config['output_dir']}{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.reads.out",
        parsed_tsv=f"{config['output_dir']}/{wildcards.sample}/celltag/v{config['celltag_version']}.celltag.parsed.tsv",
        barcodes_10X=f"{config['output_dir']}{wildcards.sample}/10X_barcodes.tsv",
    output:
        output_files=expand("{config['output_dir']}/{sample}/celltags/{prefix}{filetypes}",
              prefix='CT',
              filetypes=[
                  '.celltag.stats.txt',
                  '.matrix.tsv',
                  '.celltag.matrix.Rds'
              ],
              allow_missing=True),
    conda:
        'envs/DB_Qinyu-multiome_snakemake_R.yaml',
    params:
        prefix=f"{config['output_dir']}{wildcards.sample}/celltag/CT",
    shell:
        "Rscript BiddyetalWorkflow/scripts/matrix.count.celltags.R {input.barcodes_10X} {input.parsed_tsv} {params.prefix}"


# rule seurat_clustering:
#     input:
#         prefiltered_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_prefiltering/{wildcards.observation}/prefiltering-{wildcards.sample}.rds"
#     output:
#         seurat_object=config['output_dir'] +
#         "/seurat_clustered/{observation}/{sample}/seurat-object.rds",
#         hvgs_csv=report(config['output_dir'] +
#         "/seurat_clustered/{observation}/{sample}/hvgs-{sample}.csv",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-clustering_hvgs-csv.rst'),
#             category=steps['1-2'],
#             subcategory='{sample}'
#         ),
#         post_filtering_plot=report(config['output_dir'] +
#         "/seurat_clustered/{observation}/{sample}/post-filtering-plot.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-clustering_post-filtering-plot.rst'),
#             category=steps['1-2'],
#             subcategory='{sample}'
#         ),
#         umap_plot=report(config['output_dir'] +
#         "/seurat_clustered/{observation}/{sample}/{sample}-UMAP.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-clustering_umap-plot.rst'),
#             category=steps['1-2'],
#             subcategory='{sample}'
#         ),
#     params:
#         what='clustering'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-clustering_smk.R"
# 
# rule seurat_all_integration:
#     input:
#         seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/all/seurat-object.rds"
#     output:
#         umap_plot=report(config['output_dir'] +
#                                "/seurat_all-integration/{observation}/umap.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_unintegrated-umap.rst'),
#             category=steps['1-5'],
#             subcategory='Before integration'
#         ),
#         seurat_object=config['output_dir'] + "/seurat_all-integration/{observation}/seurat-object.rds",
#     params:
#         what='all-integration'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-all-integration_smk.R"
# 
# rule seurat_rpca_integration:
#     input:
#         resistant_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_prefiltering/{wildcards.observation}/prefiltering-resistant.rds",
#         sensitive_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_prefiltering/{wildcards.observation}/prefiltering-sensitive.rds"
#     output:
#         umap_plot=report(config['output_dir'] +
#                                "/seurat_all-integration/{observation}/RPCA-umap.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_umap.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated',
#         ),
#         umap_samples_plot=report(config['output_dir'] +
#                                "/seurat_all-integration/{observation}/RPCA-umap-samples.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_umap-samples.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated',
#         ),
#         seurat_object=config['output_dir'] +
#         "/seurat_all-integration/{observation}/RPCA-seurat-object.rds",
#     params:
#         what='rpca-integration'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-rpca-integration_smk.R"
# 
# rule seurat_integrated_diff_exp:
#     input:
#         seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_all-integration/{wildcards.observation}/{wildcards.method}-seurat-object.rds",
#         sensitive_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/sensitive/seurat-object.rds",
#         resistant_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/resistant/seurat-object.rds"
#     output:
#         # dot_plot=report(config['output_dir'] +
#                                # "/seurat_all-integration/{observation}/RPCA-umap.svg",
# # #             caption=os.path.join(config['docs_dir'],
# # #                                  'seurat-diff-exp-testing.rst'),
#             # category='Seurat, integration of all samples',
#             # subcategory='{wildcards.method}',
#         # ),
#         umap_samples_plot=report(config['output_dir'] +
#                                "/seurat_all-integration/{observation}/{method}-umap-samples-clusters.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_samples-clusters.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated',
#         ),
#         clusters_heatmap=report(config['output_dir'] +
#                                "/seurat_all-integration/{observation}/{method}-clusters-heatmap.svg",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_heatmap.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated',
#         ),
#         # conserved_clusters_csvs=expand(config['output_dir'] +
#         # "/seurat_all-integration/{observation}/{method}-{sens_cluster}-conserved.csv",
#                                        # sens_cluster=['s1', 's2', 's5', 'all'],
#                                        # allow_missing=True,
#                                       # ),
#         differential_expression_csvs=report(expand(config['output_dir'] +
#             "/seurat_all-integration/{observation}/{method}-{sens_cluster}-vs-res-diff-exp.csv",
#                                        sens_cluster=['s1', 's2', 's5', 'all'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_csv.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, CSV files (differential expression)',
#         ),
#         feature_plots_enriched_sens=report(expand(config['output_dir'] +
#             "/seurat_all-integration/{observation}/{method}-{sens_cluster}-vs-res-enriched.svg",
#                                        sens_cluster=['s1', 's2', 's5', 'all'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_featureplot-sens-enriched.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, FeaturePlots (differential expression)',
#         ),
#         feature_plots_enriched_res=report(expand(config['output_dir'] +
#             "/seurat_all-integration/{observation}/{method}-res-vs-{sens_cluster}-enriched.svg",
#                                        sens_cluster=['s1', 's2', 's5', 'all'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_featureplot-res-enriched.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, FeaturePlots (differential expression)',
#         ),
#         feature_plots_similar=report(expand(config['output_dir'] +
#             "/seurat_all-integration/{observation}/{method}-{sens_cluster}-vs-res-similar.svg",
#                                        sens_cluster=['s1', 's2', 's5', 'all'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_featureplot-similar.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, FeaturePlots (differential expression)',
#         ),
#         seurat_object=config['output_dir'] +
#         "/seurat_all-integration/{observation}/{method}-seurat-object-w-sample-clusters.rds",
#     params:
#         what='integration-diff-exp',
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-all-integrated-diff-exp_smk.R"
# 
# rule seurat_integrated_concurrent_exp:
#     input:
#         seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_all-integration/fl/RPCA-seurat-object-w-sample-clusters.rds",
#         hvgs_resistant=config['output_dir'] +
#         "/seurat_clustered/fl/resistant/hvgs-resistant.csv",
#         hvgs_sensitive=config['output_dir'] +
#         "/seurat_clustered/fl/sensitive/hvgs-sensitive.csv"
#     output:
#         hvgs_approach_csv=report(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-concurrent-expression-hvgs-approach.csv",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-hvg-csv.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#         marker_approach_csv=report(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-concurrent-expression-marker-approach.csv",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-marker-csv.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#         plots_hvg_shape=report(expand(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-similar-hvg-approach-shape-{markers}.svg",
#                                        markers=['1to4', '5to8', '9to12',
#                                                 '13to16'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-hvg-shape-plot.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#         plots_hvg_split=report(expand(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-similar-hvg-approach-split-{markers}.svg",
#                                        markers=['1to2', '3to4', '5to6', '7to8'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-hvg-split-plot.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#         plots_marker_shape=report(expand(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-similar-marker-approach-shape-{markers}.svg",
#                                        markers=['1to4', '5to8', '9to12',
#                                                 '13to16'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-marker-shape-plot.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#         plots_marker_split=report(expand(config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-similar-marker-approach-split-{markers}.svg",
#                                        markers=['1to2', '3to4', '5to6', '7to8'],
#                                        allow_missing=True,
#                                       ),
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-integration_concurrent-marker-split-plot.rst'),
#             category=steps['1-5'],
#             subcategory='Integrated, concurrent expression',
#             ),
#     params:
#         what='integration-concurrent-exp'
#     conda:
#         "../../../envs/snakes-and-pirates.yaml"
#     notebook:
#         '../../../notebooks/k562-rna-seq/integration-concurrent-exp.ipynb'
# 
# rule seurat_diff_exp_testing:
#     input:
#         seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/{wildcards.sample}/seurat-object.rds"
#     output:
#         de_genes=report(config['output_dir'] + "/seurat_diff-exp-testing/{observation}/{sample}/DE-genes.csv",
#             caption=os.path.join(config['docs_dir'],
#                                  'seurat-diff-exp-testing.rst'),
#             category=steps['1-3'],
#             subcategory='{sample}'
#         ),
#     params:
#         what='diff-exp-testing'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-diff-expr-testing_smk.R"
# 
# rule seurat_mapping:
#     input:
#         query = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/{wildcards.query}/seurat-object.rds",
#         ref = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/{wildcards.ref}/seurat-object.rds"
#     output:
#         mapped_seurat_object=config['output_dir'] + "/seurat_mapping/{observation}/{query}_on_{ref}/{reduction}/mapped_seurat-object.rds",
#         anchors_rds=config['output_dir'] + "/seurat_mapping/{observation}/{query}_on_{ref}/{reduction}/transfer_anchors.rds",
#     params:
#         what='mapping'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-mapping_smk.R"
# 
# rule export_mappings:
#     input:
#         mapped_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/mapped_seurat-object.rds"
#     output:
#         mapping_csv=config['output_dir'] + "/seurat_mapping/{observation}/{query}_on_{ref}/{reduction}/mappings_clusters.csv"
#     params:
#         what='export_mappings'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-export-mappings_smk.R"
# 
# rule relative_enrichment_csv:
#     input:
#         mappings_csv = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/mappings_clusters.csv",
#         reference_csv = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.ref}_on_{wildcards.query}/{wildcards.reduction}/mappings_clusters.csv"
#     output:
#         csv=report(config['output_dir'] +
#                    "/seurat_mapping/{observation}/{query}_on_{ref}/{reduction}/relative_enrichment.csv",
#               caption=os.path.join(config['docs_dir'],
#                                    'seurat-mapping_enrichment-csv.rst'),
#               category=steps['1-4'] + ' ({query} on {ref})',
#               subcategory='{reduction}'),
#     params:
#         what='relative_enrichment_csv'
#     container:
#         "docker://razofz/k562_scrna-seq_py:0.1"
#     script:
#         "relative_enrichment_csv_smk.py"
# 
# rule relative_enrichment_plots:
#     input:
#         mappings_csv = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/mappings_clusters.csv",
#         enrichment_csv = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/relative_enrichment.csv"
#     output:
#         enrichment_plot=report(config['report_dir'] +
#                                "/mapping/{observation}/{query}_on_{ref}/{reduction}/relative_enrichment.svg",
#               category=steps['1-4'] + ' ({query} on {ref})',
#               caption=os.path.join(config['docs_dir'],
#                                    'seurat-mapping_enrichment-plot.rst'),
#               subcategory='{reduction}'),
#         parallel_plot=report(config['report_dir'] + "/mapping/{observation}/{query}_on_{ref}/{reduction}/parallel_categories.svg",
#               caption=os.path.join(config['docs_dir'],
#                                    'seurat-mapping_parallel-plot.rst'),
#               category=steps['1-4'] + ' ({query} on {ref})',
#               subcategory='{reduction}'),
#         percent_plot=report(config['report_dir'] + "/mapping/{observation}/{query}_on_{ref}/{reduction}/percent_in_clusters.svg",
#               caption=os.path.join(config['docs_dir'],
#                                    'seurat-mapping_percent-plot.rst'),
#               category=steps['1-4'] + ' ({query} on {ref})',
#               subcategory='{reduction}'),
#     params:
#         what='relative_enrichment_plots'
#     container:
#         "docker://razofz/k562_scrna-seq_py:0.1"
#     script:
#         "relative_enrichment_plots_smk.py"
# 
# rule mapping_sanity_plots:
#     input:
#         anchorset_rds = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/transfer_anchors.rds",
#         mapped_seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_mapping/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.reduction}/mapped_seurat-object.rds"
#     output:
#         query_plot=report(config['report_dir'] + "/mapping/{observation}/{query}_on_{ref}/{reduction}/query_plot.svg",
#                caption=os.path.join(config['docs_dir'],
#                                     'seurat-mapping_query-plot.rst'),
#               category=steps['1-4'] + ' ({query} on {ref})',
#               subcategory='{reduction}'),
#         samples_plot=report(config['report_dir'] + "/mapping/{observation}/{query}_on_{ref}/{reduction}/samples_plot.svg",
#                caption=os.path.join(config['docs_dir'],
#                                     'seurat-mapping_samples-plot.rst'),
#               category=steps['1-4'] + ' ({query} on {ref})',
#               subcategory='{reduction}'),
#     params:
#         what='mapping_sanity_plots'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "mapping-sanity-plots_smk.R"
# 
# ################################################################################
# # Nabo part
# ################################################################################
# 
# rule nabo_export_seurat_metadata:
#     input:
#         seurat_object = lambda wildcards:
#             f"{config['output_dir']}/seurat_clustered/{wildcards.observation}/{wildcards.ref if wildcards.mappee == 'reference' else wildcards.query}/seurat-object.rds",
#     output:
#         metadata=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_metadata.csv",
#         hvgs=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_hvgs.csv",
#         feature_loadings=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_featureloadings.csv",
#         cell_embeddings=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_cellembeddings.csv",
#         cell_embeddings_umap=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_cellembeddings_umap.csv"
#     wildcard_constraints:
#         mappee="reference|query"
#     params:
#         what='nabo_export_seurat_metadata'
#     container:
#         "docker://razofz/k562_scrna-seq:0.1"
#     script:
#         "seurat-export-metadata_smk.R"
# 
# rule nabo_make_h5:
#     input:
#         observation_dir = lambda wildcards:
#             f"{config['observations_dirs'][wildcards.observation]}",
#         metadata = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.mappee}_metadata.csv",
#         hvgs = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.mappee}_hvgs.csv",
#         feature_loadings = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.mappee}_featureloadings.csv",
#         cell_embeddings = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/{wildcards.mappee}_cellembeddings.csv",
#     output:
#         mappee_h5=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}.h5",
#         mappee_pca_h5=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/{mappee}_hvg_pca.h5"
#     wildcard_constraints:
#         mappee="reference|query"
#     params:
#         what='nabo_make_h5'
#     container:
#         "docker://razofz/k562_scrna-seq_nabo:0.2"
#     script:
#         "nabo_make_h5.py"
# 
# rule nabo_mapping:
#     input:
#         reference_pca_h5 = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/reference_hvg_pca.h5",
#         query_pca_h5 = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/query_hvg_pca.h5",
#         query_metadata = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/query_metadata.csv",
#         reference_metadata = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/reference_metadata.csv",
#         reference_umap = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/reference_cellembeddings_umap.csv",
#     output:
#         mapping_h5=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/mapping.h5",
#         graph_gml=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/graph.gml",
#         query_mapping_scores=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/query_mapping_scores.json",
#         query_cells_clusters=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/query_cells_clusters.json",
#         query_seurat_clusters=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/query_seurat_clusters.json",
#     params:
#         what='nabo_mapping'
#     container:
#         "docker://razofz/k562_scrna-seq_nabo:0.2"
#     script:
#         "nabo_mapping_smk.py"
# 
# rule nabo_mapping_plots:
#     input:
#         query_seurat_clusters = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/query_seurat_clusters.json",
#         query_cells_clusters = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/query_cells_clusters.json",
#         query_mapping_scores = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/{wildcards.query}_on_{wildcards.ref}/query_mapping_scores.json",
#         graph_gml=config['output_dir'] + "/nabo/{observation}/{query}_on_{ref}/graph.gml",
#     output:
# #         plots=report(expand("{report_dir}/nabo/{observation}/{query}_on_{ref}/{plot}.png",
# #                      report_dir=config['report_dir'],
# #                      plot=['reference_umap',
# #                            'query_projection_all_cells',
# #                            'query_mapped_percent_in_clusters_all_cells',
# #                            'query_mapped_enrichment_in_clusters_all_cells',
# #                           ],
# #                       allow_missing=True,
# #                  ),
# #                  category='Nabo mapping ({query} on {ref})',
# #                  caption=os.path.join(config['docs_dir'], 'nabo-mapping.rst'),
# #                  ),
# #         plot_directories=report(directory(expand("{report_dir}/nabo/{observation}/{query}_on_{ref}/{dirs}",
# #                              report_dir=config['report_dir'],
# #                              dirs=['query_mapped_percent_in_clusters_per_query_cluster',
# #                                     'query_projection_per_query_cluster'],
# #                              allow_missing=True,
# #                                )),
# #                  category='Nabo mapping ({query} on {ref})',
# #                  subcategory='Per (query) cluster',
# #                  caption=os.path.join(config['docs_dir'], 'nabo-mapping.rst'),
# #                  patterns=['{cluster_id}.png'],
# #                  ),
#         reference_umap=report(config['report_dir'] + "/nabo/{observation}/{query}_on_{ref}/reference_umap.png",
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='All cells',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_reference-umap.rst'),
#                  ),
#         query_projection_all_cells=report(config['report_dir'] + "/nabo/{observation}/{query}_on_{ref}/query_projection_all_cells.png",
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='All cells',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-projection-all-cells.rst'),
#                  ),
#         query_mapped_percent_in_clusters_all_cells=report(config['report_dir'] +
#         "/nabo/{observation}/{query}_on_{ref}/query_mapped_percent_in_clusters_all_cells.png",
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='All cells',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-mapped-percent-all-cells.rst'),
#                  ),
#         query_mapped_enrichment_in_clusters_all_cells=report(config['report_dir'] +
#         "/nabo/{observation}/{query}_on_{ref}/query_mapped_enrichment_in_clusters_all_cells.png",
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='All cells',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-mapped-enrichment-all-cells.rst'),
#                  ),
#         query_projection_per_query_cluster=report(directory(config['report_dir'] +
#         "/nabo/{observation}/{query}_on_{ref}/query_projection_per_query_cluster"),
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='Per (query) cluster',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-projection-per-cluster.rst'),
#                  patterns=['{cluster_id}.png'],
#                  ),
#         query_mapped_percent_in_clusters_per_query_cluster=report(directory(config['report_dir'] +
#         "/nabo/{observation}/{query}_on_{ref}/query_mapped_percent_in_clusters_per_query_cluster"),
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='Per (query) cluster',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-mapped-percent-per-cluster.rst'),
#                  patterns=['{cluster_id}.png'],
#                  ),
#         query_mapped_enrichment_in_clusters_per_query_cluster=report(directory(config['report_dir'] +
#                                                              "/nabo/{observation}/{query}_on_{ref}/query_mapped_enrichment_in_clusters_per_query_cluster"),
#                  category=steps['2'] + ' ({query} on {ref})',
#                  subcategory='Per (query) cluster',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_query-mapped-enrichment-per-cluster.rst'),
#                  patterns=['{cluster_id}.png'],
#                  ),
# #         query_mapped_enrichment_in_clusters_per_query_cluster=directory(config['report_dir'] +
# # query_mapped_enrichment_per_cluster
# #         "/nabo/{observation}/{query}_on_{ref}/query_mapped_enrichment_in_clusters_per_query_cluster"),
#     params:
#         what='nabo_mapping_plots'
#     container:
#         "docker://razofz/k562_scrna-seq_nabo:0.2"
#     script:
#         "nabo_mapping_plots_smk.py"
# 
# rule nabo_resistant_NAs_plot:
#     input:
#         query_umap = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/resistant_on_sensitive/query_cellembeddings_umap.csv",
#         query_cells_clusters = lambda wildcards:
#             f"{config['output_dir']}/nabo/{wildcards.observation}/resistant_on_sensitive/query_cells_clusters.json",
#     output:
#         resistant_NAs_plot=report(config['report_dir'] +
#                                   "/nabo/{observation}/resistant_on_sensitive/resistant_NAs_plot.png",
#                  category=steps['2'] + ' (resistant on sensitive)',
#                  subcategory='All cells',
#                  caption=os.path.join(config['docs_dir'],
#                                       'nabo-mapping_resistant_NAs-plot.rst'),
#                  ),
#     params:
#         what='nabo_resistant_NAs_plot'
#     container:
#         "docker://razofz/k562_scrna-seq_nabo:0.2"
#     script:
#         "nabo_resistant_cells_in_enriched_sensitive_clusters_smk.py"
# 
# ################################################################################
# # Gather all final outputs, to enable topmost 'all' rule
# ################################################################################
# 
# rule gather_seurat:
#     input:
#         mapping_plots=expand("{report_dir}/mapping/{observation}/{query_on_ref}/{reduction}/{plot}.svg",
#                              report_dir=config['report_dir'],
#                              observation=config['observations'][1],
#                              query_on_ref=[f"{config['samples'][0]}_on_{config['samples'][1]}",
#                                            f"{config['samples'][1]}_on_{config['samples'][0]}"],
#                              reduction=['cca', 'pcaproject'],
#                              plot=['query_plot',
#                                    'samples_plot',
#                                    'percent_in_clusters',
#                                    'parallel_categories',
#                                    'relative_enrichment',
#                                   ]
#                          ),
#         HVGs_csvs=expand(config['output_dir'] +
#                        "/seurat_clustered/{observation}/{sample}/hvgs-{sample}.csv",
#                              observation=config['observations'][1],
#                              sample=config['samples'],
#                       ),
#         DE_csvs=expand(config['output_dir'] +
#                        "/seurat_diff-exp-testing/{observation}/{sample}/DE-genes.csv",
#                              observation=config['observations'][1],
#                              sample=config['samples'],
#                       ),
#         umap_plots=expand(config['output_dir'] +
#                           "/seurat_clustered/{observation}/{sample}/{sample}-UMAP.svg",
#                           observation=config['observations'][1],
#                           sample=config['samples'],
#                       ),
#         umap_plot=expand(config['output_dir'] + "/seurat_all-integration/{observation}/umap.svg",
#                          observation=config['observations'][1],
#                         ),
#         differential_expression_csvs=expand(config['output_dir'] +
#             "/seurat_all-integration/{observation}/{method}-{sens_cluster}-vs-res-diff-exp.csv",
#                           observation=config['observations'][1],
#                           sens_cluster=['s1', 's2', 's5', 'all'],
#                           method='RPCA',
#                           allow_missing=True,
#                       ),
#         hvgs_approach_csv=config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-concurrent-expression-hvgs-approach.csv",
#         marker_approach_csv=config['output_dir'] +
#             "/seurat_all-integration/fl/RPCA-concurrent-expression-marker-approach.csv",
#     output:
#         touch('seurat.done')
# 
# rule gather_nabo:
#     input:
# #         rules.nabo_mapping_plots.output
#         plots=expand("{report_dir}/nabo/{observation}/{query_on_ref}/{plots}.png",
#                              report_dir=config['report_dir'],
#                              observation=config['observations'][1],
#                              query_on_ref=[f"{config['samples'][0]}_on_{config['samples'][1]}",
#                                            f"{config['samples'][1]}_on_{config['samples'][0]}"],
#                              plots=['reference_umap',
#                                     'query_projection_all_cells',
#                                     'query_mapped_percent_in_clusters_all_cells',
#                                     'query_mapped_enrichment_in_clusters_all_cells',
#                                    ]
#                              ),
#         resistant_NAs_plot=expand("{report_dir}/nabo/{observation}/resistant_on_sensitive/{plot}.png",
#                              report_dir=config['report_dir'],
#                              observation=config['observations'][1],
#                              plot=['resistant_NAs_plot']
#                          ),
#         plot_directories=expand("{report_dir}/nabo/{observation}/{query_on_ref}/{dirs}",
#                              report_dir=config['report_dir'],
#                              observation=config['observations'][1],
#                              query_on_ref=[f"{config['samples'][0]}_on_{config['samples'][1]}",
#                                            f"{config['samples'][1]}_on_{config['samples'][0]}"],
#                              dirs=['query_mapped_percent_in_clusters_per_query_cluster',
#                                     'query_projection_per_query_cluster',
#                                     'query_mapped_enrichment_in_clusters_per_query_cluster',
#                                    ]
#                                ),
# #         query_projection_per_query_cluster=directory(config['report_dir'] +
# #         "/nabo/{observation}/{query}_on_{ref}/"),
# #         query_mapped_percent_in_clusters_per_query_cluster=directory(config['report_dir'] +
# #         "/nabo/{observation}/{query}_on_{ref}/query_mapped_percent_in_clusters_per_query_cluster"),
# #         query_mapped_enrichment_in_clusters_per_query_cluster=directory(config['report_dir'] +
# #         "/nabo/{observation}/{query}_on_{ref}/"),
#     output:
#         touch('nabo.done')
# 


# rule gather_all:
#     input:
#     #         nabo='nabo.done',
#     #         seurat='seurat.done',
#         'start_over.marker'
#     output:
#         touch('all.done')
