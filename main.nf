#!/usr/bin/env nextflow 

// read query data from 10x directory into SCE object 
QUERY_10X_DIR = Channel.fromPath(params.query_10x_dir)
process read_query_sce {
    conda "${baseDir}/envs/dropletutils.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 5
    memory { 16.GB * task.attempt }

    input:
        file(query_10x_dir) from QUERY_10X_DIR

    output:
        file("query_sce.rds") into QUERY_SCE

    """
    dropletutils-read-10x-counts.R\
        --samples ${query_10x_dir}\
        --col-names ${params.col_names}\
        --output-object-file query_sce.rds
    """
}
QUERY_SCE = QUERY_SCE.first()

// pre-process query dataset 
process preprocess_query_sce {
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }
    input:
        file(query_sce) from QUERY_SCE

    output:
        file("query_sce_processed.rds") into QUERY_SCE_PROC 

    """
    scmap-preprocess-sce.R --input-object ${query_sce} --output-sce-object query_sce_processed.rds
    """
}

projection_method = params.projection_method
QUERY_SCE_CLUST = Channel.create()
QUERY_SCE_CELL = Channel.create()

// mapping to re-direct data into correct channel 
channels = ["cluster":0, "cell":1]
QUERY_SCE_PROC.choice(QUERY_SCE_CLUST, QUERY_SCE_CELL){channels[projection_method]}

// create name-file tuples for index objects
SCMAP_CLUST_IDX = Channel.fromPath(params.clust_idx).map{ f -> tuple("${f.simpleName}", f) }
SCMAP_CELL_IDX = Channel.fromPath(params.cell_idx).map{ f -> tuple("${f.simpleName}", f) }

process run_scmap_cluster {
  conda "${baseDir}/envs/scmap.yaml"

  errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
  maxRetries 10
  memory { 16.GB * task.attempt }

  input:
    set val(acc), file(scmap_clust_idx) from SCMAP_CLUST_IDX
    file(query_sce_clust) from QUERY_SCE_CLUST.first()
   
  output:
     set val(acc), file("${acc}_cluster_proj.txt") into CLUSTER_PROJ

  """
  scmap-scmap-cluster.R\
        --index-object-file ${scmap_clust_idx}\
        --projection-object-file ${query_sce_clust}\
        --threshold ${params.pred_threshold}\
        --output-text-file  ${acc}_cluster_proj.txt\
        --output-object-file cluster_projection_result_sce.rds
  """
}

process run_scmap_cell {
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        set val(acc), file(scmap_cell_idx) from SCMAP_CELL_IDX
        file(query_sce_cell) from QUERY_SCE_CELL.first()

    output:
        set val(acc), file("${acc}_cell_proj.txt") into CELL_PROJ

    """
    scmap-scmap-cell.R\
                 --index-object-file ${scmap_cell_idx}\
                 --projection-object-file ${query_sce_cell}\
                 --cluster-col ${params.cluster_col}\
                 --output-clusters-text-file ${acc}_cell_proj.txt\
                 --closest-cells-text-file closest_cells.txt\
                 --closest-cells-similarities-text-file closest_cell_similarity.txt\
                 --output-object-file cell_projection_result_object.rds
    """
}


PROJECTION_OUTPUT = CLUSTER_PROJ.concat(CELL_PROJ)

process get_final_tables {
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        set val(acc), file(scmap_proj_table) from PROJECTION_OUTPUT

    output:
        file("${acc}_final_labs.tsv") into STD_OUTPUTS

    """
    scmap_get_std_output.R\
                --predictions-file ${scmap_proj_table}\
                --output-table ${acc}_final_labs.tsv\
                --include-scores

    """
}

process scmap_combine_labels{
    input:
        file(labels) from STD_OUTPUTS.collect()

    output:
        file("${params.label_dir}") into SCMAP_LABELS_DIR

    """
    mkdir -p ${params.label_dir}
    for file in ${labels}
    do
        mv \$file ${params.label_dir}
    done
    """
}

process select_top_labs {
    conda "${baseDir}/envs/cell_types_analysis.yaml" 
    publishDir "${params.results_dir}", mode: 'copy'

    input:
        file(labels_dir) from SCMAP_LABELS_DIR
    output:
        file("scpred_output.txt") into SCMAP_TOP_LABS

    """
    combine_tool_outputs.R\
        --input-dir ${labels_dir}\
        --scores\
        --output-table scpred_output.txt
    """

