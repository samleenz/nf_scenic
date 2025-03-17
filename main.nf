#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* This pipeline executes the three pySCENIC steps:
*  1. GRN Inference (pySCENIC grn)
*  2. Motif Enrichment / Context (pySCENIC ctx)
*  3. Regulon Activity Scoring (pySCENIC aucell)
*
* It uses an Apptainer container specified by the 'params.container' parameter.
* Modify the input files and command options to suit your dataset and pySCENIC version.
*/


// Process for GRN Inference using pySCENIC grn
process GRN {
    container params.container
    cpus 16
    time "1h"

    input:
    path(expr)
    path(genes)
    val(seed)

    output:
    path("grn_output.tsv")

    script:
    """
    pyscenic grn \
        -o grn_output.tsv \
        --num_workers ${task.cpus} \
        --method grnboost2 \
        --seed ${seed} \
        ${expr} \
        ${genes} 
    """
}

// Process for motif enrichment / context using pySCENIC ctx
process CTX {
    container params.container
    cpus 16
    time "1h"


    input:
    path(expr)
    path(grn_file)
    path(gene_motifs)
    path(motifs)

    output:
    path("ctx_output.tsv")

    script:
    """
    pyscenic ctx \
        ${grn_file} \
        ${genes_motifs} \
        --annotations_fname ${motifs} \
        --expression_mtx_fname ${expr}\
        -o ctx_output.tsv \
        --mask_dropouts \
        --num_workers ${task.cpus}
    """
}

// Process for scoring regulon activity using pySCENIC aucell
process AUCell {
    container params.container
    cpus 12
    time "1h"
    publishDir "results/${params.project}/${seed}", mode: 'copy'

    input:
    path(expr)
    path(ctx_file)
    path(seed)

    output:
    path(output)

    script:
    output = "auc_output.tsv"
    """
    pyscenic aucell \
        ${expr} \
        ${ctx_file} \
        -o ${output} \
        --num_workers ${task.cpus}
    """
}

// Define the workflow by chaining the processes
workflow {
    grn_ch = GRN(file(params.expr), file(params.genes), params.myseed)   // Run GRN inference
    
    ctx_ch = CTX(file(params.expr), grn_ch, file(params.genes_motifs), file(params.motifs))   // Run motif enrichment using the GRN output
    
    auc_ch = AUCell(file(params.expr), ctx_ch, params.myseed) // Run regulon activity scoring using the expression matrix and CTX output
}
